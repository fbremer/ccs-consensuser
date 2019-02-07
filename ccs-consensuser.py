#!/usr/bin/env python
# coding=utf-8
"""
If running as a script, run with the -h flag for usage documentation.
"""

import argparse
import copy
import logging
import multiprocessing as mp
import os
import re
import statistics
import sys
from glob import glob

import Bio
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# logging ##############################################################################################################

class OneLineExceptionFormatter(logging.Formatter):
    """custom logging configuration for exceptions to
    format on a single line, but with <newline> markers
    for easy reformatting in a text editor"""

    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        """Format error message to fit on a single line.
        """
        result = super().format(record)
        if record.exc_text:
            result = result.replace(r"\n", "<newline>")
        return result


handler = logging.StreamHandler()
formatter = OneLineExceptionFormatter(logging.BASIC_FORMAT)
handler.setFormatter(formatter)
root = logging.getLogger()
root.addHandler(handler)

root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
log = logging.getLogger("ccs-consensuser.py")


# external helper functions used with permission #######################################################################
# from: https://github.com/chapmanb/bcbb/blob/master/align/adaptor_trim.py under the MIT license

def _remove_adaptor(seq, region, right_side=True):
    """This function adapted from https://github.com/chapmanb/bcbb/blob/master/align/adaptor_trim.py
    Remove an adaptor region and all sequence to the right or left.
    """
    if right_side:
        try:
            pos = seq.find(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.find(region)
        return seq[:pos]
    else:
        try:
            pos = seq.rfind(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.rfind(region)
        return seq[pos + len(region):]


def trim_adaptor(seq, adaptor, primer_mismatch, right_side=True):
    """Trim the given adaptor sequence from a starting sequence.
    * seq can be either of:
       - string
       - Seq
    * adaptor is a string sequence
    * primer_mismatch specifies how many errors are allowed in the match between
    adaptor and the base sequence. Matches with more than this number of errors
    are not allowed.
    """
    gap_char = '-'
    exact_pos = str(seq).find(adaptor)
    if exact_pos >= 0:
        seq_region = str(seq[exact_pos:exact_pos + len(adaptor)])
        adapt_region = adaptor
    else:
        aligns = pairwise2.align.localms(str(seq), str(adaptor),
                                         5.0, -4.0, -9.0, -0.5, one_alignment_only=True,
                                         gap_char=gap_char)
        if len(aligns) == 0:
            adapt_region, seq_region = ("", "")
        else:
            seq_a, adaptor_a, score, start, end = aligns[0]
            adapt_region = adaptor_a[start:end]
            seq_region = seq_a[start:end]
    matches = sum((1 if s == adapt_region[i] else 0) for i, s in
                  enumerate(seq_region))
    # too many errors -- no trimming
    if (len(adaptor) - matches) > primer_mismatch:
        return seq
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq, seq_region.replace(gap_char, ""),
                               right_side)


def trim_adaptor_w_qual(seq, qual, adaptor, primer_mismatch, right_side=True):
    """Trim an adaptor with an associated quality string.
    Works like trimmed adaptor, but also trims an associated quality score.
    """
    assert len(seq) == len(qual)
    tseq = trim_adaptor(seq, adaptor, primer_mismatch, right_side=right_side)
    if right_side:
        pos = seq.find(tseq)
    else:
        pos = seq.rfind(tseq)
    tqual = qual[pos:pos + len(tseq)]
    assert len(tseq) == len(tqual)
    return tseq, tqual


# modified summary_align.gap_consensus() ###############################################################################
"""Modified such that if consensus_ignore_mask_char flag is set, the sequence with mask_char is ignored in residue 
calculations at the masked position
"""


def gap_consensus(summary_align, threshold=.7, mask_char="N", consensus_ambiguous_char="X",
                  consensus_alpha=None, require_multiple=False, consensus_ignore_mask_char=False):
    """Output a fast consensus sequence of the alignment, allowing gaps.

    Same as dumb_consensus(), but allows gap on the output.

    Things to do:
     - Let the user define that with only one gap, the result
       character in consensus is gap.
     - Let the user select gap character, now
       it takes the same as input.

    """
    # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
    consensus = ''

    # find the length of the consensus we are creating
    con_len = summary_align.alignment.get_alignment_length()

    # go through each seq item
    for n in range(con_len):
        # keep track of the counts of the different atoms we get
        atom_dict = {}
        num_atoms = 0

        for record in summary_align.alignment:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            if n < len(record.seq):
                if consensus_ignore_mask_char and (record.seq[n] == mask_char):
                    continue
                if record.seq[n] not in atom_dict:
                    atom_dict[record.seq[n]] = 1
                else:
                    atom_dict[record.seq[n]] += 1

                num_atoms += 1

        max_atoms = []
        max_size = 0

        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)

        if require_multiple and num_atoms == 1:
            consensus += consensus_ambiguous_char
        elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
            consensus += max_atoms[0]
        else:
            consensus += consensus_ambiguous_char

    # we need to guess a consensus alphabet if one isn't specified
    if consensus_alpha is None:
        # noinspection PyProtectedMember
        consensus_alpha = summary_align._guess_consensus_alphabet(consensus_ambiguous_char)

    return Seq(consensus, consensus_alpha)


# helper functions #####################################################################################################

def get_unique_dir(path, width=3):
    # if it doesn't exist, create
    if not os.path.isdir(path):
        log.debug("Creating new directory - {}".format(path))
        os.makedirs(path)
        return path

    # if it's empty, use
    if not os.listdir(path):
        log.debug("Using empty directory - {}".format(path))
        return path

    # otherwise, increment the highest number folder in the series

    def get_trailing_number(search_text):
        serch_obj = re.search(r"([0-9]+)$", search_text)
        if not serch_obj:
            return 0
        else:
            return int(serch_obj.group(1))

    dirs = glob(path + "*")
    next_num = sorted([get_trailing_number(d) for d in dirs])[-1] + 1
    new_path = "{0}_{1:0>{2}}".format(path, next_num, width)
    os.makedirs(new_path)

    return new_path

def trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=False):
    if reverse_complement:
        rc = seq_rec.reverse_complement()
        un_trimed_seq = rc.seq
        un_trimed_qual = rc.letter_annotations["phred_quality"]
    else:
        un_trimed_seq = seq_rec.seq
        un_trimed_qual = seq_rec.letter_annotations["phred_quality"]

    # primer A,B found
    found_a = False
    found_b = False

    half_trimed_seq, half_trimed_qual = trim_adaptor_w_qual(un_trimed_seq,
                                                            un_trimed_qual,
                                                            adaptor=primer_a,
                                                            primer_mismatch=primer_mismatch,
                                                            right_side=False)
    if len(half_trimed_seq) < len(un_trimed_seq):
        found_a = True

    full_trimed_seq, full_trimed_qual = trim_adaptor_w_qual(half_trimed_seq,
                                                            half_trimed_qual,
                                                            adaptor=primer_b,
                                                            primer_mismatch=primer_mismatch,
                                                            right_side=True)
    if len(full_trimed_seq) < len(half_trimed_seq):
        found_b = True

    if found_a and found_b:
        trimed_seq_rec = copy.deepcopy(seq_rec)
        del trimed_seq_rec.letter_annotations["phred_quality"]
        trimed_seq_rec.seq = full_trimed_seq
        trimed_seq_rec.letter_annotations["phred_quality"] = full_trimed_qual
        return trimed_seq_rec
    else:
        return None


def mask_seq_record(seq_rec, min_score, mask_char="N", inplace=False):
    if not inplace:
        masked_seq_req = copy.deepcopy(seq_rec)
    else:
        masked_seq_req = seq_rec

    base_list = list(masked_seq_req.seq)
    for loc in range(len(base_list)):
        if masked_seq_req.letter_annotations["phred_quality"][loc] < min_score:
            base_list[loc] = mask_char

    masked_seq_req.seq = Seq("".join(base_list), alphabet=IUPAC.ambiguous_dna)

    return masked_seq_req


def trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch, min_base_score, basename, mask_char="N",
                              min_seq_score=None):
    for seq_rec in records:
        trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=False)

        # primers found in forword direction
        if trimed_seq_rec is not None:
            avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
            if min_seq_score and (avg_score < min_seq_score):
                log.info("seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                              basename, seq_rec.id))
                continue
            yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

        # primers not found in forword direction
        else:
            trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=True)
            if trimed_seq_rec is not None:  # primers found in reverse direction
                avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
                if min_seq_score and (avg_score < min_seq_score):
                    log.info(
                        "seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                             basename, seq_rec.id))
                    continue
                yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

            # primers not found in either direction
            else:
                log.info("seq excluded - primers not found - {} {}".format(basename, seq_rec.id))
                continue


# core logic ###########################################################################################################

def param_dict_generator(args):
    output_dir = get_unique_dir(args.out_dir)

    if args.in_file_list is None:
        in_file_list = [args.in_file]
    else:
        with open(args.in_file_list, "r") as f:
            in_file_list = [os.path.join(args.in_dir, l.strip()) for l in f.readlines()]

    for fn in in_file_list:
        yield {"input_fn": fn,
               "output_dir": output_dir,
               "primer_mismatch": args.primer_mismatch,
               "min_base_score": args.min_base_score,
               "min_seq_score": args.min_seq_score,
               "min_seq_count": args.min_seq_count,
               "max_len": args.max_len,
               "aligner": args.aligner,
               "sequence_max_mask": args.sequence_max_mask,
               "alignment_max_amb": args.alignment_max_amb,
               "max_len_delta": args.max_len_delta,
               "expected_length": args.expected_length,
               "consensus_threshold": args.consensus_threshold,
               "consensus_require_multiple": args.consensus_require_multiple,
               "mask_char": args.mask_char,
               "consensus_ambiguous_char": args.consensus_ambiguous_char,
               "consensus_ignore_mask_char": args.consensus_ignore_mask_char}


def spawn(param_dict):
    process_fastq(**param_dict)


def process_fastq(input_fn, output_dir, primer_mismatch, min_base_score, min_seq_score=None, min_seq_count=1,
                  max_len=None, aligner="muscle", sequence_max_mask=None, alignment_max_amb=None, max_len_delta=None,
                  expected_length=None, consensus_threshold=0.7, consensus_require_multiple=False,
                  mask_char="N", consensus_ambiguous_char="X", consensus_ignore_mask_char=False):
    # parse filename
    basename = os.path.splitext(os.path.basename(input_fn))[0]
    # noinspection PyTypeChecker
    primer_a = basename.split("_")[-1].split(".")[2]
    # noinspection PyTypeChecker
    primer_b = basename.split("_")[-1].split(".")[3]

    # parse fastq file
    if max_len is None:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna))
    else:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna) if len(r) < max_len)

    clean_records = list(trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch,
                                                   min_base_score, basename, mask_char, min_seq_score))

    if len(clean_records) < min_seq_count:
        log.info("alignment excluded - seq_count:{} < min_seq_count:{} after trim_and_mask_seq_records - {}".format(
            len(clean_records), min_seq_count, basename))
        return

    # define filter functions
    def mask_count_filter(r, _sequence_max_mask, _basename):
        _n_count = r.seq.upper().count(mask_char)
        if _n_count > _sequence_max_mask:
            log.info("seq excluded - mask_count:{} > sequence_max_mask:{} - {} {}".format(_n_count, _sequence_max_mask,
                                                                                          _basename, r.id))
            return False
        else:
            return True

    def len_variance_filter(r, _typical_len, _max_len_delta, _basename):
        len_delta = abs(len(r.seq) - _typical_len)
        if len_delta > _max_len_delta:
            log.info("seq excluded - len_delta:{} > max_len_delta:{} - {} {}".format(len_delta, _max_len_delta,
                                                                                     _basename, r.id))
            return False
        else:
            return True

    # apply mask_count_filter
    if sequence_max_mask is not None:
        clean_records = [r for r in clean_records if mask_count_filter(r, sequence_max_mask, basename)]

    if len(clean_records) < min_seq_count:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_count:{} after mask_count_filter - {}".format(
                len(clean_records),
                min_seq_count, basename))
        return

    # determin typical_len for len_variance_filter
    # should be last per seq filter b/c requires stats on full seq list
    if max_len_delta is not None:
        if expected_length is not None:
            typical_len = expected_length
        else:
            try:
                typical_len = statistics.mode([len(r) for r in clean_records])
            except statistics.StatisticsError as _e:
                log.info("failover from mode to mean - {} - {}".format(_e, basename))
                try:
                    typical_len = statistics.mean([len(r) for r in clean_records])
                except statistics.StatisticsError as __e:
                    log.info("alignment excluded - {} {} - {}".format(_e, __e, basename))
                    return
        # apply len variance filter
        clean_records = [r for r in clean_records
                         if len_variance_filter(r, typical_len, max_len_delta, basename)]

    if len(clean_records) < min_seq_count:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_count:{} after len_variance_filter - {}".format(
                len(clean_records),
                min_seq_count, basename))
        return

    # write clean fasta files
    with open(os.path.join(output_dir, basename) + ".fasta", "wt") as f:
        if len(clean_records) == 1:
            r = clean_records[0]
            mask_count = r.seq.upper().count(mask_char)
            r.description = "seq_count:1 mask_count:{} seq_len:{}".format(mask_count, len(r.seq))
            pass
        else:
            for r in clean_records:
                f.write(r.format("fasta"))

    # align clean fasta
    if aligner == "muscle":
        cline = MuscleCommandline(input=os.path.join(output_dir, basename) + ".fasta",
                                  out=os.path.join(output_dir, basename) + ".aln", clw=True)
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as _e:
            log.info("alignment failed - {} - {}".format(_e, basename))
            return

    elif aligner == "clustalw":
        cline = ClustalwCommandline("clustalw2", infile=os.path.join(output_dir, basename) + ".fasta")
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as _e:
            log.info("alignment failed - {} - {}".format(_e, basename))
            return

    else:
        log.info("alignment failed - invalid aligner: {} - {}".format(aligner, basename))
        return

    # parse aligned fasta
    with open(os.path.join(output_dir, basename) + ".aln", "r") as f:
        alignment = AlignIO.read(f, "clustal", alphabet=IUPAC.ambiguous_dna)

    # take consensus of aligned fasta
    summary_align = AlignInfo.SummaryInfo(alignment)

    # write consensus fasta
    consensus = gap_consensus(summary_align, threshold=consensus_threshold, mask_char=mask_char,
                              consensus_ambiguous_char=consensus_ambiguous_char, consensus_alpha=IUPAC.ambiguous_dna,
                              require_multiple=consensus_require_multiple,
                              consensus_ignore_mask_char=consensus_ignore_mask_char)
    amb_count = consensus.upper().count(mask_char)
    if (alignment_max_amb is not None) and (amb_count > alignment_max_amb):
        log.info(
            "alignment excluded - amb_count:{} > alignment_max_amb:{} - {}".format(amb_count, alignment_max_amb,
                                                                                   basename))
        return
    seq = Seq(data=str(consensus), alphabet=IUPAC.ambiguous_dna)
    description = "seq_count:{} amb_count:{} seq:len:{}".format(len(alignment), amb_count, len(consensus))
    # noinspection PyTypeChecker
    seq_rec = SeqRecord(seq=seq, id=basename.split("prime_")[0] + "prime", description=description)
    with open(os.path.join(output_dir, basename) + ".{}.consensus.fasta".format(len(alignment)), "wt") as f:
        fasta_entry = seq_rec.format("fasta").strip().split("\n")
        fasta_entry = fasta_entry[0] + "\n" + "".join(fasta_entry[1:]) + "\n"
        f.write(fasta_entry)


# cli and multiprocessing ##############################################################################################

class HelpAndQuitOnFailParser(argparse.ArgumentParser):
    """custom argparse configuration
    if error parsing, prints help and exits"""

    def error(self, message):
        sys.stderr.write('error: {}\n'.format(message))
        self.print_help()
        sys.exit(2)


def main():
    parser = HelpAndQuitOnFailParser()

    # settings/options
    parser.add_argument('-l', '--aligner', help='the alignment software to use. {clustalw, muscle}', default="muscle")
    parser.add_argument('--mask_char', help='char used to mask bases with quality below threshold', default="N")

    # files/directories
    parser.add_argument('-i', '--in_file', help='path to input file',
                        default=('Final.HQpolish99nomaxmin600.ccs.BX140414_001.5prime_'
                                 'CACAGAGACACGCACA.'
                                 'CGTCTCTATCTCTCTA.'
                                 'GCAGTCGAACATGTAGCTGACTCAGGTCACTCGCCTAAACTTCAGCCATT.'
                                 'TGATTYTTTGGACACCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCA.'
                                 'fastq'))
    parser.add_argument('-o', '--out_dir', help='path to output dir', default="output")

    # batch mode filelist settings
    parser.add_argument('--in_file_list', help='path to text file w/ an in_file filename on each line', default="")
    parser.add_argument('--in_dir', help='input dir of files named in in_file_list', default=None)

    # base filters
    parser.add_argument('-n', '--min_base_score', type=int, help='score below which a nuc is masked', default="60")

    # sequence filters
    parser.add_argument('-e', '--primer_mismatch', type=int, help='number of errors allowed in primer match',
                        default="2")
    parser.add_argument('-p', '--sequence_max_mask', type=int,
                        help='number of mask_char allowed in sequences to be aligned', default="5")
    parser.add_argument('-d', '--max_len_delta', type=int, help='allowed variation from mode of sequence length',
                        default="5")
    # optional sequence filters
    parser.add_argument('--expected_length', type=int, help='optional, replaces average in max_len_delta filter',
                        default=None)
    parser.add_argument('--max_len', type=int, help='max length overwhich a seq is excluded', default=None)
    parser.add_argument('--min_seq_score', type=int, help='avg score below which a seq is excluded', default=None)

    # alignment filters
    parser.add_argument('-s', '--min_seq_count', type=int, help='min count of sequences before alignment excluded',
                        default="5")
    parser.add_argument('-t', '--alignment_max_amb', type=int,
                        help='number of consensus_ambiguous_char allowed in final consensus sequences', default="5")

    # consensus options
    parser.add_argument('--consensus_threshold', type=float, help='proportion threshold for consensus to call a base.',
                        default=".7")
    parser.add_argument('--consensus_require_multiple', help='require multiple instanses for consensus to call a base',
                        action='store_true')
    parser.add_argument('--consensus_ambiguous_char', help='char representing bases under consensus threshold',
                        default="n")
    parser.add_argument('--consensus_ignore_mask_char', help='discount mask char when calculating consensus',
                        action='store_true')

    args = parser.parse_args()

    processor_count = min(len(args.sample_names), mp.cpu_count())
    with mp.Pool(processor_count) as p:
        data = p.map(spawn, param_dict_generator(args))

    print(data)


if __name__ == "__main__":
    # run main
    try:
        exit(main())
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
