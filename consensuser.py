#!/usr/bin/env python
# coding=utf-8
"""
If running as a script, run with the -h flag for usage documentation.
"""

import argparse
import copy
import logging
import multiprocessing
import os
import re
import statistics
import sys
import time
from glob import glob

import Bio
import numpy as np
import pysam
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cutadapt.adapters import LinkedAdapter
from cutadapt.seqio import Sequence


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


# modified summary_align.gap_consensus() ###############################################################################

# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

# https://github.com/biopython/biopython/blob/master/LICENSE.rst
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

    log.debug("Creating new incremented directory - {}".format(new_path))
    os.makedirs(new_path)
    return new_path


def read_bam(fn, bc_whitelist=None, return_dict=False):
    assert os.path.isfile(fn)

    samobj = pysam.AlignmentFile(fn, "rb", check_sq=False)
    rec_itr = samobj.fetch(until_eof=True)

    for rec in rec_itr:
        bc = str(sorted([int(i) for i in rec.get_tag("bc").tolist()]))
        if bc_whitelist and bc not in bc_whitelist:
            continue

        rec_id = rec.query_name
        length = rec.query_length
        rq = rec.get_tag("rq")
        np = rec.get_tag("np")

        if return_dict:
            yield {"rec_id": rec_id,
                   "length": length,
                   "np": np,
                   "rq": rq,
                   "bc": bc,
                   "sequence": rec.query_sequence,
                   "quality": quality_string_map(rec.query_qualities.tolist())}

        else:
            descr = "{} len:{} np:{} rq:{} bc:{}".format(rec_id, length, np, rq, bc)
            record = SeqRecord(Seq(rec.query_sequence, single_letter_alphabet),
                               id=rec_id, name=rec_id, description=descr)
            record.letter_annotations["phred_quality"] = rec.query_qualities.tolist()
            yield record


def quality_string_map(qualities, qual_str_to_list=False):
    """By default, this function takes a list of quality score integers
    and returns a phred quality string. If qual_str_to_list is set to True,
    it will convert in the oppisite direction.
    """
    # precompute dictionaries on first run
    if ("int_to_char" not in quality_string_map.__dict__
            and "char_to_int" not in quality_string_map.__dict__):
        offset = ord("!")  # 33, ! will have q of 0
        quality_string_map.int_to_char = dict()
        quality_string_map.char_to_int = dict()

        for q in range(0, 94):
            char = chr(q + offset)
            quality_string_map.int_to_char[q] = char
            quality_string_map.char_to_int[char] = q

    if qual_str_to_list:
        char_to_int = quality_string_map.char_to_int
        return [char_to_int[q] for q in qualities]

    # qual_list_to_str (default)
    int_to_char = quality_string_map.int_to_char
    return "".join([int_to_char[q] for q in qualities])


def trim_both_ends(seq_rec, primer_a, primer_b, reverse_complement=False):
    # primer_b = str(Seq(primer_b).reverse_complement())
    if reverse_complement:
        rc = seq_rec.reverse_complement()
        un_trimed_seq = str(rc.seq)
        un_trimed_qual = quality_string_map(rc.letter_annotations["phred_quality"])
    else:
        un_trimed_seq = str(seq_rec.seq)
        un_trimed_qual = quality_string_map(seq_rec.letter_annotations["phred_quality"])

    adapter = LinkedAdapter(primer_a, primer_b,
                            front_restriction=None,
                            back_restriction=None,
                            require_both=True)
    sequence = Sequence(seq_rec.id, un_trimed_seq, un_trimed_qual)
    match = adapter.match_to(sequence)

    if match:
        trimed_seq_rec = copy.deepcopy(seq_rec)
        del trimed_seq_rec.letter_annotations["phred_quality"]
        trimed_seq_rec.seq = Seq(match.trimmed().sequence, alphabet=IUPAC.ambiguous_dna)
        trimed_seq_rec.letter_annotations["phred_quality"] = quality_string_map(match.trimmed().qualities,
                                                                                qual_str_to_list=True)
        return trimed_seq_rec
    else:
        return None


def mask_seq_record(seq_rec, min_score, mask_char="N", inplace=False):
    if not inplace:
        masked_seq_rec = copy.deepcopy(seq_rec)
    else:
        masked_seq_rec = seq_rec

    base_list = list(masked_seq_rec.seq)
    for loc in range(len(base_list)):
        if masked_seq_rec.letter_annotations["phred_quality"][loc] < min_score:
            base_list[loc] = mask_char

    masked_seq_rec.seq = Seq("".join(base_list), alphabet=IUPAC.ambiguous_dna)

    return masked_seq_rec


def trim_and_mask_seq_records(records, primer_a, primer_b, min_base_score, basename, mask_char="N",
                              min_seq_score=None):
    for seq_rec in records:
        trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, reverse_complement=False)

        # primers found in forword direction
        if trimed_seq_rec is not None:
            qualities = trimed_seq_rec.letter_annotations["phred_quality"]
            if len(qualities) == 0:  # todo: only check if we care
                log.info("empty qualities list, excluding - {} {}".format(basename, seq_rec.id))
                continue
            avg_score = np.mean(qualities)
            if min_seq_score and (avg_score < min_seq_score):
                log.info("seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                              basename, seq_rec.id))
                continue
            yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

        # primers not found in forword direction
        else:
            trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, reverse_complement=True)
            if trimed_seq_rec is not None:  # primers found in reverse direction
                qualities = trimed_seq_rec.letter_annotations["phred_quality"]
                if len(qualities) == 0:
                    log.info("empty qualities list, excluding - {} {}".format(basename, seq_rec.id))
                    continue
                avg_score = np.mean(qualities)
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
    time.sleep(15)  # todo: this should be optional and a user set wait time

    if args.in_file_list is None:
        in_file_list = [args.in_file]
    else:
        with open(args.in_file_list, "r") as f:
            lines = f.readlines()
            in_file_list = []
            for line in lines:
                in_file = line.strip()
                if os.path.isfile(in_file):
                    in_file_list.append(in_file)
                else:
                    log.info("file not found - {}".format(in_file))

    for fn in in_file_list:
        yield {
            # files/directories
            "input_fn": fn,
            "output_dir": output_dir,
            # settings/options
            "aligner": args.aligner,
            "mask_char": args.mask_char,
            # base filters
            "min_base_score": args.min_base_score,
            # sequence filters
            "sequence_max_mask": args.sequence_max_mask,
            "max_len_delta": args.max_len_delta,
            # optional sequence filters
            "expected_length": args.expected_length,
            "max_len": args.max_len,
            "min_seq_score": args.min_seq_score,
            # alignment filters
            "min_seq_count": args.min_seq_count,
            "alignment_max_amb": args.alignment_max_amb,
            # consensus options
            "consensus_threshold": args.consensus_threshold,
            "consensus_require_multiple": args.consensus_require_multiple,
            "consensus_ambiguous_char": args.consensus_ambiguous_char,
            "consensus_ignore_mask_char": args.consensus_ignore_mask_char}


def spawn(param_dict):
    process_fastq(**param_dict)


def process_fastq(input_fn,
                  output_dir,
                  aligner="muscle",
                  mask_char="N",
                  min_base_score=60,
                  sequence_max_mask=None,
                  max_len_delta=None,
                  expected_length=None,
                  max_len=None,
                  min_seq_score=None,
                  min_seq_count=1,
                  alignment_max_amb=None,
                  consensus_threshold=0.7,
                  consensus_require_multiple=False,
                  consensus_ambiguous_char="X",
                  consensus_ignore_mask_char=False):

    # parse filename
    basename = os.path.splitext(os.path.basename(input_fn))[0]
    serch_obj = re.search(r"_([^._]+?)\.([^.]+?)\.([^.]+?)\.([^.]+?)$", basename)  # todo: regex should be user set
    if not serch_obj:
        log.info("couldn't parse adapters from filename, skipping - {}".format(input_fn))
        return
    else:
        primer_a = serch_obj.group(3)
        primer_b = serch_obj.group(4)

    # parse fastq file
    if max_len is None:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna))
    else:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna) if len(r) < max_len)

    clean_records = list(trim_and_mask_seq_records(records, primer_a, primer_b,
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

    # files/directories
    parser.add_argument('-i', '--in_file',
                        default=('all.ccs.3347.COI_ATAGCGACGCGATATA'
                                 '.AGCGTCTCGCATCATG'
                                 '.TYTCAACDAAYCAYAAAGATATTGA'
                                 '.TAATATGGCAGATTAGTGCAATGGA'
                                 '.fastq'),
                        help='path to input file')

    parser.add_argument('-o', '--out_dir', default="output",
                        help='path to output dir')

    # batch mode filelist settings
    parser.add_argument('--in_file_list', default=None,
                        help='path to text file w/ an in_file filename on each line')

    # settings/options
    parser.add_argument('-l', '--aligner', default="muscle",
                        help='the alignment software to use. {clustalw, muscle}')

    parser.add_argument('--mask_char', default="N",
                        help='char used to mask bases with quality below threshold')

    # base filters
    parser.add_argument('-n', '--min_base_score', type=int, default="60",
                        help='score below which a nuc is masked')

    # sequence filters
    parser.add_argument('-p', '--sequence_max_mask', type=int, default="5",
                        help='number of mask_char allowed in sequences to be aligned')

    parser.add_argument('-d', '--max_len_delta', type=int, default="5",
                        help='allowed variation from mode of sequence length')

    # optional sequence filters
    parser.add_argument('--expected_length', type=int, default=None,
                        help='optional, replaces average in max_len_delta filter')

    parser.add_argument('--max_len', type=int, default=None,
                        help='max length overwhich a seq is excluded')

    parser.add_argument('--min_seq_score', type=int, default=None,
                        help='avg score below which a seq is excluded')

    # alignment filters
    parser.add_argument('-s', '--min_seq_count', type=int, default="5",
                        help='min count of sequences before alignment excluded')

    parser.add_argument('-t', '--alignment_max_amb', type=int, default="5",
                        help='number of consensus_ambiguous_char allowed in final consensus sequences')

    # consensus options
    parser.add_argument('--consensus_threshold', type=float, default=".7",
                        help='proportion threshold for consensus to call a base.')

    parser.add_argument('--consensus_require_multiple', action='store_true',
                        help='require multiple instanses for consensus to call a base')

    parser.add_argument('--consensus_ambiguous_char', default="n",
                        help='char representing bases under consensus threshold')

    parser.add_argument('--consensus_ignore_mask_char', action='store_true',
                        help='discount mask char when calculating consensus')

    args = parser.parse_args()

    if args.in_file_list:
        processor_count = min(len(args.in_file_list), multiprocessing.cpu_count())
    else:
        processor_count = multiprocessing.cpu_count()

    with multiprocessing.Pool(processor_count) as p:
        data = p.map(spawn, param_dict_generator(args))

    print(data)


if __name__ == "__main__":
    # run main
    try:
        exit(main())
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
