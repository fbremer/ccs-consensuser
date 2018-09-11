#!/usr/bin/env python
# coding=utf-8
"""
If running as a script, run with the -h flag for usage documentation.
"""

import argparse
import copy
import itertools
import logging
import multiprocessing as mp
import os
import sys

import Bio
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# logging configuration
# noinspection PyMissingOrEmptyDocstring
class OneLineExceptionFormatter(logging.Formatter):
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
root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
root.addHandler(handler)

log = logging.getLogger("ccs-consensuser.py")


# some functions from: https://github.com/chapmanb/bcbb/blob/master/align/adaptor_trim.py
# _remove_adaptor
# trim_adaptor
# trim_adaptor_w_qual

def _remove_adaptor(seq, region, right_side=True):
    """Remove an adaptor region and all sequence to the right or left.
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


def trim_adaptor(seq, adaptor, num_errors, right_side=True):
    """Trim the given adaptor sequence from a starting sequence.
    * seq can be either of:
       - string
       - Seq
    * adaptor is a string sequence
    * num_errors specifies how many errors are allowed in the match between
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
    if (len(adaptor) - matches) > num_errors:
        return seq
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq, seq_region.replace(gap_char, ""),
                               right_side)


def trim_adaptor_w_qual(seq, qual, adaptor, num_errors, right_side=True):
    """Trim an adaptor with an associated quality string.
    Works like trimmed adaptor, but also trims an associated quality score.
    """
    assert len(seq) == len(qual)
    tseq = trim_adaptor(seq, adaptor, num_errors, right_side=right_side)
    if right_side:
        pos = seq.find(tseq)
    else:
        pos = seq.rfind(tseq)
    tqual = qual[pos:pos + len(tseq)]
    assert len(tseq) == len(tqual)
    return tseq, tqual


def create_unique_dir(path, limit=99):
    """Return a path to an empty directory. Either the dir at path, or a dir of the form 'path + _01'
    :param path: The initial path to use
    :param limit: The maximum number of directory variations this function will attempt to create.
    :return: A path to an empty directory.
    """
    width = len(str(limit))
    original = path.rstrip(os.sep)
    if len(os.listdir(original)) == 0:
        log.info("Using output directory - {}".format(path))
        return original  # folder empty, let's use it
    count = 1
    while count < limit:
        try:
            os.mkdir(path)
            log.info("Creating output directory - {}".format(path))
            return path
        except OSError as path_error:
            if path_error.errno == 17:  # file exists
                path = "{0}_{1:0>{2}}".format(original, count, width)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))


def trim_both_ends(seq_rec, primer_a, primer_b, num_errors, reverse_complement=False):
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
                                                            num_errors=num_errors,
                                                            right_side=False)
    if len(half_trimed_seq) < len(un_trimed_seq):
        found_a = True

    full_trimed_seq, full_trimed_qual = trim_adaptor_w_qual(half_trimed_seq,
                                                            half_trimed_qual,
                                                            adaptor=primer_b,
                                                            num_errors=num_errors,
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


def mask_seq_record(seq_rec, min_score, inplace=False):
    if not inplace:
        masked_seq_req = copy.deepcopy(seq_rec)
    else:
        masked_seq_req = seq_rec

    base_list = list(masked_seq_req.seq)
    for loc in range(len(base_list)):
        if masked_seq_req.letter_annotations["phred_quality"][loc] < min_score:
            base_list[loc] = 'N'

    masked_seq_req.seq = Seq("".join(base_list), alphabet=IUPAC.ambiguous_dna)

    return masked_seq_req


def trim_and_mask_seq_records(records, primer_a, primer_b, num_errors, min_nuc_score, min_avg_score, basename):
    for seq_rec in records:
        trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, num_errors, reverse_complement=False)

        # primers found in forword direction
        if trimed_seq_rec is not None:
            avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
            if avg_score < min_avg_score:
                log.info("seq excluded - avg_score:{:4.2f} < min_avg_score:{} - {} {}".format(avg_score, min_avg_score,
                                                                                              basename, seq_rec.id))
                continue
            yield mask_seq_record(trimed_seq_rec, min_nuc_score, inplace=True)

        # primers not found in forword direction
        else:
            trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, num_errors, reverse_complement=True)
            if trimed_seq_rec is not None:  # primers found in reverse direction
                avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
                if avg_score < min_avg_score:
                    log.info(
                        "seq excluded - avg_score:{:4.2f} < min_avg_score:{} - {} {}".format(avg_score, min_avg_score,
                                                                                             basename, seq_rec.id))
                    continue
                yield mask_seq_record(trimed_seq_rec, min_nuc_score, inplace=True)

            # primers not found in either direction
            else:
                log.info("seq excluded - primers not found - {} {}".format(basename, seq_rec.id))
                continue


def process_fastq(input_fn, output_dir, num_errors, min_nuc_score, min_avg_score, min_seqs=1, max_len=None,
                  aligner="muscle", pre_align_max_n=None, post_align_max_n=None):
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

    clean_records = trim_and_mask_seq_records(records, primer_a, primer_b, num_errors,
                                              min_nuc_score, min_avg_score, basename)

    clean_records = list(clean_records)
    if len(clean_records) < min_seqs:
        log.info(
            "alignment excluded - seq_count:{} < min_seqs:{} - {}".format(len(clean_records), min_seqs, basename))
        return

    # write clean fasta files
    with open(os.path.join(output_dir, basename) + ".clean.fasta", "wt") as f:
        for r in clean_records:
            n_count = r.seq.upper().count('N')
            if (pre_align_max_n is not None) and (n_count > pre_align_max_n):
                log.info("seq excluded - n_count:{} > pre_align_max_n:{} - {} {}".format(n_count, pre_align_max_n,
                                                                                         basename, r.id))
                continue
            f.write(r.format("fasta"))

    # align clean fasta
    if aligner == "muscle":
        cline = MuscleCommandline(input=os.path.join(output_dir, basename) + ".clean.fasta",
                                  out=os.path.join(output_dir, basename) + ".clean.aln", clw=True)
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as e:
            log.info("alignment failed - {} - {}".format(e, basename))
            return

    elif aligner == "clustalw":
        cline = ClustalwCommandline("clustalw2", infile=os.path.join(output_dir, basename) + ".clean.fasta")
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as e:
            log.info("alignment failed - {} - {}".format(e, basename))
            return

    else:
        log.info("alignment failed - invalid aligner: {} - {}".format(aligner, basename))
        return

    # parse aligned fasta
    with open(os.path.join(output_dir, basename) + ".clean.aln", "r") as f:
        alignment = AlignIO.read(f, "clustal", alphabet=IUPAC.ambiguous_dna)

    # take consensus of aligned fasta
    summary_align = AlignInfo.SummaryInfo(alignment)

    # write consensus fasta
    consensus = summary_align.gap_consensus(ambiguous="N")
    n_count = consensus.upper().count('N')
    if (post_align_max_n is not None) and (n_count > post_align_max_n):
        log.info("alignment excluded - n_count:{} > post_align_max_n:{} - {}".format(n_count, pre_align_max_n, basename))
        return
    seq = Seq(data=str(consensus), alphabet=IUPAC.ambiguous_dna)
    seq_rec = SeqRecord(seq=seq, id=basename.split("prime_")[0] + "prime", description="{}".format(len(alignment)))
    with open(os.path.join(output_dir, basename) + ".{}.consensus.fasta".format(len(alignment)), "wt") as f:
        f.write(seq_rec.format("fasta"))


def process_file_list(in_file_list, output_dir, num_errors, min_nuc_score, min_avg_score, min_seqs=1, max_len=None,
                      aligner="muscle", pre_align_max_n=None, post_align_max_n=None):
    # create pool
    with mp.Pool(min(len(in_file_list), mp.cpu_count())) as p:
        p.starmap(process_fastq, zip(in_file_list,
                                     itertools.repeat(output_dir),
                                     itertools.repeat(num_errors),
                                     itertools.repeat(min_nuc_score),
                                     itertools.repeat(min_avg_score),
                                     itertools.repeat(min_seqs),
                                     itertools.repeat(max_len),
                                     itertools.repeat(aligner),
                                     itertools.repeat(pre_align_max_n),
                                     itertools.repeat(post_align_max_n)
                                     ))


def main(_args):
    # if output_path does not exist, create it
    output_dir = _args.out_dir
    if not os.path.isdir(output_dir):
        log.info("Creating output directory - {}".format(output_dir))
        os.makedirs(output_dir)
    else:  # If it exists and is empty, use it. If it exists and is not empty, create and use new unique directory
        output_dir = create_unique_dir(output_dir)

    # call function based on batch or single mode
    if _args.in_file_list is None:
        log.info("Single file mode - {}".format(_args.in_file))
        process_fastq(_args.in_file, output_dir, _args.num_errors, _args.min_nuc_score, _args.min_avg_score,
                      min_seqs=_args.min_seqs, max_len=_args.max_len, aligner=_args.aligner,
                      pre_align_max_n=_args.pre_align_max_n, post_align_max_n=_args.post_align_max_n)
    else:
        log.info("Batch mode - {}".format(_args.in_file_list))
        with open(_args.in_file_list, "r") as f:
            in_file_list = [os.path.join(_args.in_dir, l.strip()) for l in f.readlines()]

        process_file_list(in_file_list, output_dir, _args.num_errors, _args.min_nuc_score, _args.min_avg_score,
                          min_seqs=_args.min_seqs, max_len=_args.max_len, aligner=_args.aligner,
                          pre_align_max_n=_args.pre_align_max_n, post_align_max_n=_args.post_align_max_n)


if __name__ == "__main__":

    # argparse configuration
    # noinspection PyMissingOrEmptyDocstring
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)


    parser = MyParser()
    parser.add_argument('-i', '--in_file', help='path to input file',
                        default=('Final.HQpolish99nomaxmin600.ccs.BX140414_001.5prime_'
                                 'CACAGAGACACGCACA.'
                                 'CGTCTCTATCTCTCTA.'
                                 'GCAGTCGAACATGTAGCTGACTCAGGTCACTCGCCTAAACTTCAGCCATT.'
                                 'TGATTYTTTGGACACCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCA.'
                                 'fastq'))
    parser.add_argument('-f', '--in_file_list', help='path to text file w/ an in_file filename on each line',
                        default=None)
    parser.add_argument('-b', '--in_dir', help='input dir of files named in in_file_list', default=None)
    parser.add_argument('-o', '--out_dir', help='path to output dir', default="output")
    parser.add_argument('-e', '--num_errors', type=int, help='number of errors allowed in primer match', default="2")
    parser.add_argument('-n', '--min_nuc_score', type=int, help='score below which a nuc is masked', default="60")
    parser.add_argument('-a', '--min_avg_score', type=int, help='avg score below which a seq is excluded', default="80")
    parser.add_argument('-s', '--min_seqs', type=int, help='min count of sequences before alignment excluded',
                        default=1)
    parser.add_argument('-m', '--max_len', type=int, help='max length overwhich a seq is excluded', default=None)
    parser.add_argument('-l', '--aligner', help='the alignment software to use. {clustalw, muscle}', default="muscle")
    parser.add_argument('-p', '--pre_align_max_n', type=int, help='number of Ns allowed in sequences to be aligned',
                        default="5")
    parser.add_argument('-t', '--post_align_max_n', type=int, help='number of Ns allowed in final aligned sequences',
                        default="5")

    args = parser.parse_args()

    # run main
    try:
        exit(main(args))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
