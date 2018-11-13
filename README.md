# ccs-consensuser

ccs-consensuser.py provides a number of options alowing a user combine multiple ccs fastq results.

features include:

 - uses primers from filename to locate and orient amplicon
 - masks low quality bases to exclude from consensus calculations
 - filters outlying reads based on several paramaters
 - emits log entries for every base or read filtered
 - aligns reads
 - calls consensus
 - designed for single or batch processing


 ----------------

 ## example

```
ccs-consensuser.py \
--in_dir ../fastq \
--in_file_list ../3prime_filelists/${FILELIST} \
--out_dir output \
--overwrite \
--primer_mismatch 2 \
--min_base_score 60 \
--mask_char N \
--sequence_max_mask 5 \
--min_seq_count 1 \
--aligner muscle \
--expected_length 812 \
--max_len_delta 9 \
--consensus_ignore_mask_char \
--consensus_threshold .51 \
--consensus_ambiguous_char n \
--alignment_max_amb 5
```


 ----------------

 ## help

 usage: ccs-consensuser.py [-h] [-l ALIGNER] [--mask_char MASK_CHAR]
                          [-i IN_FILE] [-o OUT_DIR] [--overwrite]
                          [--in_file_list IN_FILE_LIST] [--in_dir IN_DIR]
                          [-n MIN_BASE_SCORE] [-e PRIMER_MISMATCH]
                          [-p SEQUENCE_MAX_MASK] [-d MAX_LEN_DELTA]
                          [--expected_length EXPECTED_LENGTH]
                          [--max_len MAX_LEN] [--min_seq_score MIN_SEQ_SCORE]
                          [-s MIN_SEQ_COUNT] [-t ALIGNMENT_MAX_AMB]
                          [--consensus_threshold CONSENSUS_THRESHOLD]
                          [--consensus_require_multiple]
                          [--consensus_ambiguous_char CONSENSUS_AMBIGUOUS_CHAR]
                          [--consensus_ignore_mask_char]

optional arguments:
  -h, --help            show this help message and exit
  -l ALIGNER, --aligner ALIGNER
                        the alignment software to use. {clustalw, muscle}
  --mask_char MASK_CHAR
                        char used to mask bases with quality below threshold
  -i IN_FILE, --in_file IN_FILE
                        path to input file
  -o OUT_DIR, --out_dir OUT_DIR
                        path to output dir
  --overwrite           set this flag to disable creating new out dir
  --in_file_list IN_FILE_LIST
                        path to text file w/ an in_file filename on each line
  --in_dir IN_DIR       input dir of files named in in_file_list
  -n MIN_BASE_SCORE, --min_base_score MIN_BASE_SCORE
                        score below which a nuc is masked
  -e PRIMER_MISMATCH, --primer_mismatch PRIMER_MISMATCH
                        number of errors allowed in primer match
  -p SEQUENCE_MAX_MASK, --sequence_max_mask SEQUENCE_MAX_MASK
                        number of mask_char allowed in sequences to be aligned
  -d MAX_LEN_DELTA, --max_len_delta MAX_LEN_DELTA
                        allowed variation from mode of sequence length
  --expected_length EXPECTED_LENGTH
                        optional, replaces average in max_len_delta filter
  --max_len MAX_LEN     max length overwhich a seq is excluded
  --min_seq_score MIN_SEQ_SCORE
                        avg score below which a seq is excluded
  -s MIN_SEQ_COUNT, --min_seq_count MIN_SEQ_COUNT
                        min count of sequences before alignment excluded
  -t ALIGNMENT_MAX_AMB, --alignment_max_amb ALIGNMENT_MAX_AMB
                        number of consensus_ambiguous_char allowed in final
                        consensus sequences
  --consensus_threshold CONSENSUS_THRESHOLD
                        proportion threshold for consensus to call a base.
  --consensus_require_multiple
                        require multiple instanses for consensus to call a
                        base
  --consensus_ambiguous_char CONSENSUS_AMBIGUOUS_CHAR
                        char representing bases under consensus threshold
  --consensus_ignore_mask_char
                        discount mask char when calculating consensus
