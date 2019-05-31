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
consensuser.py \
--in_file_list filelists/${FILELIST} \
--out_dir output \
--sleep_time 15 \
--aligner muscle \
--mask_char N \
--min_base_score 60 \
--sequence_max_mask 5 \
--min_seq_count 3 \
--consensus_diploid_mode \
--consensus_threshold .7 \
--consensus_diploid_threshold .3 \
--consensus_ignore_mask_char \
--consensus_ambiguous_char n \
--alignment_max_amb 5
```


 ----------------

 ## help
 
```
usage: consensuser.py [-h] [-i IN_FILE] [-o OUT_DIR]
                      [--in_file_list IN_FILE_LIST] [--sleep_time SLEEP_TIME]
                      [-l ALIGNER] [--mask_char MASK_CHAR] [-n MIN_BASE_SCORE]
                      [-p SEQUENCE_MAX_MASK] [-d MAX_LEN_DELTA]
                      [--expected_length EXPECTED_LENGTH] [--max_len MAX_LEN]
                      [--min_seq_score MIN_SEQ_SCORE] [-s MIN_SEQ_COUNT]
                      [-t ALIGNMENT_MAX_AMB] [--consensus_diploid_mode]
                      [--consensus_threshold CONSENSUS_THRESHOLD]
                      [--consensus_diploid_threshold CONSENSUS_DIPLOID_THRESHOLD]
                      [--consensus_require_multiple]
                      [--consensus_ambiguous_char CONSENSUS_AMBIGUOUS_CHAR]
                      [--consensus_ignore_mask_char]
                      [--primer_regex PRIMER_REGEX]

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --in_file IN_FILE
                        path to input file (default: diploid.test_ATAGCGACGCGA
                        TATA.AGCGTCTCGCATCATG.TYTCAACDAAYCAYAAAGATATTGA.TAATAT
                        GGCAGATTAGTGCAATGGA.fastq)
  -o OUT_DIR, --out_dir OUT_DIR
                        path to output dir (default: output)
  --in_file_list IN_FILE_LIST
                        path to text file w/ an in_file filename on each line
                        (default: None)
  --sleep_time SLEEP_TIME
                        max length overwhich a seq is excluded (default: None)
  -l ALIGNER, --aligner ALIGNER
                        the alignment software to use. {clustalw, muscle}
                        (default: muscle)
  --mask_char MASK_CHAR
                        char used to mask bases with quality below threshold
                        (default: N)
  -n MIN_BASE_SCORE, --min_base_score MIN_BASE_SCORE
                        score below which a nuc is masked (default: 60)
  -p SEQUENCE_MAX_MASK, --sequence_max_mask SEQUENCE_MAX_MASK
                        number of mask_char allowed in sequences to be aligned
                        (default: 5)
  -d MAX_LEN_DELTA, --max_len_delta MAX_LEN_DELTA
                        allowed variation from mode of sequence length
                        (default: 5)
  --expected_length EXPECTED_LENGTH
                        optional, replaces average in max_len_delta filter
                        (default: None)
  --max_len MAX_LEN     max length overwhich a seq is excluded (default: None)
  --min_seq_score MIN_SEQ_SCORE
                        avg score below which a seq is excluded (default:
                        None)
  -s MIN_SEQ_COUNT, --min_seq_count MIN_SEQ_COUNT
                        min count of sequences before alignment excluded
                        (default: 5)
  -t ALIGNMENT_MAX_AMB, --alignment_max_amb ALIGNMENT_MAX_AMB
                        number of consensus_ambiguous_char allowed in final
                        consensus sequences (default: 5)
  --consensus_diploid_mode
                        use diploid mode for consensus calling (default:
                        False)
  --consensus_threshold CONSENSUS_THRESHOLD
                        proportion threshold for consensus to call a base.
                        (default: .7)
  --consensus_diploid_threshold CONSENSUS_DIPLOID_THRESHOLD
                        proportion threshold for consensus to call a base.
                        (default: .3)
  --consensus_require_multiple
                        require multiple instanses for consensus to call a
                        base (default: False)
  --consensus_ambiguous_char CONSENSUS_AMBIGUOUS_CHAR
                        char representing bases under consensus threshold
                        (default: n)
  --consensus_ignore_mask_char
                        discount mask char when calculating consensus
                        (default: False)
  --primer_regex PRIMER_REGEX
                        regex to parse primers from filename. default assumes
                        sampleName_barcodeA.bardoceB.primerA.primerB.extension
                        (default: _[^._]+?\.[^.]+?\.([^.]+?)\.([^.]+?)$)

```
