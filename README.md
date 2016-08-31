My bioinformatics scripts
=========================

A collection of useful bioinformatics scripts

starStats.py
------------

**Description:** Extract alignment statistics from STAR log files.

**Language:** Python

**Usage:** `starStats.py [-h] -o OUT LOGFILE [LOGFILE ...]`

alignStats.py
-------------

**Description:** Extract alignment statistics from SAM/BAM files.

**Language:** Python

**Details:** Adapted from the `Celloline` [stats script](https://github.com/Teichlab/celloline/blob/master/lib/stats.py).
More details [here](http://lazappi.id.au/extracting-alignment-statistics-using-python/).

**Usage**: `alignStats.py [-h] -o OUT [-g GTF] [-i INDEX] -t {sam,bam}
           [-p PARALLEL] SAM/BAM [SAM/BAM ...]`

sraDownload.R
-------------

**Description:** Download files from SRA using ASCP.

**Language:** R

**Usage:** `sraDownload.R [-h] [-d DATABASE] [-o OUT] [-a ASCPCMD] [-t TYPE]
           SRA [SRA ...]`

sampleFastq.R
-------------

**Description:** Sample reads from FASTQ files.

**Language:** R

**Usage:** `sampleFastq [-h] -1 PAIR1 [PAIR1 ...] [-2 [PAIR2 [PAIR2 ...]]]
           [-n NREADS] [-s SEED] [-o OUTPATH]`

mergeFastqs.py
-------------

**Description:** Merge FASTQ files based on a filename pattern.

**Language:** Python

**Usage:** `mergeFastqs [-h] [-s SEPARATOR] [-o OUTDIR] -p PATTERN FASTQ
           [FASTQ ...]`

kallistoMulti.py
----------------

**Description:** Run kallisto over multiple FASTQ files.

**Language:** Python

**Usage:** `kallistoMulti [-h] -i INDEX -o OUTPUT_DIR [-1 READ1 [READ1 ...]]
            [-2 READ2 [READ2 ...]]`

kallistoMerge.R
---------------

**Description:** Merge output from multiple kallisto runs.

**Language:** R

**Usage:** `kallistoMerge [-h] -d DIRECTORY -t TX2GENE [-o OUTPATH]`
