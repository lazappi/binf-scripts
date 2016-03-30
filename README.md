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

**Language**: Python

**Details:** Adapted from the `Celloline` [stats script](https://github.com/Teichlab/celloline/blob/master/lib/stats.py).

**Usage**: `alignStats.py [-h] -o OUT [-g GTF] [-i INDEX] -t {sam,bam}
           [-p PARALLEL] SAM/BAM [SAM/BAM ...]`

sraDownload.R
-------------

**Description:** Download files from SRA using ASCP. 

**Language**: R

**Usage**: `sraDownload.R [-h] [-d DATABASE] [-o OUT] [-a ASCPCMD] [-t TYPE]
           SRA [SRA ...]`
