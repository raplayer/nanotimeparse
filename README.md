# nanotimeparse

This is the repository for nanotimeparse, which can parse an Oxford Nanopore fastq file by read generation times.

Get subsets of (-i) Oxford Nanopore Technologies (ONT) basecalled fastq reads in slices of (-s) minutes, over a period of (-p) minutes. Input fastq file is output in n fasta files (n = p/s). Using INT for minutes is preferrable, although you may get away with using FLOAT if you're lucky (particularly if p/s = INT, e.g. n = 4/0.25 = 16). For best results, concatenate all fastq files from a single ONT flowcell

# Installation

Clone this repo with:

`git clone https://github.com/player5858/

# Dependencies

*GNU CoreUtils

# Usage

bash nanotimeparse.sh -t 5 -i </absolute/path/to/nanopore.fastq> -s <minutes> -p <minutes>

Output root:

/<input directory path>/nanotimeparse-$(basename $fastq)/

OPTIONS:

	-h      help		show this message
	-t	INT		number of threads to GNU parallel over
	-i      FASTQ		input basecalled nanopore fastq
	-s	INT or float	time slice in (minutes)
	-p	INT or float	period of time to slice up since start of sequencing run (minutes)



________________________________________________________________________________

References:

1. O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.

