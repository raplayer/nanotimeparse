# nanotimeparse

This is the repository for nanotimeparse. It parses an Oxford Nanopore fastq file on read sequencing start times.

## Description

	Get subsets of (-i) Oxford Nanopore Technologies (ONT) basecalled fastq reads in slices of (-s) minutes, over a period of (-p) minutes. Input fastq file is output as 2 sets of n fasta files (n = p/s). Set 1 is n fasta files, and each file contains reads generated from the start of the ONT run to each time slice. Set 2 is also n fasta files, but each file contains only newly generated reads between each time slice.

## Notes

	- example output, see \`sandbox/nanotimeparse-test.fq/\`, and check log for details
	- dependencies: GNU Parallel [1], GNU Core Utils; see README.md for full list
	- using INT for minutes is preferrable, FLOAT may work if p/s = INT, e.g. n = 4/0.25 = 16
	- for best results, concatenate all fastq files from a single ONT flowcell

## Usage

'nanotimeparse.sh -t <threads> -i </absolute/path/to/nanopore_basecalled.fastq> -s <minutes> -p <minutes>'

	-h	help			help message
	-t	INT				number of threads to GNU parallel over
	-i	FASTQ			input basecalled nanopore fastq
	-s	INT or float	time slice in (minutes)
	-p	INT or float	period of time to slice up since start of sequencing run (minutes)


## Dependencies

GNU CoreUtils: cat, mkdir, sed, sort, cut, date, paste, basename, printf, comm, grep, awk, uniq, dirname, split, find, parallel [1], head, bc, seq, tail, wc


## Installation

Clone this repo with:

`git clone https://github.com/raplayer/nanotimeparse.git`

If you'd like to call the tool globally, symbolically link the shell script into a $PATH path.
Example:

`sudo ln -s $PWD/nanotimepare/nanotimeparse.sh /usr/local/bin`


## References
	1. O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.

