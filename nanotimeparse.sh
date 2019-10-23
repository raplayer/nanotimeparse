#!/bin/bash
# author:	Robert Player


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF

Help message for \`nanotimeparse.sh\`:

	Get subsets of (-i) Oxford Nanopore Technologies (ONT) basecalled fastq reads in slices of (-s) minutes, over a period of (-p) minutes. Input fastq file is output as 2 sets of n fasta files (n = p/s). Set 1 is n fasta files, and each file contains reads generated from the start of the ONT run to each time slice. Set 2 is also n fasta files, but each file contains only newly generated reads between each time slice.

NOTES:
	- for best results, p should be evenly divisible by s (i.e. p/s = INT)
	- using INT for minutes is preferrable, however, FLOAT is fine if p/s = INT
	- example inputs/outputs may be found in \`sandbox/\` (check log for details)
	- be sure to concatenate all fastq files output from a single ONT flowcell
	- dependencies: GNU Parallel [1], GNU Core Utils; see README.md for full list

USAGE:
	bash nanotimeparse.sh -t <threads> -i </absolute/path/to/nanopore_basecalled.fastq> -s <minutes> -p <minutes>

OPTIONS:
	-h      help		show this message
	-t	INT		number of threads to GNU parallel over
	-i      FASTQ		input basecalled nanopore fastq
	-s	INT or float	time slice in (minutes)
	-p	INT or float	period of time to slice up since start of sequencing run (minutes)


________________________________________________________________________________
References:
	1. O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.

EOF
}

# Function description:	make_since1970
#	arguments:	arg1 - one of the split fasta files derived from input fastq
#	actions:	Pulls start datetime of read generation from header of each read in the fasta file.
#				Converts each to a format usable by the 'date' GNU CoreUtil, stores datetimes in 'datetime.list'.
#					example, 'start_time=2018-08-05T06:35:49Z' > '2018-08-05 06:35:49'
#				Call 'date' function on the list, converts datetimes to Unix Epoch time (seconds since 01/01/1970).
#	output:		A concatenated list of all 'since1970' times for all reads in input.
#				These times are used as headers for reads in final script output, and are utilized in the make_diff function.
make_since1970()
{
	mkdir "$1-tmp"
	sed -e 's/.*start_time=//' -e 's/Z.*\t/\t/' -e 's/T/ /' "$1" | sort -T "$outdir/tmp/sort" > "$1-tmp/datetime.seq"
	cut -f1 "$1-tmp/datetime.seq" > "$1-tmp/datetime.list"
	date --file="$1-tmp/datetime.list" +%s > "$1-tmp/since1970.list"
	cat "$1-tmp/since1970.list" >> "$outdir/tmp/since1970.list"
	paste <(sort -T "$outdir/tmp/sort" "$1-tmp/since1970.list") <(sort -T "$outdir/tmp/sort" "$1-tmp/datetime.seq" | cut -f2) >> "$outdir/tmp/since1970.array"

}
export bin outdir
export -f make_since1970

# Function description:	make_diff
#	arguments:	arg1 - fasta file containing reads from $basetime to since1970 time in file name
#				arg2 - same fasta type but containing reads from $basetime to the next time slice
#	actions:	Finds non-overlapping reads in arg1 and arg2 fastas via 'comm -13' function.
#					comm:	Compare sorted files FILE1 and FILE2 line by line.
#					  -1              suppress column 1 (lines unique to FILE1)
#					  -2              suppress column 2 (lines unique to FILE2)	<- omitted
#					  -3              suppress column 3 (lines that appear in both files)
#	output:		A fasta file containing only newly generated reads between times indicated in the file names.
make_diff()
{
	t1=$(basename "$1" | sed -e 's/.*-//' -e 's/\.fasta//')
	t2=$(basename "$2" | sed -e 's/.*-//' -e 's/\.fasta//')
	printf "processing: $t1\t$t2\n"
	comm -13 <(grep "^>" "$1" | sort -T "$outdir/tmp/sort") <(grep "^>" "$2" | sort -T "$outdir/tmp/sort") > "$outdir/tmp/diff.from$t1-$t2"
	grep -A1 -f "$outdir/tmp/diff.from$t1-$t2" "$2"  | grep -v "\-\-" > "$outdir/from$t1-$t2.fasta"
}
export t1 t2 outdir
export -f make_diff



#	DEFAULTS & INPUTS & CHECKS
#===============================================================================
#	notes:
#		echo $? (0 = successful execution)
# absolute path to script dir
absolute_path_x="$(readlink -fn -- "$0"; echo x)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
if [[ $? != 0 ]]; then
	echo "Please locate the function 'dirname' and symlink it to this script's bin."
	echo "example: ln -s /usr/bin/dirname /full/path/to/nanoparse/bin/dirname"
	exit
fi


# parse args
while getopts "ht:i:s:p:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		i) INPUT=$OPTARG ;;
		s) SLICE=$OPTARG ;;
		p) PERIOD=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# check args
if [[ -z "$THREADS" ]]; then printf "%s\n" "Please specify number of threads (-t)."; exit; fi
if [[ -z "$INPUT" ]]; then printf "%s\n" "Please specify input fastq (-i)."; exit; fi
if [[ ! -f "$INPUT" ]]; then printf "%s\n" "The input (-i) $INPUT file does not exist."; exit; fi
if [[ -z "$SLICE" ]]; then printf "%s\n" "Please specify slice time (-s)."; exit; fi
if [[ -z "$PERIOD" ]]; then printf "%s\n" "Please specify period time (-p)."; exit; fi

# setup other variables
runtime=$(date +"%Y%m%d%H%M%S%N")
window_sec=$(printf "%s" "$SLICE" | awk '{print($0*60)}')
indir=$(dirname "$INPUT")
inbase=$(basename "$INPUT")
outdir="$indir/nanotimeparse-i_$inbase.s_$SLICE.p_$PERIOD"
if [[ -d "$outdir" ]]; then
	printf "Warning: It appears this input has already been parsed with these parameters.\nProcess terminating.\n"
	exit
else
	mkdir -p "$outdir/tmp/sort"
fi


# quick and dirty fastq checks
#	check 1 ensures that every 4th line starting with line 1 starts with the '@' character
#	check 2 is a very rough check of ONT format by searching for 'start_time=' in headers
#			if total reads == number of times 'start_time=' is found, we should be good to parse
fqcheck=$(awk 'NR % 4 == 1' "$INPUT" | cut -c1 | sort -T "$outdir/tmp/sort" | uniq -c | sed -e 's/ \+//' -e 's/ /\t/')
char=$(printf "$fqcheck" | cut -f2)
if [[ "$char" != "@" ]]; then
	printf "%s\n%s\n%s\n%s\t%s\n%s\n%s\n" "Check 1 error, all headers do not start with '@'." "The input (-i) does not appear to be fastq formatted." "Here's what we found:" "COUNT" "FIRST CHAR" "$fqcheck" "Process terminating."
	exit
fi
totalreads=$(awk 'END{print(NR/4)}' "$INPUT")
starttimes=$(grep -o "start_time=" "$INPUT" | wc -l)
if [[ "$totalreads" != "$starttimes" ]]; then
	printf "%s\n%s\n%s\n" "Check 2 error, there are more start times ($starttimes) than number of reads ($totalreads)." "The input (-i) does not appear to be ONT formatted." "Process terminating."
	exit
fi
printf "%s\n" "Checks completed successfully."

# print commands etc to log
log="$outdir/nanotimeparse.log-runtime_$runtime"
printf "\n%s\n" "Runtime $runtime" > "$log"
printf "%s\n" "nanotimeparse.sh -t $THREADS -i $INPUT -s $SLICE -p $PERIOD" >> "$log"
printf "%s\t%s\n" "Input fastq:" "$INPUT" >> "$log"
printf "%s\t%s\n" "Window(mins):" "$SLICE" >> "$log"
printf "%s\t%s\n" "window(sec):" "$window_sec" >> "$log"
printf "%s\t%s\n" "Period(mins):" "$PERIOD" >> "$log"


#	MAIN
#===============================================================================
# convert fastq to fasta2col, then split into 4000 reads per file
sed $'$!N;s/\\\n/\t/' "$INPUT" | sed $'$!N;s/\\\n/\t/' | cut -f1,2 > "$outdir/tmp/input.fasta2col"
sed -i 's/^@/>/' "$outdir/tmp/input.fasta2col"
split -l 10000 "$outdir/tmp/input.fasta2col" "$outdir/tmp/split_fasta2col."

# get list of all start timepoints for all reads
if [[ ! -f "$outdir/tmp/since1970.list" ]]; then
	printf "Making read start time list...\n"
	find "$outdir/tmp/" -name "split_fasta2col.*" -print0 | parallel -0 -n 1 -P "$THREADS" -I '{}' make_since1970 '{}'
else
	printf "The file 'since1970.list' already exists @ '$outdir/'\n"
fi

# find earliest time point
basetime=$(sort -T "$outdir/tmp/sort" "$outdir/tmp/since1970.list" | head -1)
printf "\trun start:\t$basetime\n"
printf "%s\t%s\n" "basetime:" "$basetime" >> "$log"

#	setup cut times
cut_start=$(printf "%0.0f\n" $(bc -l <<< "$basetime+$window_sec"))
printf "%s\t%s\n" "cut_start:" "$cut_start" >> "$log"
printf "\tcut start:\t$cut_start\n"
cut_stop=$(printf "%0.0f\n" $(bc -l <<< "$basetime+($PERIOD*60)"))
printf "%s\t%s\n" "cut_stop:" "$cut_stop" >> "$log"
printf "\tcut stop:\t$cut_stop\n"

# make file of each time slice
if [[ ! -f "$outdir/set1.complete" ]]; then
	printf "Generating set 1 files...\n"
	seq "$cut_start" "$window_sec" "$cut_stop" > "$outdir/tmp/since1970.slicetimes"
	# make separate fasta for each 'time slice'
	while read time; do
		printf "processing: $basetime\t$time\n"
		awk -F'\t' -v time="$time" '{
			fnrheader[NR]=$1
			fnrseq[NR]=$2
		}END{
			for(i=1;i<length(fnrheader);i++){
				if(fnrheader[i]<time){
					printf(">%s\n%s\n",fnrheader[i],fnrseq[i])
				}
			}
		}' "$outdir/tmp/since1970.array" > "$outdir/fromStart-$time.fasta"
	done < "$outdir/tmp/since1970.slicetimes"
	touch "$outdir/set1.complete"
else
	printf "Set 1 already generated, moving on.\n"
fi


# only output reads produced during the time slice
# from 'basetime' to 'cut_start' (the first fasta)
if [[ ! -f "$outdir/set2.complete" ]]; then
	printf "Generating set 2 files...\n"
	printf "processing: $basetime\t$cut_start\n"
	cat "$outdir/fromStart-$cut_start.fasta" > "$outdir/from$basetime-$cut_start.fasta"
	# and the rest
	parallel --xapply --jobs="$THREADS" make_diff \
		::: $(find "$outdir" -maxdepth 1 -name "fromStart*fasta" | sort -T "$outdir/tmp/sort" | head -n-1) \
		::: $(find "$outdir" -maxdepth 1 -name "fromStart*fasta" | sort -T "$outdir/tmp/sort" | tail -n+2)
	touch "$outdir/set2.complete"
else
	printf "Set 2 already generated, moving on.\n"
fi


printf "Checking total output files...\n"
expected=$(printf "%0.0f\n" $(bc -l <<< "2*($PERIOD/$SLICE)"))
actual=$(find $outdir -mindepth 1 -maxdepth 1 -type f -name "from*.fasta" | wc -l)
if [[ "$expected" == "$actual" ]]; then
	printf "\t2*(p/s) == expected == actual :: 2*($PERIOD/$SLICE) == $expected == $actual\n"
	printf "\tnanotimeparse.sh completed successfully, yay!\n"
else
	printf "\t2*(p/s) != expected != actual :: 2*($PERIOD/$SLICE) != $expected != $actual\n"
	printf "\tsomething terrible happened, or p/s != INT and this output is what you expected :)\n"
fi
