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
	- example output, see \`sandbox/nanotimeparse-test.fq/\`, and check log for details
	- dependencies: GNU Parallel [1], GNU Core Utils; see README.md for full list
	- using INT for minutes is preferrable, FLOAT may work if p/s = INT, e.g. n = 4/0.25 = 16
	- for best results, concatenate all fastq files from a single ONT flowcell

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

make_since1970()
{
	$bin/mkdir "$1-tmp"
	$bin/sed -e 's/.*start_time=//' -e 's/Z flow_cell.*\t/\t/' -e 's/T/ /' "$1" | $bin/sort -T "$outdir/tmp/sort" > "$1-tmp/datetime.seq"
	$bin/cut -f1 "$1-tmp/datetime.seq" > "$1-tmp/datetime.list"
	$bin/date --file="$1-tmp/datetime.list" +%s > "$1-tmp/since1970.list"
	$bin/cat "$1-tmp/since1970.list" >> "$outdir/tmp/since1970.list"
	paste <($bin/sort -T "$outdir/tmp/sort" "$1-tmp/since1970.list") <($bin/sort -T "$outdir/tmp/sort" "$1-tmp/datetime.seq" | $bin/cut -f2) >> "$outdir/tmp/since1970.array"

}
export bin outdir
export -f make_since1970

make_diff()
{
	t1=$($bin/basename "$1" | $bin/sed -e 's/.*-//' -e 's/\.fasta//')
	t2=$($bin/basename "$2" | $bin/sed -e 's/.*-//' -e 's/\.fasta//')
	printf "processing: $t1\t$t2\n"
	$bin/comm -13 <($bin/grep "^>" "$1" | $bin/sort -T "$outdir/tmp/sort") <($bin/grep "^>" "$2" | $bin/sort -T "$outdir/tmp/sort") > "$outdir/tmp/diff.from$t1-$t2"
	$bin/grep -A1 -f "$outdir/tmp/diff.from$t1-$t2" "$2"  | $bin/grep -v "\-\-" > "$outdir/from$t1-$t2.fasta"
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
bin="$scriptdir/bin"

# create bin of paths to dependencies
if [[ ! -f "$scriptdir/install.complete" ]]; then
	echo "installing..."
	if [[ -d "$bin" ]]; then rm -r "$bin"; fi
	mkdir "$bin"
	whereis "cat" "mkdir" "sed" "sort" "cut" "date" "paste" "basename" "printf" "comm" "grep" "awk" "uniq" "dirname" "split" "find" "parallel" "head" "bc" "seq" "tail" "wc" | cut -f2 -d' ' | while read path; do
		util=$(basename "$path")
		ln -s "$path" "$bin/$util"
	done
	touch "$scriptdir/install.complete"
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
if [[ -z "$THREADS" ]]; then $bin/printf "%s\n" "Please specify number of threads (-t)."; exit; fi
if [[ -z "$INPUT" ]]; then $bin/printf "%s\n" "Please specify input fastq (-i)."; exit; fi
if [[ ! -f "$INPUT" ]]; then $bin/printf "%s\n" "The input (-i) $INPUT file does not exist."; exit; fi
if [[ -z "$SLICE" ]]; then $bin/printf "%s\n" "Please specify slice time (-s)."; exit; fi
if [[ -z "$PERIOD" ]]; then $bin/printf "%s\n" "Please specify period time (-p)."; exit; fi

# quick and dirty fastq checks
#	check 1 ensures that every 4th line starting with line 1 starts with the '@' character
#	check 2 is a very rough check of ONT format by searching for 'start_time=' in headers
#			if total reads == number of times 'start_time=' is found, we should be good to parse
fqcheck=$($bin/awk 'NR % 4 == 1' "$INPUT" | $bin/cut -c1 | $bin/sort -T "$outdir/tmp/sort" | $bin/uniq -c | $bin/sed -e 's/ \+//' -e 's/ /\t/')
char=$($bin/printf "$fqcheck" | $bin/cut -f2)
if [[ "$char" != "@" ]]; then
	$bin/printf "%s\n%s\n%s\n%s\n%s\n" "Check 1 error, all headers do not start with '@'." "The input (-i) does not appear to be fastq formatted." "Here's what we found (<count>\t<first character>):" "$fqcheck" "Process terminated."
	exit
fi
totalreads=$($bin/awk 'END{print(NR/4)}' "$INPUT")
starttimes=$($bin/grep -o "start_time=" "$INPUT" | $bin/wc -l)
if [[ "$totalreads" != "$starttimes" ]]; then
	$bin/printf "%s\n%s\n%s\n" "Check 2 error, there are more start times ($starttimes) than number of reads ($totalreads)." "The input (-i) does not appear to be ONT formatted." "Process terminated."
	exit
fi
$bin/echo "Checks completed successfully."
# if checks ok, rename some inputs
fastq="$INPUT"
window_min="$SLICE"
period="$PERIOD"

# setup other variables
runtime=$($bin/date +"%Y%m%d%H%M%S%N")
window_sec=$($bin/printf "%s" "$window_min" | $bin/awk '{print($0*60)}')
indir=$($bin/dirname "$fastq")
inbase=$($bin/basename "$fastq")
outdir="$indir/nanotimeparse-$inbase"
$bin/mkdir -p "$outdir/tmp/sort"

# print commands etc to log
$bin/printf "\n%s\n" "Runtime $runtime" >> "$outdir/nanotimeparse.log"
$bin/printf "%s\n" "nanotimeparse.sh -t $THREADS -i $INPUT -s $SLICE -p $PERIOD" >> "$outdir/nanotimeparse.log"
$bin/printf "%s\t%s\n" "Input fastq:" "$fastq" >> "$outdir/nanotimeparse.log"
$bin/printf "%s\t%s\n" "Window(mins):" "$window_min" >> "$outdir/nanotimeparse.log"
$bin/printf "%s\t%s\n" "window(sec):" "$window_sec" >> "$outdir/nanotimeparse.log"
$bin/printf "%s\t%s\n" "Period(mins):" "$period" >> "$outdir/nanotimeparse.log"


#	MAIN
#===============================================================================
# convert fastq to fasta2col, then split into 4000 reads per file
$bin/sed $'$!N;s/\\\n/\t/' "$fastq" | $bin/sed $'$!N;s/\\\n/\t/' | $bin/cut -f1,2 > "$outdir/tmp/input.fasta2col"
$bin/split -l 10000 "$outdir/tmp/input.fasta2col" "$outdir/tmp/split_fasta2col."

# get list of all start timepoints for all reads
if [[ ! -f "$outdir/tmp/since1970.list" ]]; then
	$bin/printf "Making read start time list...\n"
	$bin/find "$outdir/tmp/" -name "split_fasta2col.*" -print0 | parallel -0 -n 1 -P "$THREADS" -I '{}' make_since1970 '{}'
else
	$bin/printf "The file 'since1970.list' already exists @ '$outdir/'\n"
fi

# find earliest time point
basetime=$($bin/sort -T "$outdir/tmp/sort" "$outdir/tmp/since1970.list" | $bin/head -1)
$bin/printf "\trun start:\t$basetime\n"
$bin/printf "%s\t%s\n" "basetime:" "$basetime" >> "$outdir/nanotimeparse.log"
#printf "%s\t%s\n" "lasttime:" "$lasttime" >> "$outdir/nanotimeparse.log"

#	setup cut times
# iterate starting at $basetime+window_sec ($2*60), out to $basetime+period_sec ($3*60)
cut_start=$($bin/printf "%0.0f\n" $($bin/bc -l <<< "$basetime+$window_sec"))
$bin/printf "%s\t%s\n" "cut_start:" "$cut_start" >> "$outdir/nanotimeparse.log"
$bin/printf "\tcut start:\t$cut_start\n"
cut_stop=$($bin/printf "%0.0f\n" $($bin/bc -l <<< "$basetime+($period*60)"))
$bin/printf "%s\t%s\n" "cut_stop:" "$cut_stop" >> "$outdir/nanotimeparse.log"
$bin/printf "\tcut stop:\t$cut_stop\n"

# make file of each time slice
if [[ ! -f "$outdir/set1.complete" ]]; then
	$bin/printf "Generating set 1 files...\n"
	$bin/seq "$cut_start" "$window_sec" "$cut_stop" > "$outdir/tmp/since1970.slicetimes"
	# make separate fasta for each 'time slice'
	while read time; do $bin/awk -F'\t' -v time="$time" '{fnrheader[NR]=$1; fnrseq[NR]=$2}END{for(i=1;i<length(fnrheader);i++){if(fnrheader[i]<time){printf(">%s\n%s\n",fnrheader[i],fnrseq[i])}}}' "$outdir/tmp/since1970.array" > "$outdir/fromStart-$time.fasta"; done < "$outdir/tmp/since1970.slicetimes"
	touch "$outdir/set1.complete"
else
	$bin/printf "Set 1 already generated, moving on.\n"
fi


# only output reads produced during the time slice
# from 'basetime' to 'cut_start' (the first fasta)
if [[ ! -f "$outdir/set2.complete" ]]; then
	$bin/printf "Generating set 2 files...\n"
	$bin/cat "$outdir/fromStart-$cut_start.fasta" > "$outdir/from$basetime-$cut_start.fasta"
	# and the rest
	$bin/parallel --xapply --jobs="$THREADS" make_diff \
		::: $($bin/find "$outdir" -maxdepth 1 -name "fromStart*fasta" | $bin/sort -T "$outdir/tmp/sort" | $bin/head -n-1) \
		::: $($bin/find "$outdir" -maxdepth 1 -name "fromStart*fasta" | $bin/sort -T "$outdir/tmp/sort" | $bin/tail -n+2)
	touch "$outdir/set2.complete"
else
	$bin/printf "Set 2 already generated, moving on.\n"
fi


$bin/printf "Checking total output files...\n"
expected=$($bin/printf "%0.0f\n" $($bin/bc -l <<< "$PERIOD/$SLICE"))
actual=$($bin/find $outdir -maxdepth 1 -name "fromStart*fasta" | $bin/wc -l)
if [[ "$expected" == "$actual" ]]; then
	$bin/printf "\t(p/s) == expected == actual :: ($PERIOD/$SLICE) == $expected == $actual\n"
	$bin/printf "\tnanotimeparse.sh completed successfully, yay!\n"
else
	$bin/printf "\t(p/s) != expected != actual :: ($PERIOD/$SLICE) != $expected != $actual\n"
	$bin/printf "\tsomething's not right, booo...\n"
fi
