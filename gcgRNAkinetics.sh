#!/usr/bin/env bash

###############################################################################
# Runs a complete RNA kinetics study for a given RNA sequence, i.e.
# - enumerating the (partial) micro-state energy landscape using indel-moves
# - computing the gradient-basing partitioning (level 1) using barriers
# - computing further coarse grainings (level >1) via gcgBarriers.pl
# - computing probability trajectories using treekin for each coarse graining
# - generating pdf plots of the trajectories using R
#
# arguments:
# (1) RNA sequence of file that contains a SINGLE RNA sequence
# (2) (optional) maximal energy of RNA structures to be considered
# (3) (optional) maximal number of time steps for trajectory computations
#
# author Martin Raden - 2018
###############################################################################

# maximal memory available for 'sort' command. 
SORTMAXMEM=20G # if exceeded temporary files in current directory are used.

# maximal absolute energy (in kcal/mol) of micro-states to be considered 
MAXE=9999; 

# maximal time for kinetics computation
MAXTIME=100000000000; # maximal time (step) to compute

# ID of the local minimum that is populated 100% at trajectory start (p0)
# a value of 0 indicates the open chain 
STARTSTATE=0;

# exit function
die() { echo -e "$*" 1>&2 ; exit 1; }

# get RNA sequence
RNA=""
PREFIX=$RNA
[[ ($# -eq 0 || $# -gt 3) ]] && die "\n usage : gcgRNAkinetics.sh <RNA sequence (file)> <MAXE=$MAXE> <MAXTIME=$MAXTIME>\n";
if [ -f $1 ]; then
	# get file content: should contain only one RNA (consecutive) sequence
	# trim leading/trailing whitespaces
	PREFIX=$1
	RNA=$(cat $PREFIX | awk '{if(NF>0){printf $1;for(i=2;i<=NF;i++){printf " "$i}}}');
else
	# get RNA sequence from argument
	RNA=$1;
	PREFIX=$RNA
fi
# check sequence input
[[ ( $RNA =~ ^[ACGU]+$ ) ]] || die "\nERROR : RNA sequence not valid (single sequence of [ACGU]* expected)\n";
# check command line arguments
# (optional) check and set maximal energy of micro-states for kinetics computation
if [ $# -gt 1 ]; then
	# check
	[[ "x$2" =~ ^x[-+]?[0-9]+\.?[0-9]*$ ]] || die "\nERROR : MAXE (2nd argument) is not a floating point number like '5.4'\n"
	# store
	MAXE=$2
	# upper bound check (to avoid enumeration problems with RNAsubopt)
	MAXMAXE=10000
	[[ ${MAXE%%.*} -lt $MAXMAXE ]] || die "\nERROR : MAXE (2nd argument) has to be < $MAXMAXE\n"
fi
# (optional) check and set maximal time for kinetics computation 
if [ $# -gt 2 ]; then
	# check
	[[ "x$3" =~ ^x[0-9]+$ ]] || die "\nERROR : maximum number of trajectory time steps (3rd argument) has to be an integer\n"
	# store
	MAXTIME=$3
fi

# whether or not R is available for plotting
Ravailable=1;
# ensure required binaries are in PATH
[[ $(type -P "RNAfold") ]] || die "\nERROR : RNAfold not in PATH\n";
[[ $(type -P "RNAsubopt") ]] || die "\nERROR : RNAsubopt not in PATH\n";
[[ $(type -P "barriers") ]] || die "\nERROR : barriers not in PATH\n";
[[ $(type -P "treekin") ]] || die "\nERROR : treekin not in PATH\n";
[[ $(type -P "gcgBarriers.pl") ]] || die "\nERROR : gcgBarriers.pl not in PATH or not executable\n";
[[ $(type -P "gcgBarriersMfeComponent.pl") ]] || die "\nERROR : gcgBarriersMfeComponent.pl not in PATH or not executable\n";
if [ ! $(type -P "R") ]; then
	echo -e "\nWARNING : R not in PATH => skipping figure generation\n"; 
	Ravailable=0;
fi

[[ $STARTSTATE =~ ^[0-9]+$ ]] || die "\nERROR : index of start state has to be >= 0\n";

echo "##############  MFE COMPUTATION  #################"

# get mfe to compute deltaE for RNAsubopt call
MFE=$(echo $RNA | RNAfold --noPS | awk 'NR==2{print $NF}' | tr -d "()" | tr -d " ");
# compute and check (floating point) deltaE for RNAsubopt call (has to be positive)
[[ $(echo "if( $MAXE < $MFE ) 1 else 0" | bc) == 1 ]] && die "\nERROR : MAXE < MFE\n";
DELTAE=$(echo "$MAXE - $MFE" | bc)

echo "##############  LEVEL 0 ENUMERATION  #################"

if [ ! -f $PREFIX.RNAsubopt.bz2 ]; then

# enumerate (and sort) all secondary structures (first by energy using structure string for tie breaking)
echo $RNA | RNAsubopt --deltaEnergy=$DELTAE | sort -k 2,2n -k 1,1dr -S $SORTMAXMEM -T $PWD | bzip2 > $PREFIX.RNAsubopt.bz2;

fi

# check if state enumeration was successful
[[ $(bzcat $PREFIX.RNAsubopt.bz2 | grep -m 1 "A" -c) == 1 ]] || die "\nERROR : RNAsupopt output $PREFIX.RNAsubopt.zip seems to be empty\n";


echo "##############  LEVEL 1 COARSE GRAINING  #################"

if [ ! -f $PREFIX.barriers.out ]; then

# create temporary subdir to avoid parallel file generation
CURPWD=$PWD
TMPDIR=$(mktemp -d -p $CURPWD)
echo -e "\nINFO : creating barriers temporary directory $TMPDIR for $PREFIX";
cd $TMPDIR  
# run barriers to compute coarse graining level 1
bzcat $CURPWD/$PREFIX.RNAsubopt.bz2 | barriers --rates -G RNA -M noShift --bsize --max=999999 --minh=0 > $CURPWD/$PREFIX.barriers.out;
# store rate matrix generated by barriers
mv rates.out $CURPWD/$PREFIX.barriers.rates;
# go back to previous working directoy
cd $CURPWD
# cleanup obsolete barriers files and temporary directory
rm -rf $TMPDIR

fi

# get open chain ID
OCID=$(grep -P "^\\s*\\d+\\s+[\\.]+\\s+0" $PREFIX.barriers.out | awk 'NR==1 {print $1}');

echo "##############  MFE COMPONENT REDUCTION  #################"

if [ ! -f $PREFIX.barriers.out.all ]; then

# get connected component of mfe structure
gcgBarriersMfeComponent.pl $PREFIX.barriers.out $PREFIX.barriers.rates;
# rename files for further processing
mv -f $PREFIX.barriers.out $PREFIX.barriers.out.all
ln -s $PREFIX.barriers.out.mfeComp.barriers $PREFIX.barriers.out
mv -f $PREFIX.barriers.rates $PREFIX.barriers.rates.all
ln -s $PREFIX.barriers.out.mfeComp.rates $PREFIX.barriers.rates

fi

echo "##############  LEVEL >1 COARSE GRAINING  #################"

if [ ! -f $PREFIX.barriers.out.gcgBarriers.out ]; then

# compute generalized coarse grainings
gcgBarriers.pl $PREFIX.barriers.out $PREFIX.barriers.rates > $PREFIX.barriers.out.gcgBarriers.out

fi

echo "################  TRAJECTORY COMPUTATION AND PLOTTING  ######################"


# define variables for the following shell function
# shell function (tested in bash) to compute one plot via treekin 
function runTreekin {
	CURFILE=$1 # base file name without extension (.barriers .rates)
	CURLVL=$2  # the coarse graining level
	CURSTARTID=$3 # id of start state to be taken from according barriers output file
	CURSTARTROW=$4 # row number of the start state

echo "################  TRAJECTORY LEVEL $CURLVL : $CURFILE $CURSTARTID $CURSTARTROW ###############"

 # output file prefix 
OUTFILE=$CURFILE.treekin.p0-$CURSTARTID.t8-$MAXTIME
# check if output file exists already (do not replace)
if [ ! -f $OUTFILE.out.bz2 ]; then
 # create temporary subdir to avoid parallel file generation
 CURPWD=$PWD
 TMPDIR=$(mktemp -d -p $CURPWD)
#echo -e "\nINFO: creating treekin temporary directory $TMPDIR for $PREFIX";
 cd $TMPDIR  
 # ensure file naming for treekin call
 ln -s ../$CURFILE.rates rates.out;
 # call treekin
treekin -m I --p0 $CURSTARTROW=1  --t8=$MAXTIME < ../$CURFILE.barriers > $CURPWD/$OUTFILE.out;
 # go back to previous working directoy
 cd $CURPWD
 # cleanup temporary directory
 rm -rf $TMPDIR
 # compress treekin output 
 bzip2 $OUTFILE.out; 
 rm -rf $OUTFILE.out
 # generate output figure in pdf format using R
fi # treekin output exists
if [ $Ravailable == "1" ]; then
R --vanilla --silent -e "k <- read.table(\"$OUTFILE.out.bz2\", header=F, sep=\"\", comment.char=\"#\", skip=8);pdf(\"$OUTFILE.pdf\");matplot(k[,1], k[,2:ncol(k)], main=\"level = $CURLVL, p0 = $CURSTARTID\", xlab=\"time (arbitrary units)\", ylab=\"state probability\", ylim=c(0,1), log=\"x\", type=\"l\");dev.off(); q();"
fi # R available
}


# set start state for trajectory computation
STARTID=STARTSTATE
[[ $STARTSTATE == 0 ]] && STARTID=$OCID;
# ensure STARTSTATE is within mfe component
[[ $(grep -m 1 -c -P "^\\s*$STARTID\\s+" $PREFIX.barriers.out) == 1 ]] || die "\nERROR : start state $STARTID is not within mfe component\n";

# handle level 1
ln -s -f $PREFIX.barriers.out $PREFIX.barriers.out.1.barriers
ln -s -f $PREFIX.barriers.rates $PREFIX.barriers.out.1.rates
runTreekin $PREFIX.barriers.out.1 1 $STARTID $STARTID

# get maximal level
MAXLVL=$(ls $PREFIX.*.gradient | awk -F "." '{print $4}' | tail -n 1)
# reduce by one (last level has only one state)
MAXLVL=$((MAXLVL-1))

# iterate over each level
for LVL in $(seq 2 $MAXLVL); do
	CURPREFIX=$PREFIX.barriers.out.$LVL
	# get current STARTID via gradient mapping from last level (last column)
	LASTSTARTID=$STARTID
	STARTID=$(grep -m 1 -P "^\\s*$LASTSTARTID\\s+" $CURPREFIX.gradient | awk '{print $NF}')
	# get row of STARTID within barriers file (treekin input)
	STARTROW=$(awk -v ocid=$STARTID '$1==ocid{print (NR-1)}' $CURPREFIX.barriers)
	# run treekin
	runTreekin $CURPREFIX $LVL $STARTID $STARTROW
done # iterate all LVL 

