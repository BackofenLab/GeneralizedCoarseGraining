#!/usr/bin/env bash

MAXMEM=4G

# maximal time for kinetics computation
MAXTIME=100000000000; # maximal time (step) to compute

# exit function
die() { echo -e "$*" 1>&2 ; exit 1; }

# get RNA sequence
RNA=""
PREFIX=$RNA
if [ "$#" -gt "0" ]; then
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
fi
# check command line arguments
[[ ($# -eq 1 || ($# -eq 2 && $2 =~ ^[0-9]+$ )) && ($RNA =~ ^[ACGU]+$ ) ]] || die "\n usage : gcgRNAkinetics.sh <RNA sequence (file)> <MAXTIME=$MAXTIME>";
# (optional) get maximal time for kinetics computation 
[[ $# -eq 2 ]] && MAXTIME=$2

# whether or not R is available for plotting
Ravailable=1;
# ensure required binaries are in PATH
[[ $(type -P "RNAsubopt") ]] || die "\nERROR : RNAsubopt not in PATH\n";
[[ $(type -P "barriers") ]] || die "\nERROR : barriers not in PATH\n";
[[ $(type -P "treekin") ]] || die "\nERROR : treekin not in PATH\n";
[[ $(type -P "gcgBarriers.pl") ]] || die "\nERROR : gcgBarriers.pl not in PATH or not executable\n";
if [ ! $(type -P "R") ]; then
	echo -e "\nWARNING : R not in PATH => skipping figure generation\n"; 
	Ravailable=0;
fi

##############  LEVEL 0 ENUMERATION  #################

if [ ! -f $PREFIX.RNAsubopt.zip ]; then

# enumerate (and sort) all secondary structures (first by energy using structure string for tie breaking)
echo $RNA | RNAsubopt --deltaEnergy=99999 | sort -k 2,2n -k 1,1dr -S $MAXMEM | zip $PREFIX.RNAsubopt.zip -;

fi

##############  LEVEL 1 COARSE GRAINING  #################

if [ ! -f $PREFIX.barriers.out ]; then

# create temporary subdir to avoid parallel file generation
CURPWD=$PWD
TMPDIR=$(mktemp -d -p $CURPWD)
echo -e "\nINFO: creating barriers temporary directory $TMPDIR for $PREFIX";
cd $TMPDIR  
# run barriers to compute coarse graining level 1
unzip -p $CURPWD/$PREFIX.RNAsubopt.zip | barriers --rates -G RNA -M noShift --bsize --max=999999 --minh=0 > $CURPWD/$PREFIX.barriers.out;
# store rate matrix generated by barriers
mv rates.out $CURPWD/$PREFIX.barriers.rates;
# go back to previous working directoy
cd $CURPWD
# cleanup obsolete barriers files and temporary directory
rm -rf $TMPDIR

fi

##############  LEVEL >1 COARSE GRAINING  #################

if [ ! -f $PREFIX.gcgBarriers.out ]; then

(
# get micro-state landscape size
printf "#states level 0 = "; unzip -p $PREFIX.RNAsubopt.zip | grep -c -v $RNA;
# compute generalized coarse grainings
gcgBarriers.pl $PREFIX.barriers.out $PREFIX.barriers.rates
) > $PREFIX.gcgBarriers.out

fi

################  PLOTTING  ######################


# define variables for the following shell function
# shell function (tested in bash) to compute one plot via treekin 
function runTreekin {
 # output file prefix 
 OUTFILE=$FILE.treekin.p0-$OCID.t8-$MAXTIME
# check if output file exists already (do not replace)
if [ ! -f $OUTFILE.out.bz2 ]; then
	CURFILE=$1 # base file name without extension (.barriers .rates)
	CURLVL=$2  # funnel of open chain state to be taken from according barriers output file
	CUROCID=$3 # the coarse graining level
 # create temporary subdir to avoid parallel file generation
 CURPWD=$PWD
 TMPDIR=$(mktemp -d -p $CURPWD)
 echo -e "\nINFO: creating treekin temporary directory $TMPDIR for $PREFIX";
 cd $TMPDIR  
 # ensure file naming for treekin call
 ln -s ../$CURFILE.rates rates.out;
 # call treekin
treekin -m I --p0 $CUROCID=1  --t8=$MAXTIME < ../$CURFILE.barriers > $CURPWD/$OUTFILE.out;
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
R --vanilla --silent -e "k <- read.table(\"$OUTFILE.out.bz2\", header=F, sep=\"\", comment.char=\"#\", skip=8);pdf(\"$OUTFILE.pdf\");matplot(k[,1], k[,2:ncol(k)], main=\"level = $CURLVL, p0 = $CUROCID\", xlab=\"time (arbitrary units)\", ylab=\"state probability\", ylim=c(0,1), log=\"x\", type=\"l\");dev.off(); q();"
fi # R available
}


# handle level 1
ln -s -f $PREFIX.barriers.out $PREFIX.barriers.out.1.barriers
ln -s -f $PREFIX.barriers.rates $PREFIX.barriers.out.1.rates
# get open chain ID
OCID=$(grep -P "^\\s*\\d+\\s+[\\.]+\\s+0" $PREFIX.barriers.out | awk 'NR==1 {print $1}');
runTreekin $PREFIX.barriers.out.1 1 $OCID

# get maximal level
MAXLVL=$(ls $PREFIX.*.gradient | awk -F "." '{print $4}' | tail -n 1)
# reduce by one (last level has only one state)
MAXLVL=$((MAXLVL-1))

# iterate over each level
for LVL in $(seq 2 $MAXLVL); do
	# get current OCID via gradient mapping from last level (last column)
	LASTOCID=$((OCID+1))
	OCID=$(grep -m 1 -P "^\\s*$LASTOCID\\+" $FILE.gradient | awk '{print $NF}')
	# run treekin
	runTreekin $PREFIX.barriers.out.$LVL $LVL $OCID
done # iterate all LVL 

