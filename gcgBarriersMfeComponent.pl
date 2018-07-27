#!/usr/bin/env perl

###############################################################################
# Identifies the mfe-(state 1)-containing connected component of the given 
# gradient-basing energy landscape (as produced by barriers) and prints a 
# reduced pseudo-barriers output and rates matrix.
#
# This is needed to ensure ergodicity of the rate matrix as assumed by treekin. 
#
# author Martin Raden - 2018
###############################################################################

use strict;
use List::Util qw(sum max min);

# check argument number
if ( ($#ARGV + 1) != 2 ) {
    print "\nUsage: ".$0." <barriersOutputFile> <barriersRatesFile>\n";
    exit -1;
}

# get input files
my $barriersOutputFile = $ARGV[0];
my $barriersRatesFile = $ARGV[1];

my $fileNamePrefix="";

# check file existence
die "ERROR: can not find barriers output file '".$barriersOutputFile."'" unless -e $barriersOutputFile;
die "ERROR: can not find barriers reates file '".$barriersRatesFile."'" unless -e $barriersRatesFile;

# original local min index
my @locMinIdx = ();

# read partition functions from barriers output
my @barriersOutput = ();
open( my $fBarOut, '<', $barriersOutputFile ) or die "ERROR: can not open barriers output file '".$barriersOutputFile."'";
my $sequence = <$fBarOut>; # ignore first line = sequence
while (my $row = <$fBarOut>) {
  $row =~ s/^\s+|\s+$//g;
  push (@barriersOutput,$row);
  my @cols = split(/\s+/,$row);
  push @locMinIdx, $cols[0];
}
close($fBarOut);


##############################################
# generates the hash key for a given index pair
# @param = ($from, $to, $dimension)
# @return = hash key
sub getKey {
	my ($from,$to,$dim) = @_;
	die "ERROR: from $from out of range [0,".($dim-1)."]" if $from > $dim || $from < 0;
	die "ERROR: to $to out of range [0,".($dim-1)."]" if $to > $dim || $to < 0;
	return ($from * $dim) + $to;
}

##############################################
# rate for a given transition (by indices)
# @param = ($from, $to, \%rates, $dimension)
# @return = rate for given transition
sub getRate {
	my ($from,$to,$rates,$dim) = @_;
	my $key = getKey($from,$to,$dim);
	# get rate if something is stored
	if ( exists $rates->{ $key } ) {
		return $rates->{ $key };
	} else {
		# otherwise zero
		return 0.0;
	}
}

##############################################
# set rate for a given transition (by indices)
# @param = ($from, $to, $rate, \%rates, $dim)
# @return = updates @rates
sub setRate {
	my ($from,$to,$rate,$rates,$dim) = @_;
	# store rate if non-zero
	if ($rate > 0) {
		$rates->{ getKey($from,$to,$dim) } = $rate;
	}
}


# read rate matrix 
open( my $fBarRate, '<', $barriersRatesFile ) or die "ERROR: can not open barriers rates file '".$barriersRatesFile."'";
my %rates;
my $ratesDim=0;
my $from=0;
while (my $row = <$fBarRate>) {
  $row =~ s/^\s+|\s+$//g;
  my @cols = split(/\s+/,$row);
  $ratesDim=($#cols)+1;
  for (my $to=0; $to<$ratesDim; $to++) {
  	setRate( $from, $to, $cols[$to], \%rates, $ratesDim );
  }
  $from++;
}
close($fBarRate);

# correct diagonal
for (my $i=0; $i<$ratesDim; $i++) {
	setRate($i,$i,0,\%rates,$ratesDim);
}



######################  IDENTIFY MFE COMPONENT  ####################

my @mfeCompIds = (); # ids of all loc mins of the mfe component
my $mfeCompIdString = ","; # concatenation of all mfeCompIds for regex lookup
my @mfeCompIdsToCheck = (0); # ids of loc mins still to process

# check if something left for processing
while( scalar(@mfeCompIdsToCheck) > 0 ) {
	
	# get element to be processed
	my $from = pop @mfeCompIdsToCheck;
	
	# skip if already processed
	if ($mfeCompIdString =~ /,$from,/) {
		print "skipping ".($from+1)." (already processed)\n";
		next;
	}
	print "processing ".($from+1)."\n";
	
	# store as processed
	push @mfeCompIds, $from;
	$mfeCompIdString .= $from.",";
	
	# check all neighbors if already part of the component
	for (my $to=0; $to<=$#locMinIdx; $to++) {
		if ($to == $from) { next; }
		# check if from and to are neighbored
		my $outRate = getRate($from,$to,\%rates,$ratesDim);
		if ($outRate > 0) {
			# check if to was already processed
			if ( $mfeCompIdString !~ /,$to,/ ) {
				push @mfeCompIdsToCheck, $to;	
			}
		}
	}
	
} # end : while not all processed 

###########  WRITE MFE COMPONENT  #################

# sort ids of mfe component
my @mfeCompIdsSorted = sort { $a <=> $b } @mfeCompIds;

# generate reduced barriers output
open( $fBarOut, ">", "$barriersOutputFile.mfeComp.barriers" ) or die "ERROR: cannot open output file '$barriersOutputFile.mfeComp.barriers'";
print $fBarOut $sequence;
for (my $f=0; $f<scalar(@mfeCompIdsSorted); $f++) {
	# print according line
	print $fBarOut $barriersOutput[$mfeCompIdsSorted[$f]]."\n";
}
close($fBarOut);

# print rates between mfe component ids
my $newDim = scalar(@mfeCompIdsSorted);
open( $fBarRate, ">", "$barriersOutputFile.mfeComp.rates" ) or die "ERROR: cannot open output file '$barriersOutputFile.mfeComp.rates'";
for( my $i=0; $i<$newDim; $i++) {
	for ( my $j=0; $j<$newDim; $j++) {
		# get rate between according minima of sorted mfe comp ID set
		my $curRate = getRate($mfeCompIdsSorted[$i],$mfeCompIdsSorted[$j],\%rates,$ratesDim);
		if ($curRate == 0) {
			print $fBarRate sprintf(" %9s ",$curRate)
		} elsif ($curRate > 0.0001) {
			print $fBarRate sprintf(" %.7f ",$curRate)
		} else {
			print $fBarRate sprintf(" %.3e ",$curRate);
		}
	}
	print $fBarRate "\n";
}
close($fBarRate);


