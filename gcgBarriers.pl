#!/usr/bin/env perl

###############################################################################
# Iteratively compresses an energy landscape encoded in barries output format
# using generalized coarse graining based on gradient-neighbor aggregation.
#
# author Martin Raden - 2018
###############################################################################

use strict;
use List::Util qw(sum max min);

# check argument number
if ( ($#ARGV + 1) < 2 ) {
    print "\nUsage: ".$0." <barriersOutputFile> <barriersRatesFile> [list of IDs not to be merged]\n";
    exit -1;
}

# get input files
my $barriersOutputFile = $ARGV[0];
my $barriersRatesFile = $ARGV[1];

my @basinsNotToMerge = ();
my $fileNamePrefix="";
if ( ($#ARGV + 1) > 2 ) {
    $fileNamePrefix .=".keep";
    for (my $i=2; $i <= $#ARGV; $i++) {
    	if (! $ARGV[$i] =~ /^\d+$/ or $ARGV[$i] < 1 ) {
    		print "\nArgument $i = '".$ARGV[$i]."' is no valid basin ID (>0) \n";
    		print "\nUsage: ".$0." <barriersOutputFile> <barriersRatesFile> [list IDs not to be merged]\n";
    		exit -1;
    	}
    	push ( @basinsNotToMerge, ($ARGV[$i])-1 );
	    $fileNamePrefix .="-".$ARGV[$i];
    }
}

use constant T => 310.15; # temperature in Kelvin

# check file existence
die "ERROR: can not find barriers output file '".$barriersOutputFile."'" unless -e $barriersOutputFile;
die "ERROR: can not find barriers rates file '".$barriersRatesFile."'" unless -e $barriersRatesFile;

# original local min index
my @locMinIdx = ();
# partition function array
my @Z = ();
# number of states array
my @states = ();
# define 1/R constant
use constant beta => 1/(8.314472 * 0.000239 * T); # in kcal/mol

print("##### reading partition functions from barriers output #####\n");

my @barriersOutput = ();
open( my $fBarOut, '<', $barriersOutputFile ) or die "ERROR: can not open barriers output file '".$barriersOutputFile."'";
my $sequence = <$fBarOut>; # ignore first line = sequence
while (my $row = <$fBarOut>) {
  $row =~ s/^\s+|\s+$//g;
  push (@barriersOutput,$row);
  my @cols = split(/\s+/,$row);
  push @locMinIdx, $cols[0];
  push @states,$cols[8];
  push @Z, ( exp(-$cols[9]*beta) );
}
close($fBarOut);
# print join(" ",@Z);

# check if basins not to merged are among the parsed basins
for (my $i=0; $i<=$#basinsNotToMerge; $i++) {
	if ($basinsNotToMerge[$i] > $#Z) {
    	print "\nBasin index '".$basinsNotToMerge[$i]."' is no valid basin ID\n";
    	print "\nUsage: ".$0." <barriersOutputFile> <barriersRatesFile> [list IDs not to be merged]\n";
    	exit -1;
	}
}

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


print("##### read rate matrix #####\n");
 
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


##############################################
# gradient neighbor for a given state index
# @param = ($from, \%rates, $dim, \@Z)
# @return = the index of the gradient neighbor or $from if no gradient neighbor exists
sub gradientNeighbor {
	my ($from,$rates,$dim,$Z) = @_;
	# identify maximal out-rate
	my $maxOutRate = -1;
	for (my $to=0; $to<$dim; $to++) {
		if ($to == $from) { next; }
		my $outRate = getRate($from,$to,$rates,$dim);
		# check for maximal rate
		if ($outRate>0 && $outRate>$maxOutRate) {
			# update maximal out rate
			$maxOutRate = $outRate;
		}
	}
	# find first adaptive neighbors with maximal out-rate
	for (my $to=0; $to<$dim; $to++) {
		if ($to == $from) { next; }
		my $outRate = getRate($from,$to,$rates,$dim);
		my $inRate = getRate($to, $from,$rates,$dim);
		# check if adaptive neighbor (use Z for final tie break, ie in- and outRate are equal)
		if ($outRate>0 && $outRate==$maxOutRate && $inRate>0 && $outRate>=$inRate && $Z->[$from] < $Z->[$to]) {
			# return first gradient-neighbor (keep barriers input order)
			return $to;
		}
	}
	# no gradient neighbor available
	return $from;
}

# correct diagonal
for (my $i=0; $i<$ratesDim; $i++) {
	setRate($i,$i,0,\%rates,$ratesDim);
}


###############################################
## gradient neighbor for a given state index
## @param = ($from, \%rates, $dim)
## @return = the index of the gradient neighbor or $from if no gradient neighbor exists
#sub gradientNeighbor {
#	my ($from,$rates,$dim) = @_;
#	# collect all adaptive neighbors
#	my @adaptive = ();
#	my @adaptiveOut = ();
#	my @adaptiveIn = ();
#	for (my $to=0; $to<$dim; $to++) {
#		if ($to == $from) { next; }
#		my $outRate = getRate($from,$to,$rates,$dim);
#		my $inRate = getRate($to, $from,$rates,$dim);
#		# check if adaptive neighbor
#		if ($outRate>0 && $inRate>0 && $outRate>=$inRate) {
#			# store adaptive neighbor
#			push (@adaptive, $to);
#			push (@adaptiveOut, $outRate);
#			push (@adaptiveIn, $inRate);
#		}
#	}
#	# check if no adaptive neighbors available
#	if (scalar(@adaptive)==0) {
#		return $from;
#	}
#	# get maximal out rate (for pre-gradient neighbor filtering)
#	my $maxOut = max(@adaptiveOut);
#	# find minimal back-rate among all pre-gradient neighbors
#	my $minIn = max(@adaptiveIn);
#	for (my $i=0; $i<scalar(@adaptive); $i++) {
#		# check if pre-gradient neighbor
#		if ($adaptiveOut[$i]==$maxOut) {
#			# update minimal back-rate for all gradient-neighbors
#			if ($minIn > $adaptiveIn[$i]) {
#				$minIn = $adaptiveIn[$i];
#			}
#		}
#	}
#	# return first gradient-neighbor (keep barriers input order)
#	for (my $i=0; $i<scalar(@adaptive); $i++) {
#		# check if gradient neighbor
#		if ($adaptiveOut[$i]==$maxOut && $adaptiveIn[$i]==$minIn) {
#			return $adaptive[$i];
#		}
#	}	
#}

######################  ITERATE FUNNEL GENERATION  ####################
my $abstractionLevel = 1;
do {

print "#states level $abstractionLevel = ".scalar(@Z)."\n";

	# increase abstraction level
	$abstractionLevel++;
	
	print("##### compile local funnel assignment for each state #####\n");
	my @funnel = ();
	# get gradient neighbor for each state
	for (my $from=0; $from<=$#Z; $from++) {
		push (@funnel, gradientNeighbor($from,\%rates,$ratesDim,\@Z));
	}
	#print join (" ",@funnel)."\n";
	
	# overwrite funnel assignment for basins not to be merged
	for (my $i=0; $i<=$#basinsNotToMerge; $i++) {
		$funnel[ $basinsNotToMerge[$i] ] = $basinsNotToMerge[$i];
	}
	#print join (" ",@funnel)."\n";
	
	print("##### find global funnel assignment for each state #####\n");
	for (my $from=0; $from<=$#funnel; $from++) {
		my $to = $funnel[$from];
		while ($to != $funnel[$to]) {
			$to = $funnel[$to];
		}
		$funnel[$from] = $to;
	}
	#print join (" ",@funnel)."\n";
	
	my @funnelCentersUnsrt = do { my %seen; grep { !$seen{$_}++ } @funnel };
	# sort funnel centers numerically
	my @funnelCenters = sort { $a <=> $b } @funnelCentersUnsrt;
	#print join(" ",@funnelCenters)."\n";
	
	# check if we get less funnels than we have basins so far
	if ($#funnelCenters == $#Z) {
		# no reduction -> stop
		exit 0;
	}
	
	print("##### store funnel assignment #####\n");
	open( $fBarOut, ">", "$barriersOutputFile.$abstractionLevel$fileNamePrefix.gradient" ) or die "ERROR: cannot open output file '$barriersOutputFile.$abstractionLevel$fileNamePrefix.gradient'";
	print $fBarOut $sequence;
	for (my $i=0; $i<=$#funnel; $i++) {
		print $fBarOut $barriersOutput[$i]." ".($locMinIdx[$funnel[$i]])."\n";
		#print $fBarOut $barriersOutput[$i]." ".($funnel[$i]+1)."\n";
	}
	close($fBarOut);

	# generate new data
	my @newZ = ();
	my @newStates = ();
	my @newBarriersOutput = ();
	my @newBasinsNotToMerge = ();
	my @newLocMinIdx = ();
	print("##### generate fake barriers output #####\n");
	open( $fBarOut, ">", "$barriersOutputFile.$abstractionLevel$fileNamePrefix.barriers" ) or die "ERROR: cannot open output file '$barriersOutputFile.$abstractionLevel$fileNamePrefix.barriers'";
	print $fBarOut $sequence;
	for (my $f=0; $f<=$#funnelCenters; $f++) {
		# get index of the funnel center
		my $i = $funnelCenters[$f];
		# collect partition function for funnel i
		my $Zi = 0;
		my $statesi = 0;
		for (my $j=0; $j<=$#Z; $j++) {
			if ($funnel[$j]==$i) {
				$Zi += $Z[$j];
				$statesi += $states[$j];
			}
		}
		push (@newZ, $Zi);
		push (@newLocMinIdx, $locMinIdx[$i]);
		push (@newStates, $statesi);
		push (@newBarriersOutput, $barriersOutput[$i]);
		# check if this funnel center is not to be merged
		if ( grep( /^$i$/, @basinsNotToMerge ) ) {
			# dont merge this funnel in next iteration
			push (@newBasinsNotToMerge, $#newZ);
		}
		# generate new output line with updated energy
		my @outputLine = split(/\s+/,$barriersOutput[$i]);
		# overwrite old data
		$outputLine[3] = 0;
		$outputLine[4] = "0.00";
		$outputLine[5] = 0;
		$outputLine[6] = 0;
		$outputLine[7] = 0;
		$outputLine[8] = $statesi;
		$outputLine[9] = -log($Zi) / beta;
		print $fBarOut join(" ",@outputLine)."\n";
	}
	close($fBarOut);
	
	print("##### generate new rates #####\n");
	my $newDim = $#funnelCenters +1;
	my %newRates;
	for( my $i=0; $i+1<$newDim; $i++) {
		my $f1 = $funnelCenters[$i];
		for ( my $j=$i+1; $j<$newDim; $j++) {
			my $f2 = $funnelCenters[$j];
			my $ratef1f2 = 0;
			my $ratef2f1 = 0;
			# collect all rates
			for (my $k=0; $k+1<=$#Z; $k++) {
				# check if in funnels of interest
				if ($funnel[$k]!=$f1 && $funnel[$k]!=$f2) {
					next;
				} 
				for (my $l=$k+1; $l<=$#Z; $l++) {
					# check if in funnels of interest
					if ($funnel[$l]!=$f1 && $funnel[$l]!=$f2) {
						next;
					}
					# check if in different funnels
					if ($funnel[$l]==$funnel[$k]) {
						next;
					}
					# update rates
					if ($funnel[$k]==$f1) {
						$ratef1f2 += $Z[$k] * getRate($k,$l,\%rates,$ratesDim);
						$ratef2f1 += $Z[$l] * getRate($l,$k,\%rates,$ratesDim);
					} else {
						$ratef1f2 += $Z[$l] * getRate($l,$k,\%rates,$ratesDim);
						$ratef2f1 += $Z[$k] * getRate($k,$l,\%rates,$ratesDim);
					}
				}
			}
			# normalize rates to get probability multiplier
			if ($ratef1f2 > 0) {
				$ratef1f2 /= $newZ[$i];
				$ratef2f1 /= $newZ[$j];
			}
			# store rate
			setRate($i,$j,$ratef1f2,\%newRates,$newDim);
			setRate($j,$i,$ratef2f1,\%newRates,$newDim);
		}
	}
	
	print("##### print new rates #####\n");
	open( $fBarRate, ">", "$barriersOutputFile.$abstractionLevel$fileNamePrefix.rates" ) or die "ERROR: cannot open output file '$barriersOutputFile.$abstractionLevel$fileNamePrefix.rates'";
	for( my $i=0; $i<$newDim; $i++) {
		for ( my $j=0; $j<$newDim; $j++) {
			my $curRate = getRate($i,$j,\%newRates,$newDim);
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
	
	
	# copy all data for next iteration
	@Z = @newZ;
	@states = @newStates;
	@barriersOutput = @newBarriersOutput;
	%rates = %newRates;
	$ratesDim = $newDim;
	@basinsNotToMerge = @newBasinsNotToMerge;
	@locMinIdx = @newLocMinIdx;
	
} while ( scalar(@Z) > 1 );

if (scalar(@Z) == 1) {
	print "#states level $abstractionLevel = ".scalar(@Z)."\n";
}
