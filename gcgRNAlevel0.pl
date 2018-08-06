#!/usr/bin/env perl

###############################################################################
# Computes the number of level 0 states for a given RNA sequence using a
# dynamic programming scheme as introduced by Waterman and Smith (1970).
#
# author Martin Raden - 2018
###############################################################################

use strict;

# check argument number
if ( ($#ARGV + 1) != 1 ) {
    print "\nUsage: ".$0." <RNASEQUENCE>\n";
    exit;
}

# get input files
my $seq = $ARGV[0];
my $n = length($seq);

# minimal loop length
use constant MINLOOP => 3;

# check sequence
die "ERROR: no RNA sequence {ACGU}* given : '".$seq."'" unless $seq =~ m/^[ACGU]+$/;
# convert to integers
$seq =~ tr/ACGU/0123/;
my @nt = split("", $seq);

# base pairing of nts
my @bp = (
# 	A	C	G	U
#A
	0,	0,	0,	1,
#C
	0,	0,	1,	0,
#G
	0,	1,	0,	1,
#U
	1,	0,	1,	0,
);


########### init counting matrix  ############
my @C = (1) x ($n * $n);

##############################################
# value of cell C[i,j]
# @param = ($i,$j)
sub getC {
	my ($i,$j) = @_;
	if ($i>$j || $i > $n || $i < 0 || $j > $n || $j < 0) {
		return 1;
	}
	return $C[ ($i * $n) + $j ];
}

##############################################
# set value of cell C[i,j]
# @param = ($from, $to, $value)
sub setC {
	my ($i,$j,$value) = @_;
	die "ERROR: from $i out of range [0,".($n-1)."]" if $i > $n || $i < 0;
	die "ERROR: from $j out of range [0,".($n-1)."]" if $j > $n || $j < 0;
	$C[ ($i * $n) + $j ] = $value;
}

# fill DP matrix
for (my $j = MINLOOP +1; $j<$n; $j++) {
	for (my $i = $j - MINLOOP - 1; $i>=0; $i--) {
		# j is unpaired
		my $value = getC($i,$j-1);
		# j is paired cases
		for (my $k=$i; $k+MINLOOP<$j; $k++) {
			if ($bp[$nt[$k]*4 + $nt[$j]] != 0) {
				$value += getC($i,$k-1) * getC($k+1,$j-1) * 1;
			}
		}
		setC($i,$j,$value);
	}
}

# print matrix
#for (my $i=0; $i<$n; $i++) { for (my $j=0; $j<$n; $j++) { printf " ".getC($i,$j); }	print "\n"; }

# print overall structure number
print getC(0,$n-1)."\n";

