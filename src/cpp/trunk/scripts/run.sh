#!/bin/perl

# This script tries to identify the best machine among
# a set for running calculation

my @m = (qw/s0 s1 s2 s3 s4 s5 s6/);


foreach $m (@m) {
	$res = `rup $m`;

	if ($res =~ /load average: (\S*),/i) {
		$m{$m} = $1;
		print $1;
	} else {
		die "unable to get load remotely";
	}

}



