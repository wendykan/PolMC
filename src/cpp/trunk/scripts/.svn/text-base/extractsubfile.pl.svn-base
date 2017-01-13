#!/usr/bin/perl 

my $spitout = 0;

my ($begin, $end) = @ARGV;

while (<STDIN>) {

	if ( m|$begin|i) {
		$spitout = 1;
		next;		
	}
	if ( m|$end|i) {
		$spitout = 0;
		
	}

	if ($spitout) {
		print;
	}

}
