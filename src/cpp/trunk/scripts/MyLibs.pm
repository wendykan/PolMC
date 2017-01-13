package MyLibs;
require Exporter;
@ISA = qw (Exporter);
@EXPORT = qw (GetDirectory GetLoadAround UpdateChildStatus SubtractFiles AddFiles MultiplyFiles DivideFiles);
use POSIX "sys_wait_h";

sub GetDirectory {
    my ($prefix, $i) = @_;
    my $dir;
    
    if ($i >= 100) {
        $dir = $prefix."_$i";
    } elsif ($i >= 10) {
        $dir = $prefix."_0$i";
    } else {
        $dir = $prefix."_00$i";
    }

    print "$dir\n";
    return $dir;
}

sub UpdateChildStatus {
    my ($host, @childs) = @_;

    for ($i = scalar @childs; $i > 0 ; $i--) {
        $pid = pop @childs;
        if ( ($err = waitpid $pid,WNOHANG) == 0) {
            # Not finished, put back in list
            unshift @childs, $pid;
        } else {
            # Done or killed. Let's move on.        
        }
    
    }
    return @childs;
}

sub GetLoadAround {
    my %machines = ( "s0" => 0.9,
                     "s1" => 0.9,
                     "s2" => 0.9,
                     "s3" => 0.9,
                     "s4" => 0.9,
                     "s5" => 0.9,
                     "s6" => 0.9,
                     "nebula" => 0.9);
                     
    foreach $machine (sort keys %machines) {
        $uptime =  `rsh $machine uptime`;
    
        if ($uptime =~ m|(\d.\d\d)$|i) {
            $load{$machine} = $1;
        }
    }

    return %load;
}


sub AddFiles
{
	my ($outfile, @files) = @_;
	 
	my $hasBeenInitialized = 0;
	my @total;
	
	while ($file = shift @files) {
		unless (open FILE, $file) {
			die "Error with $file: $!";
		}
		
		$i = 0;
		while($line =  <FILE>) {
			chomp($line);
			if ($hasBeenInitialized){
				my @row = split "\t", $line;
				if ($i > scalar @total) {
					die "Not same number of columns:", $i,scalar @total;
				}
				
				foreach $j (0..@row-1) {
					$total[$i][$j] += $row[$j];
					if ($j > scalar @{$total[$i]}) {
						die "Not same number of columns:", scalar @row," and ", scalar $$total[$i];
					}
				}		
			} else {
				my @row = split "\t", $line;
				push  @total, \@row;
			}
			$i++;
		}
		
		close (FILE);
		$hasBeenInitialized = 1;
	}	

		unless (open OUT, ">$outfile"){
			die "Error with $outfile:$!\n";
		}
		
		foreach $i (0..@total-1) {
			$t = $total[$i];
			print  OUT join "\t",@$t,"\n";
		}
		close (OUT);
}


sub SubtractFiles
{
	my ($file1, $file2, $outfile) = @_;
	 
	unless (open FILE1, $file1) {
		die "Error with $file1: $!";
	}
	
	unless (open FILE2, $file2) {
		die "Error with $file2: $!";
	}
	
	unless (open OUT, ">$outfile") {
		die "Error with $outfile: $!";
	}

	while($line1 =  <FILE1>) {
		chomp($line1);
		$line2 = <FILE2>;
		chomp($line2);


		my (@row1, @row);

		@row1 = split "\t", $line1;
		@row2 = split "\t", $line2;

		if (scalar @row1 != scalar @row2) {
			die "Not same number of columns:", scalar @row1,scalar @row2;
		}
		
		foreach $i (0..@row1-1) {
			$sub[$i] = $row1[$i]-$row2[$i];
		}		
		print OUT join("\t",@sub)."\n";
	}
	
	close (FILE1);
	close (FILE2);
	close (OUT);
	
}

sub DivideFiles
{
        my ($file1, $file2, $outfile) = @_;

		
        unless (open FILE1, $file1) {
                die "Error with $file1: $!";
        }

        unless (open FILE2, $file2) {
                die "Error with $file2: $!";
        }

        unless (open OUT, ">$outfile") {
                die "Error with $outfile: $!";
        }

        while($line1 =  <FILE1>) {
				chomp($line1);
				$line2 = <FILE2>;
				chomp($line2);

				my (@row1, @row);

                @row1 = split "\t", $line1;
                @row2 = split "\t", $line2;

				if (scalar @row1 != scalar @row2) {
					die "Not same number of columns:", scalar @row1,scalar @row2;
				}
				
                foreach $i (0..@row1-1) {
					if ($row2[$i] != 0) {
						$div[$i] = $row1[$i]/$row2[$i];
					} else {
						$div[$i] = 0;
					}
                }
                print OUT join("\t",@div)."\n";
        }

        close (FILE1);
        close (FILE2);
        close (OUT);
}


1;
