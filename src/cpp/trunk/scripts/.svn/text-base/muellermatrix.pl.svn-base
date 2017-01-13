#!/usr/local/bin/perl 

@I = (qw / Iunpol.dat Ipol.dat Ipol45.dat Ipolcirc.dat /);
@Q = (qw / Qunpol.dat Qpol.dat Qpol45.dat Qpolcirc.dat /);
@U = (qw / Uunpol.dat Upol.dat Upol45.dat Upolcirc.dat /);
@V = (qw / Vunpol.dat Vpol.dat Vpol45.dat Vpolcirc.dat /);

SubtractFiles("Ipol.dat","Iunpol.dat","M12.dat");
SubtractFiles("Ipol45.dat","Iunpol.dat","M13.dat");
SubtractFiles("Ipolcirc.dat","Iunpol.dat","M14.dat");

SubtractFiles("Qpol.dat","Qunpol.dat","M22.dat");
SubtractFiles("Qpol45.dat","Qunpol.dat","M23.dat");
SubtractFiles("Qpolcirc.dat","Qunpol.dat","M24.dat");

SubtractFiles("Upol.dat","Uunpol.dat","M32.dat");
SubtractFiles("Upol45.dat","Uunpol.dat","M33.dat");
SubtractFiles("Upolcirc.dat","Uunpol.dat","M34.dat");

SubtractFiles("Vpol.dat","Vunpol.dat","M42.dat");
SubtractFiles("Vpol45.dat","Vunpol.dat","M43.dat");
SubtractFiles("Vpolcirc.dat","Vunpol.dat","M44.dat");

DivideFiles("Qpol.dat","Ipol.dat","beta22.dat");
DivideFiles("Upol45.dat","Ipol45.dat","beta33.dat");
DivideFiles("Vpolcirc.dat","Ipolcirc.dat","beta44.dat");



sub SubtractFiles()
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
		$line2 = <FILE2>;

		@row1 = split "\t", $line1;
		@row2 = split "\t", $line2;
		
		foreach $i (1..@row1) {
			$sub[$i] = $row1[$i]-$row2[$i];
		}		
		print OUT join("\t",@sub)."\n";
	}
	
	close (FILE1);
	close (FILE2);
	close (OUT);
	
}

sub DivideFiles()
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
                $line2 = <FILE2>;

                @row1 = split "\t", $line1;
                @row2 = split "\t", $line2;

                foreach $i (1..@row1) {
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
