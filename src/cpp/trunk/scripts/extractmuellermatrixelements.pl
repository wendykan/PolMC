#!/usr/local/bin/perl 

my ($S1000, $S1100, $S1010, $S1001) = @ARGV;

`cp $S1000.I.dat M11.dat`;
SubtractFiles("$S1100.I.dat","$S1000.I.dat","M12.dat");
SubtractFiles("$S1010.I.dat","$S1000.I.dat","M13.dat");
SubtractFiles("$S1001.I.dat","$S1000.I.dat","M14.dat");

`cp $S1000.Q.dat M21.dat`;
SubtractFiles("$S1100.Q.dat","$S1000.Q.dat","M22.dat");
SubtractFiles("$S1010.Q.dat","$S1000.Q.dat","M23.dat");
SubtractFiles("$S1001.Q.dat","$S1000.Q.dat","M24.dat");

`cp $S1000.U.dat M31.dat`;
SubtractFiles("$S1100.U.dat","$S1000.U.dat","M32.dat");
SubtractFiles("$S1010.U.dat","$S1000.U.dat","M33.dat");
SubtractFiles("$S1001.U.dat","$S1000.U.dat","M34.dat");

`cp $S1000.V.dat M41.dat`;
SubtractFiles("$S1100.V.dat","$S1000.V.dat","M42.dat");
SubtractFiles("$S1010.V.dat","$S1000.V.dat","M43.dat");
SubtractFiles("$S1001.V.dat","$S1000.V.dat","M44.dat");

sub SubtractFiles()
{
	print "@_\n";

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
