#! /usr/bin/perl 


use Getopt::Std;

use vars qw($opt_s $opt_r $opt_n);

# Options for acceptance [w]idth and [h]eaders, a[n]gle

getopts('hw:n:');


if ($opt_h ) {
	$header = 1;
}

if ($opt_w ) {
	$acceptance_width = $opt_w;
} else {
	$acceptance_width = 0;
	# print "Acceptance width set to 0\n";
}


if ($opt_n){
	$acceptanceCosine = $opt_n;
} else {
	$acceptanceCosine = 0
}

my $sum = 0;


####
# Get param.dat parameters
####

my $mu_s;
my $height;
my $width;
my $rotationPerCmInClearSpace;

my $s = `grep "mu_s" param.dat`;

if ( $s =~ m|mu_s=(\d.*?)$|g ) {
	$mu_s = $1;
}

my $s = `grep "mu_a" param.dat`;

if ( $s =~ m|mu_a=(\d.*?)$|g ) {
	$mu_a = $1;
}

$s = `grep "yheight" param.dat`;

if ( $s =~ m|yheight=(.+?);?$|g ) {
	$height = $1;
}

$s = `grep "xwidth" param.dat`;

if ( $s =~ m|xwidth=(.+?);?$|g ) {
	$width = $1;
}

$s = `grep "zdepth" param.dat`;

if ( $s =~ m|zdepth=(.+?);?$|g ) {
	$depth = $1;
}

$s = `grep "rotationPerCmInClearSpace" param.dat`;

if ( $s =~ m|^rotationPerCmInClearSpace=(.+?);?$|g ) {
	$rotationPerCmInClearSpace = $1;
} else {
	$rotationPerCmInClearSpace = 0;
}

$s = `grep "g" polmc.log`;

if ( $s =~ m|g from matrix = (\d.*)|g ) {
	$anisotropy = $1;
} else {
	$anisotropy = 0;
}

push @all_titles, ("mu_s","mu_a","g","xwidth","yheight","zdepth","image_width","rotationPerCmInClearSpace");
push @all_data, ($mu_s,$mu_a,$anisotropy,$width,$height,$depth,$acceptance_width, $rotationPerCmInClearSpace);

####
# Get S1100 I info
####

($t_ref,$d_ref) = GetStats("OutS1100.dat.forward.I.$acceptanceCosine.dat","S1100_I_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1100 Q info
####

($t_ref,$d_ref) = GetStats("OutS1100.dat.forward.Q.$acceptanceCosine.dat","S1100_Q_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1100 U info
####

($t_ref,$d_ref) = GetStats("OutS1100.dat.forward.U.$acceptanceCosine.dat","S1100_U_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1100 avgpath info
####

($t_ref,$d_ref) = GetStats("OutS1100.dat.forward.avgpath.$acceptanceCosine.dat","S1100_U_avgpath");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1100 numscatter info
####

($t_ref,$d_ref) = GetStats("OutS1100.dat.forward.numscatter.$acceptanceCosine.dat","S1100_numscatter");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;



####
# Get S1010 I info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.I.$acceptanceCosine.dat","S1010_I_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 Q info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.Q.$acceptanceCosine.dat","S1010_Q_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 U info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.U.$acceptanceCosine.dat","S1010_U_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 betalin info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.betalin.$acceptanceCosine.dat","S1010_betalin_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 avgdist info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.avgpath.$acceptanceCosine.dat","S1010_avgpath_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1-100 I info
####

($t_ref,$d_ref) = GetStats("OutS1-100.dat.forward.I.$acceptanceCosine.dat","S1-100_I_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1-100 Q info
####

($t_ref,$d_ref) = GetStats("OutS1-100.dat.forward.Q.$acceptanceCosine.dat","S1-100_Q_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1-100 U info
####

($t_ref,$d_ref) = GetStats("OutS1-100.dat.forward.U.$acceptanceCosine.dat","S1-100_U_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;


####
# Get m11 info
####

($t_ref,$d_ref) = GetStats("M11.$acceptanceCosine.dat","M11_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;


####
# Get m12 info
####

($t_ref,$d_ref) = GetStats("M12.$acceptanceCosine.dat","M12_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get m21 info
####

($t_ref,$d_ref) = GetStats("M21.$acceptanceCosine.dat","M21_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get m23 info
####

($t_ref,$d_ref) = GetStats("M23.$acceptanceCosine.dat","M23_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get m32 info
####

($t_ref,$d_ref) = GetStats("M32.$acceptanceCosine.dat","M32_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get m22 info
####

($t_ref,$d_ref) = GetStats("M22.$acceptanceCosine.dat","M22_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get m33 info
####

($t_ref,$d_ref) = GetStats("M33.$acceptanceCosine.dat","M33_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 avgpath info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.avgpath.$acceptanceCosine.dat","S1010_avgpath_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 avgpath2 info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.avgpath2.$acceptanceCosine.dat","S1010_avgpath2_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

####
# Get S1010 avgpath3 info
####

($t_ref,$d_ref) = GetStats("OutS1010.dat.forward.avgpath3.$acceptanceCosine.dat","S1010_avgpath3_");

push @all_titles, @$t_ref;
push @all_data, @$d_ref;

###
# acceptance angle
###

$s = `grep "acceptanceCos" param.dat`;

if ( $s =~ m|acceptanceCos=(.+?);?$|g ) {
	$acceptanceCos = $1;
}

push @all_titles, ("acceptance");
push @all_data, ($acceptanceCos);


####
# Print all
####

if ($header) {
	print "#";
	print (join "\t",@all_titles);
	print "\n";
	}
print (join "\t",@all_data),"\n";
print "\n";


sub GetStats() 
{
	my ($file, $prefix) = @_;

	my $n_y = 0;
	my $n_x = 0;
	my   @all = ();
	
	if ( ! -e $file ) {
		return
	}
	
	unless (open INPUT, $file) {
		die "Error with $file in ".`pwd`.": $!";
	}
		
	while(<INPUT>) {
	
		@row = split "\t";
	
		if (scalar @row) {
			# If row is not empty, count it
			$n_x++;
			
			if ($n_y && scalar @row != $n_y) {
				die "Error: number of columns inconsistent";
			}
			
			$n_y = scalar @row;
		}
		
		push @all, [ @row ];
	
	}


	if ($n_x == 0 || $n_y == 0) {
		return;
	}
		
	my $weight = 0;
	my $sum = 0;
	my $normalize = 0;
	my $dx = $width / $n_x;
	my $dy = $height / $n_y;

	for ($i = 0; $i < $n_x; $i++) {
		for ($j = 0; $j < $n_y; $j++) {

			# From polmc/main.cpp:
			# long((inValue-inMin)/abs(inMax - inMin)* inN)
			
			# i = long((pos.x / width + 1.) * Nx / 2 );
            # j = long((pos.y / height + 1. ) * Ny / 2 );
			# we invert to get:

			if ($acceptance_width == 0) {
				if (($i == int $n_x/2) && ($j == int $n_y/2) ) {
					$weight = 1;
				} else {
					$weight = 0;
				}
			} else {
				$x = $dx * ($i - int $n_x/2);
				$y = $dy * ($j - int $n_y/2);
#				print "x: $x, y: $y is $i, $j\n";
				
				$weight = exp( - ($x*$x + $y*$y )/ ($acceptance_width*$acceptance_width));
			}
			$sum +=  $weight * $all[$i][$j];
			$normalize += $weight
			
		}
	 }

	$integral = $sum ;

	@titles=($prefix."center",$prefix."integr");
	@data=($all[int $n_x/2][int $n_y/2],$integral );
	
	return ([@titles],[@data]);
	
}
