#! /usr/bin/perl -w

use Getopt::Std;
use POSIX "sys_wait_h";
use MyLibs;


sub usage {
    print "\nThis program is a utility to setup, launch and analyze multiple Monte Carlo runs\n";
    print "It is not required to run Monte Carlo calculations, but might make your life easier\n";
    print "when running multiple calculations. However, it might be a little specific to the original\n";
    print "network environment for which it was written (multiple CPU cluster).\n";
    print "Typically, you provide a base directory and a few parameters and let the program cycle through all of them.\n";
    print "It can run in three modes: setup, run and analyze.  These modes and their options are:\n\n";
    
    print "-s      : setup\n";
    print "-p[dir] : prefix directory is dir\n";
    print "-f[n]   : first suffix (default n=1)\n";
    print "-l[n]   : last suffix n\n";
    print "-b[n]   : base suffix n from which param.dat is copied\n";
    print "-m[n]   : repeat same calculation n times, create subdirectories\n";
    print "-o      : will overwrite any existing files during setup (default no)\n";
    
    print "\n\n";
    print "-r      : run calculations\n";
    print "-p[dir] : prefix directory is dir\n";
    print "-f[n]   : first suffix (default n=1)\n";
    print "-l[n]   : last suffix n\n";
    print "-b[n]   : base suffix n from which param.dat is copied\n";

    print "\n\n";


}
 
# Options for [s]etup, [r]un [a]nalyse, [p]refix [b]asesuffix [f]irstsuffix [l]astsuffix
# and acceptance cosine a[n]gle [t]hreads analy[z]issuffix [c]lean [d]istclean
#-[m]ultiple submissions in directories
# 
# When setting up directories: -m[u]_s increment -optical activity[g] increment [o]verwrite
# 

use vars qw($opt_s $opt_r $opt_a $opt_c $opt_D $opt_O $opt_d $opt_o);

getopts('sra:p:b:f:l:wn:t:cdz:u:og:m:');


my $howmany;
my $prefix;
my $basedir;
my $polltime = 5;

if ($opt_c ) {
	$clean = "clean";
}

if ($opt_d ) {
	$clean = "distclean";
}

if ($opt_s ) {
	$setup = 1;
}

if ($opt_r ) {
	$run = 1;
}

if ($opt_a ) {
	$analyse = $opt_a;
} else {
	$analyse = 0;
}

if ($opt_p ) {
	$prefix = $opt_p;
} else {
    usage();
	die "You must specify a directory prefix (-p dir)\n";
}

if ($opt_b ) {
	$basesuffix = $opt_b;
} else {
	$basesuffix = 1;
	print "Base suffix ${prefix}_001 assumed\n";
}
$basedir = GetDirectory($prefix, $basesuffix);

if ($opt_f ) {
	$firstsuffix = $opt_f;
} else {
    usage();
	die "You must specify the first directory suffix (-f num)\n";
}

if ($opt_l ) {
	$lastsuffix = $opt_l;
} else {
	die "You must specify the last directory suffix (-l num)\n";
}

if ($opt_z ) {
	$analyzesuffix = $opt_z;
	print "$analyzesuffix used for analysis directory\n";
} else {
	$analyzesuffix = "_analysis";
	print "_analysis suffix assumed\n";
}

my $host = `hostname`;
chomp($host);

if ($opt_u ) {
	$mus_increment = $opt_u;
	$opticalrotation_increment = 0;
} elsif ($opt_g) {
	$mus_increment = 0;
	$opticalrotation_increment = $opt_g;
} else {
	$mus_increment = 5;
	$opticalrotation_increment = 0;
	print "mu_s increment of 5 assumed, no optical activity\n";
}


if ($opt_o ) {
	$overwrite = 1;
}

if ($opt_n){
	$acceptanceCosine = $opt_n;
} else {
	$acceptanceCosine = 0
}

if ($opt_m){
	$multipleruns = $opt_m;
} else {
	$multipleruns = 1;
}

# Number of threads that will be started
if ($opt_t ) {
	$howmany = $opt_t;
} else {
	# Try to guess from machine we are on

	$howmany = int `cat "$ENV{HOME}/.mc-threads-$host"`;
	
	if ( $howmany == 0) {
		$howmany = 7;
	}
}
# Save to file, so can be modified by user
`echo "$howmany" > "$ENV{HOME}/.mc-threads-$host"`;

if (! ($setup || $run || $analyse) ) {
	print "You probably want to do at least one thing (setup, run or analyse) -s -r or -a\n";
	exit;
}

if (! -e $basedir) {
	die "Can't even find the base directory $basedir\nMake sure you run from top-level directory\n";
}

#
# Begining of the actual program: clean, setup, run or analyze
#

#
# Clean
#

if ($clean) {
	die "not done yet";
	for ($i = $firstsuffix; $i <= $lastsuffix; $i++) {
	
		if ($i == $basesuffix) {
			next;
		}
		
		$dir = GetDirectory($prefix, $i);
		print "Cleaning up $dir\n";	
	
		if ( ! -e $dir ) {
			print "Creating directory $dir\n";
			`mkdir $dir`;
		}
	}
}


#
# Setup
#

if ($setup ) {
	for ($i = $firstsuffix; $i <= $lastsuffix; $i++) {
	
		$dir = GetDirectory($prefix, $i);
		print "Setting up $dir\n";	

		if ($i != $basesuffix) {
		
			if ( ! -e $dir ) {
				print "Creating directory $dir\n";
				`mkdir $dir`;
			}
		
			if ( -e "$dir/param.dat") {
				if (! $overwrite ) {
					print " File param.dat, will NOT overwrite (use -o option)\n";
				} else {
					print " File param.dat exists! will overwrite\n";
					`mv -f $dir/param.dat $dir/param.dat.old`;
					`cp -f $basedir/param.dat $dir/`;
				}
			} else {
				`cp $basedir/param.dat $dir/`;
			}

			if ($mus_increment) {		
				$mu_s = $mus_increment*($i-$basesuffix);
				`perl -pi -e "s/mu_s=.*/mu_s=$mu_s/g" $dir/param.dat`;
			} elsif ($opticalrotation_increment) {
				$opticalrotation = $opticalrotation_increment*($i-$basesuffix);
				`perl -pi -e "s/rotationPerCmInClearSpace=.*/rotationPerCmInClearSpace=$opticalrotation/g" $dir/param.dat`;
			}
		}
		
		if ($multipleruns > 1){
			print " Set up: multiple runs (putting a file called multipleruns.dat in directory with directory names in it)\n";
			if ( -e "$dir/multipleruns.dat") {
				if (! $overwrite ) {
					print " File multipleruns.dat, will NOT overwrite (use -o option)\n";
				} else {
					print " File multipleruns.dat exists! will overwrite\n";
					`rm -f $dir/multipleruns.dat`;
				}
			}
			for ($j = 1; $j <= $multipleruns; $j++) {
				print "Making sub-directory $dir/run_$j\n";
				if ( -e "$dir/run_$j/") {
					if (! $overwrite ) {
						print " Directory $dir/run_$j exists, will NOT overwrite (use -o option)\n";
					} else {
						print " Removing previous $dir/run_$j directories\n";
						`rm -rf $dir/run_$j/`;
					}
				}
				`mkdir -p $dir/run_$j/`;
				`cp -f $dir/param.dat $dir/run_$j/`;
				`echo "run_$j" >> $dir/multipleruns.dat`; 
			}
		}
	}
}

#
# Run: 
# We build a list of jobs to run, then we cycle through them

if ($run) {
	my @jobsDir=();

	$current = $firstsuffix;
	while ($current <= $lastsuffix) {
		$dir = GetDirectory($prefix, $current);
		if (-e "$dir/multipleruns.dat") {
			unless (open MULTILIST, "$dir/multipleruns.dat") {
				die "Can't open file: $!";
			}
			while ($theSubDir = <MULTILIST>){
				chomp($theSubDir);
				if ($theSubDir =~ /^\#(.*)/) {
					print "Directory $1 commented out.  Skipping.\n";
				} else {
					push @jobsDir,"$dir/$theSubDir";
				}
			}
			close MULTILIST;
		} else {
			push @jobsDir,"$dir/";
		}
		$current++;
	}
	
	my @childs = ();
	my $theJobDir;
	
	while (scalar  @jobsDir) {
		print "Currently ".scalar @childs."/$howmany processes running: ".join " ", @childs,"\n";
		if (scalar @childs >= $howmany) { 
			# Could put something here to check load and adjust $howmany
			sleep $polltime;
			$howmany = `cat "$ENV{HOME}/.mc-threads-$host"`;
		} else {
			$theJobDir = shift @jobsDir;
		
		   if ($pid = fork) {
				# parent
				$current++;
				unshift @childs, $pid;
			} elsif (defined $pid ) {
				# child
				print `cd $theJobDir; run_all.sh`;
				print `perl -pi -e "s|^($theJobDir)$|\\#\\1|" multipleruns.dat`;
				exit;
			} else {
				die "Can't fork: $!\n";
			}
		}
	@childs = UpdateChildStatus($host, @childs);
	}

	while (scalar @childs) {
		@childs = UpdateChildStatus($host, @childs);
		print "Waiting for completion of ".scalar @childs." processes: ".join " ", @childs;
		print "\n";
		sleep $polltime;
	}
}

#
# Analyze
#

if ($analyse) {
	if ($analyse == 1) {
		for ($i = $firstsuffix; $i <= $lastsuffix; $i++) {
			$dir = GetDirectory($prefix, $i);
			print "Analysing $dir\n";	
		
			if (-e "$dir/multipleruns.dat") {
				unless (open MULTILIST, "$dir/multipleruns.dat") {
					die "Can't open file $dir/multipleruns.dat: $!";
				}
				while ($theSubDir = <MULTILIST>){
					chomp($theSubDir);
					$theSubDir =~ s/\#//;
					print `cd $dir/$theSubDir; extractallsubfiles.sh `;
				}
				close MULTILIST;
				
				print `cd $dir/; combineMultipleRuns.pl`;
			} else {
				print `cd $dir; extractallsubfiles.sh `;
			}
		}
	}
	
	if ($analyse == 2) {
		$analysis_dir = $prefix.$analyzesuffix;
	
		if ( ! -e $analysis_dir ) {
		    print "Creating directory $analysis_dir\n";
		    `mkdir $analysis_dir`;
		} else {
		    print "Directory $analysis_dir exists\n";
		}
		
		for ($acceptanceWidth = 0; $acceptanceWidth < 0.5; $acceptanceWidth += 0.05) {
		
			print "Process acceptance: $acceptanceWidth\n";
			for ($i = $firstsuffix; $i <= $lastsuffix; $i++) {
			
				$dir = GetDirectory($prefix, $i);
		
				print "Processing $dir\n";
				if (1) { #-e "$dir/Out*.dat"
					if ($i == $firstsuffix) {
						print `cd $dir/;extractstatsfromarrays.pl -w $acceptanceWidth -n $acceptanceCosine -h > ../$analysis_dir/stats_a$acceptanceWidth.txt`;
					} else {
						print `cd $dir/;extractstatsfromarrays.pl -w $acceptanceWidth  -n $acceptanceCosine >> ../$analysis_dir/stats_a$acceptanceWidth.txt`;
					}
				} else {
					print "Can't find files\n";
				}
			}
			
		}
	}
	
}

#system('warnme "DoAllThat @ARGV just finished"');

