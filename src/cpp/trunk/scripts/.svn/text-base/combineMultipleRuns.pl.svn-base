#!/usr/local/bin/perl 

use MyLibs;

#List directories of multiple runs
my @jobsDir=();
	
if (-e "multipleruns.dat") {
	unless (open MULTILIST, "multipleruns.dat") {
		die "Can't open file multipleruns.dat: $!";
	}
	while ($theSubDir = <MULTILIST>){
		chomp($theSubDir);
		push @jobsDir,"$theSubDir";
	}
	close MULTILIST;
} else {
	die "No multipleruns.dat file\n";
}

print join "\n",@jobsDir;

# Pick files to sum.
opendir(DIR, $jobsDir[0]) || die "can't opendir $: $!";
    @all_files = readdir(DIR);
closedir DIR;

foreach $file (@all_files) {
	if ($file =~ /[IQUV]\.dat$/){
		my @filesToAdd = ();
		
		foreach $dir (@jobsDir){
			push @filesToAdd, "$dir/$file";
		}
		
		AddFiles($file, @filesToAdd);
	}
}
