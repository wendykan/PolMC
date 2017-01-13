#!/bin/sh

# Runs polmc with param.dat from the same directory,
# keeps a log of everything and emails the job owner when done.

# The notation of the commands might look funny
# to non Jedi knights:  2>&1 | mylogger is a redirect of
# stderr to the program mylogger (which is a personal program
# that keeps a log file of everything).

# Deal with the different version of echo in the Linux cluster (lc1)
if [ -f "/bin/echo" ] ; then
    echo="/bin/echo -e"
else
    echo="/usr/bin/echo"
fi


if [ -f "./polmc" ] ; then
	$echo "polmc program in current directory.  Check what version you want to use."
	exit
fi

if [ -f "param.dat" ] ; then
	# Run the actual program
	polmc param.dat 2>&1 | mylogger

	# When done, send an email to the user who owns this job,
	# hoping this is actually an email address
	$echo "Subject: Results from run_all.sh\nJob started from:\n`pwd`\ncompleted at:\n`date`" | mail `whoami` > /dev/null 
else
	$echo "Error: I expect the parameter file param.dat to be in the current directory: `pwd`"
	exit
fi





