#! /bin/sh

if [ $# = 0 ] ; then
	list="OutS1000.dat OutS1100.dat OutS1010.dat OutS1001.dat"
else
	list=$*
fi

#while getopts n: c
#	do case $c in
#		n)	acceptanceCosineIndex=$OPTARG;; 
#		\?)  echo $USAGE 
#		exit 2;; 
#esac done 

#if [ test  $acceptanceCosineIndex ]; then
#	echo "Warning: setting acceptance cosine to 0"
#	acceptanceCosineIndex="0"
#fi

for file in $list 
do
	if [ -f $file ]; then
		echo "Processing $file in `pwd`"
		java org.apache.xalan.xslt.Process -in $file -xsl $HOME/share/xml2flatfiles.xslt
	fi
done

#extractmuellermatrixelements.pl $list
