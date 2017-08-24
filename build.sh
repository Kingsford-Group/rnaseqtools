#!/bin/bash
 
libs="gtfalign gtfformat gtfsimulator gtfcuff gtfmerge gtfquant"
dir=`pwd`

for i in `echo $libs`
do
	cur="$dir/$i"
	echo $cur
	cd $cur
	aclocal
	autoconf
	autoheader
	automake -a
	./configure
	make
	cd $dir
done

