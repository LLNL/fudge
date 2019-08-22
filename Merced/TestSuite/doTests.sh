#!/bin/bash
#
# Run Merced on test input files, and compare results to baseline output files.
# If no arguments are supplied, runs all available input files.
# If one or more directories are supplied, only run tests in those directories.

if [ $# -gt 0 ]; then
	dirs=$@
else
	dirs=*
fi

for dir in $dirs; do
	if ! [[ -d $dir ]]; then
		continue
	fi

	echo Entering directory $dir:
	cd $dir;
	for fil in `ls in.*`; do
		merced $fil &> ${fil/in./}.info;
		if ! cmp ${fil/in/out} utfil >/dev/null 2>&1; then
			echo '  ' $fil output differs from baseline;
		fi
	done
	rm utfil;
	cd ..;
done
