#!/usr/bin/ksh
if [ $# -lt 1 ]; then
	echo
	echo "USAGE:"
	echo "    doEndlZA.com T yofile_1 yofile_2 ... yofile_n"
	echo
	echo "PARAMETERS:"
	echo "    T         Temperature to heat to in MeV"
	echo "    yofile_i  an ENDL cross-section file (i.e. i=0)"
	echo 
	echo "doEndlZA.com heats the cross-sections in the yofile_i file to temperature"
	echo "T (MeV) using nuc_xsec_adjust_for_heated_target_endl."
	echo
	echo "EXAMPLE:"
	echo "doEndlZA.com 1.2e-2 yo00c10i000s000"
	echo "    This command heats the cross-sections in yo00c10i000s000 to 12keV (1.2e-2MeV)"
	echo "    and write the results to file yo00c10i000s000_T12.000keV."
	echo 
	exit 0
fi
Tmp=$1
shift
if [ $# -le 0 ]; then
	l=`ls *i000*`
	for f in $l; do nuc_xsec_adjust_for_heated_target_endl $Tmp $f; done
else
	while [ $# -gt 0 ]; do
		if [ -d $1 ]; then
			l=`ls $1/*i000*`
		else
			l=$1
		fi
		for f in $l; do nuc_xsec_adjust_for_heated_target_endl $Tmp $f; done
		shift
	done
fi
