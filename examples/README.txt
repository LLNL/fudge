The 'examples' directory contains:

n-009_F_019.xml		# sample GNDS files translated from ENDF-VII.1
n-026_Fe_056.xml
n-094_Pu_239.xml

n-094_Pu_240-covar.xml	# sample of new covariance format

plotCrossSection.py	# demonstrates resonance reconstruction and new plotting capability
compareCrossSections.py # compare cross sections for a single MT from two or more sources (which may be in ENDF or GNDS format)
testMultiplicities.py	# simple checking of GNDS-formatted data
newEvaluation.py        # script to demonstrate how to create an evaluation from scratch using the evaluator toolkit

gnds.xsl			# xml 'stylesheet', for viewing GNDS files in a web browser.
	for example, try:
	>xsltproc gnds.xsl n-094_Pu_239.xml > Pu239.html
	then view Pu239.html in a web browser

