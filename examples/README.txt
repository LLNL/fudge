The 'examples' directory contains:

plotCrossSection.py		# demonstrates resonance reconstruction and new plotting capability
compareCrossSections.py # compare cross sections for a single MT from two or more sources (which may be in ENDF or GNDS format)
multiTemperature.py     # plot one cross section at various temperatures to see impact of Doppler broadening

gnds.xsl			# xml 'stylesheet', for viewing GNDS files in a web browser.
	for example, try:
	>xsltproc gnds.xsl n-094_Pu_239.xml > Pu239.html
	then view Pu239.html in a web browser

