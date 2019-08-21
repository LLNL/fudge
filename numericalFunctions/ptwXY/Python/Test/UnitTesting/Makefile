# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check clean realclean doTarget

subDirectories = thicken thinning integrate convolution Others slicing settingData binaryMath interpolation pythonOnly toFunction

default: 

check:
	$(MAKE) doTarget MODE=check

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
