# <<BEGIN-copyright>>
# <<END-copyright>>

.PHONY: clean realclean build extensions

EXTENSIONS = numericalFunctions crossSectionAdjustForHeatedTarget

default: inplace

build:
	python setup.py --quiet build

inplace:
	python setup.py --quiet build_ext --inplace build

clean:
	rm -rf build
	find . -name "*.pyc" -exec rm -f {} \;
	find . -name "*.so" -exec rm -f {} \;
	$(MAKE) extensions MODE=realclean

realclean: clean

extensions:
	SAVED_PWD=`pwd`; \
	for directory in $(EXTENSIONS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done

docs:
	epydoc-2.6 -o doc/code-ref --html --show-imports --exclude="fudge.processing.resonances.setup" --exclude="fudge.vis.gnuplot.fudge2dMultiPlot" --exclude="fudge.vis.gnuplot.endl[2-4]dplot" --no-frames --docformat="epytext" -v fudge
# --check --debug --graph="umlclasstree"

dist:
	python setup.py sdist --formats=gztar,zip
