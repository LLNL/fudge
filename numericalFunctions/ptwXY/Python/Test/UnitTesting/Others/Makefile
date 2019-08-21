# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check checke clean realclean

TARGETS = clip.py trim.py xArray.py dullEdges.py mutualify.py mod.py areDomainsMutual.py reallocate.py mergeClosePoints.py copy.py

default:

check:
	for f in $(TARGETS); do python $$f; done

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:

realclean: clean
