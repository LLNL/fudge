# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check checke clean realclean

TARGETS = setting.py

default:

check:
	for f in $(TARGETS); do python $$f; done

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:

realclean: clean
