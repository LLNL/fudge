# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: clean realclean

default: libcrossSectionAdjustForHeatedTarget.a

all: default

libcrossSectionAdjustForHeatedTarget.a: crossSectionAdjustForHeatedTarget.o
	$(AR) -cr libcrossSectionAdjustForHeatedTarget.a crossSectionAdjustForHeatedTarget.o

crossSectionAdjustForHeatedTarget.o: crossSectionAdjustForHeatedTarget.c
	$(CC) $(CFLAGS) -I. -c -o crossSectionAdjustForHeatedTarget.o crossSectionAdjustForHeatedTarget.c

clean:
	rm -f *.o *.a

realclean: clean
