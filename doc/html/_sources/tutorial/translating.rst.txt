Basic translating between ENDF and GNDS
=======================================

Fudge comes with a tools for translating ENDF-6 evaluations to and from GNDS which reside in the brownies/bin directory.

Translating from ENDF into GNDS
-------------------------------

The Python script brownies/bin/endf2gnds.py translations an ENDF-6 file to a GNDS file. This script writes a GNDS file containing 
the reactionSuite node and, if covariance data are present in the ENDF-6 file, its also writes a GNDS file containing the 
covarianceSuite node. For example, the following command translating the ENDF/B-VIII.0 file n-001_H_001.endf and outputs the 
files n-001_H_001.xml and n-001_H_001-covar.xml::

    /full/path/to/FUDGE/brownies/bin/endf2gnds.py n-001_H_001.endf n-001_H_001.xml

If using the *pip installed* version of **FUDGE**, the one only needs the::

    brownies/bin/endf2gnds.py n-001_H_001.endf n-001_H_001.xml

Translating GNDS files back into ENDF
-------------------------------------

The Python script brownies/bin/gnds2endf.py translations a GNDS file to an ENDF-6 file. This script will look for a corresponding 
covariance file and, if present, add its data to the outputted ENDF-6 file. For example, the following translates the n-001_H_001.xml 
file generated in the endf2gnds.py comannd above back into an ENDF-6 file::

    /full/path/to/FUDGE/brownies/bin/gnds2endf.py n-001_H_001.xml

Reading GNDS XML files
----------------------

If I didn't have pre-made instances of ``reactionSuite`` and ``covarianceSuite``, how would I read in the XML files?
For this purpose, both the ``fudge.reactionSuite`` and ``fudge.covariances`` have the factory function ``readXML()``.
To use them do:

    >>> from fudge import reactionSuite
    >>> from fudge.covariances import covarianceSuite
    >>> n_H1 = reactionSuite.readXML( "n-001_H_001.gnds.xml" )

This reads in the evaluation itself.  To read in the covariances, we need to tell the `covariances.readXML()` function
where the evaluation is so that it can set up the hyperlinks correctly:

    >>> myOtherCov = covarianceSuite.readXML( "n-001_H_001.gndsCov.xml", reactionSuite=myOtherEval )
