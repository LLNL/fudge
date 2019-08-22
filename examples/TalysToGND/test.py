#!/usr/local/bin/python

## Run TALYS in default mode and export to ENDL:
import geft

yi = 1 # 1 = neutron
ZA = 36077   # Kr77
#ZA = 94239   # Pu239

t = geft.talysZAClass.talysZA(ZA,yi)            # a proxy for an endlZA with some extras for running talys
#t.run() # runs talys

### New methods to access GND format of the data
### calls are a little awkward since the reaction suite is contained in the t instance,
### unlike the endl route where the endl files are listed within the a instance.
a = geft.talysGndData.CompleteGndEvaluation(t)
t.write("test")


### The old way to write endl output
### hoping to preserve backward compatibility with this old approach
#a = geft.talysGndData.CompleteEvaluation(t)     # a form of endlList, a base class list of data elements
#a.toEndl()                                      # function that controls the writing of data to endlFiles
#t.process()                                     # runs some processing and cleanups
#t.save()                                        # inherited from the endlZA
