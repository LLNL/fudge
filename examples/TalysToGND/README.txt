Generate GND file from TALYS output.

To execute, first copy the talys structure directory into 'talys' (or edit code to point to your TALYS installation)
Then, execute 'python test.py'
This will run TALYS, read the output from talys/yi01/za***, and convert into GND.

Some notes about the code:

talysGndData.CompleteGndEvaluation does most of the work:
  loops through talysZA.reactionList, adds all reactions to various lists
  (Summed, Elastic, Fission, Binary and Plain)

Features needing more work: capture, fission, adding discrete gammas to the
particle list
