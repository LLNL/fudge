Fudge Interactive Plotting
==========================

The FUDGE Multiplot Module Inheretence Diagram
----------------------------------------------

.. graphviz::
  
   digraph {
       rankdir = BT

       MultiPlotWithPyQt5 [label="MultiPlotWithPyQt5(plotTitle, xLabel, yLabel, plotData)", shape=box]


       pyQtMultiplot [label="pyQtMultiplot(plotTitle, xLabel, yLabel, plotData)", shape=box]
       pyQtMultiplot -> MultiPlotWithPyQt5
   
       plotCanvas [label="plotCanvas(parent, width, height, xLabel, yLabel, plotData, dpi)", shape=box]
       plotCanvas -> pyQtMultiplot
   
       DialogWindow [label="DialogWindow(width, height, plotObject, parent)", shape=box]
       DialogWindow -> pyQtMultiplot

       PlotOptions [label="PlotOptions(plotObject)", shape=box]
       PlotOptions -> DialogWindow   

       DataOptions [label="DataOptions(plotObject)", shape=box]
       DataOptions -> DialogWindow    

   }


Multiplot Module: Python Source Documentation 
---------------------------------------------

.. automodule:: fudge.vis.pyqt5.multiplot
   :members:
   :undoc-members:
   :show-inheritance: