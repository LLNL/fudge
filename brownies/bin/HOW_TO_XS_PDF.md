How to use xs_pdf.py
====================

Usage message:

    usage: xs_pdf.py [-h] [--source {data,eval,file}] [--endf ENDF] [--pdf PDF]
                     [--pdfFormat {PURR,default}] [--sym SYM] [--A A]
                     [--XSmin XSMIN] [--XSmax XSMAX] [--NXS NXS] [--useXSLethargy]
                     [--Emin EMIN] [--Emax EMAX] [--NE NE] [--useELethargy]
                     [--ignoreData] [--ignoreUncertainty] [--printSetStats]
                     [--saveSetStats SAVESETSTATS]
                     [--readSetModifications READSETMODIFICATIONS] [--editSets]
                     [--skipBadData] [--continuumSpectraFix] [--thicken THICKEN]
                     [--noRenormalize] [--check] [--getAveXS]
                     [--aveXSMode {quad,linear,const,pdf,endf,mughabghab}]
                     [--getXSPDF] [--getXSCorr] [--aveOutFile AVEOUTFILE]
                     [--pdfOutFile PDFOUTFILE] [--printXSPDFSlice PRINTXSPDFSLICE]
                     [--plotXSCorr] [--printXSCorr] [--plotXSPDF]
                     [--plotXSPDFSlice PLOTXSPDFSLICE] [--plotAveXS] [-v]
                     MT

Required positional arguments:

    MT                    The MT of the reaction to examine.

Misc. optional arguments:

    -h, --help            show this help message and exit
    -v                    Enable verbose output

Purpose: Get the cross section PDF $P(\sigma|E)$ given experimental data or evaluation 
and use it to do stuff.  The most obvious thing is to compute the average cross section:

        $\left<\sigma\right> = \int_0^\infty d\sigma P(\sigma|E) \sigma$

but this code can do many other (mostly random) things.

Let's go into the detailed options...

1. To read in data, first set the type of source file for the PDF using the "--source" 
   option:

        --source {data,eval,file}
                        Set where PDF comes from: data, evaluation or pre-
                        computed file
                        
   You will also have to name the actual source file.  This gets weird because the source 
   can be an ENDF file, the pre-computed URR from an ACE file processed with NJOY2016 or
   directly computed from experimental data in EXFOR.  Lets examine the three sources.

   * Using ENDF as a source is easy, use "--endf":

        --endf ENDF           The ENDF file to use (if given)

   * Using the pre-computed PDF's from ACE, use 

        --pdf PDF             The file with the URR pdf (if given)
        --pdfFormat {PURR,default}
                             The format of file with the URR pdf (if given)
                   
   * Finally, to just grab experimental data from EXFOR, we'll use the *x4i* code.  This is 
     all handled inside *xs_pdf.py*, but you'll need to tell *xs_pdf.py* which isotope you 
     want to investigate.  Do that with:                     

        --sym SYM             The symbol of the nucleus in question
        --A A                 The A of the nucleus in question

   To recap, you need to set the --source and then load the data using one of these three
   options.

   Note, if you gave and ENDF file already, *xs_pdf.py* can figure our the SYM and A from
   that.  So it's possible to give "--source data", the provide "--endf ENDF" to tell
   *xs_pdf.py* to load both an evaluation and data (useful if you want both on the same 
   plot).
  
2. You will also need to set the energy range to process.  To set the energy range, use:

        --Emin EMIN           Emin in MeV (default=1e-11 MeV)
        --Emax EMAX           Emax in MeV (default=20.0 MeV)

3. The averaging will take place in NE separate bins spaced between EMIN and EMAX.  You 
   can control the number of bins and whether the bins should be equally spaced or 
   logarithmically spaced with 

        --NE NE               Number of equal bins on the domain (Emin, Emax) (default=100)
        --useELethargy        Make the energy bins equal lethargy bins (default: False)

4. Computing the PDF itself.  First, tell *xs_pdf.py* that you want it by using 

        --getXSPDF            Compute the cross section probability distribution
   
   The PDF is computed by chopping the cross section into horizontal "bands" and then 
   averaging the cross section in that band.  To fine-tune the band placement, use

      --XSmin XSMIN         XSmin in b (default is to take from data)
      --XSmax XSMAX         XSmax in b (default is to take from data)
      --NXS NXS             Number of cross section bins (default=20)
      --useXSLethargy       Make the cross section bins equal lethargy bins (default: False)
      --noRenormalize       Do not renormalize the PDFs read in from a file

   For the averaging to be well behaved, make sure the average cross section goes right 
   through the central band and that the bands extend high and low enough (XSMIN & MXMAX) 
   to capture all the cross section fluctuations.  Also, make sure you have enough bands 
   to cover the variation of the PDF in between XSMIN and XSMAX.

5. Getting the average cross section.  This is computed using the cross section PDF which
   must be also computed, so use the --getXSPDF option too.

         --getAveXS            Compute the average cross section
         --aveXSMode {quad,linear,const,pdf,endf,mughabghab} Mode to compute the average

6. Getting the energy-energy cross section correlation R(e).  Don't worry about this option.

          --getXSCorr           Compute the cross section energy-energy correlation

7. Treatment of data

    * EXFOR data

              --ignoreData          Ignore the experimental data unless it's what's 
                                    getting averaged
              --ignoreUncertainty   Ignore the uncertainty on the experimental data
              --saveSetStats SAVESETSTATS    
                                    Store set information to this JSON file
              --readSetModifications READSETMODIFICATIONS    
                                    JSON file with potential file modifications
              --editSets            Edit the sets so that they are between Emin & Emax, etc
              
    * ENDF data

              --skipBadData         Skip bad data in the ENDF file if needed
              --continuumSpectraFix Apply the continuum spectrum fix to ENDF data if needed
              --thicken THICKEN     Thicken (reconstructed) cross sections by this number of 
                                    points between existing points (default = 0, i.e.
                                    don't thicken

8. Setting the output filenames

        --aveOutFile AVEOUTFILE The name of the output file for the average cross section
        --pdfOutFile PDFOUTFILE The name of the output file for the cross section PDF

9. What to output

          --check               Print out the PDf normalization inegral as a cross check 
                                and do other checks
          --printSetStats       Print EXFOR data set statistics
          --printXSPDFSlice PRINTXSPDFSLICE
                                Print a slice of the cross section PDF derived from the 
                                data in the bin specified
          --plotXSCorr          Plot the cross section energy-energy correlation
          --printXSCorr         Print the cross section energy-energy correlation
          --plotXSPDF           Plot cross section PDF derived from the data
          --plotXSPDFSlice PLOTXSPDFSLICE
                                Plot a slice of the cross section PDF derived from the 
                                data in the bin specified
          --plotAveXS           Plot experimental data and the average cross section
