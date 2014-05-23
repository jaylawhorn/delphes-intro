DELPHES QUICKSTART - JAY LAWHORN 10/2/13

In rootlogon.C, change file paths to your own installation of delphes.

Before running, need to setup cms environment. This can be done by 

    cd (your CMSSW_x_x_x/src directory)

    cmsenv

Run with:

    root (-l -q) skeleton.C+

You always want to use the + at the end - this forces root to compile the 
macro instead of using the sometimes sketchy interpreter. "-l" will make
root omit the splash screen on start up, and "-q" will make root quit when
it's done executing the file. The parenthesis denote that those arguments are
optional.

Delphes Resources:
Delphes-3.0.10/examples/Example1.C , Example2.C , Example3.C
https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook (particularly 2. ROOT tree description)
Delphes-3.0.10/classes/*

PDG codes: 
http://pdg.lbl.gov/2002/montecarlorpp.pdf

HHToBBTauTau Study:
http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/HHBBTauTauStudy