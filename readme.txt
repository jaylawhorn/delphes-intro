DELPHES QUICKSTART - JAY LAWHORN 10/2/13

In rootlogon.C, change file paths to your own installation of delphes.

There are three different set of scripts in this folder - 

skeleton.C is a barebones introduction to using Delphes and ROOT.

read_txt.C is an introduction to reading configuration info stored in a text file.

selectDelphes.C (and run.sh/submit.sh) is an introduction to submitting jobs to lxbtch.

Run with:

    root (-l -q) skeleton.C+

You always want to use the + at the end - this forces root to compile the 
macro instead of using the sometimes sketchy interpreter. "-l" will make
root omit the splash screen on start up, and "-q" will make root quit when
it's done executing the file. The parenthesis denote that those arguments are
optional.

Delphes Resources:
Delphes/examples/Example1.C , Example2.C , Example3.C
https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook (particularly 2. ROOT tree description)
Delphes/classes/*

PDG codes: 
http://pdg.lbl.gov/2002/montecarlorpp.pdf
