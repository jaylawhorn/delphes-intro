{
  // replace "/afs/cern.ch/work/j/jlawhorn/public/forXLL" with your directory containing delphes
  gROOT->ProcessLine(".include /afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes");  
  gROOT->ProcessLine(".include /afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes/external");

  gSystem->Load("/afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes/libDelphes.so");

}
