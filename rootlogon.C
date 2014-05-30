{
  // replace with location of your libDelphes.so
  gSystem->Load("/afs/cern.ch/user/j/jlawhorn/Delphes/libDelphes.so");
  // replace with location of your Delphes folder
  gROOT->ProcessLine(".include /afs/cern.ch/user/j/jlawhorn/Delphes");
  // replace with location of your Delphes/external folder
  gROOT->ProcessLine(".include /afs/cern.ch/user/j/jlawhorn/Delphes/external");
}
