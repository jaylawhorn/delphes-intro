{
  gROOT->ProcessLine(".include /afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes");  
  gROOT->ProcessLine(".include /afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes/external");
  // replace with location of your libDelphes.so
  //gSystem->Load("libDelphes");
  gSystem->Load("/afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes/libDelphes.so");
  // replace with location of your Delphes folder
  
  // replace with location of your Delphes/external folder
  //
}
