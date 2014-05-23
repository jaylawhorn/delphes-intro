{
  // replace with location of your libDelphes.so
  gSystem->Load("/afs/cern.ch/user/k/klawhorn/Delphes-3.0.10/libDelphes.so");
  // replace with location of your Delphes-3.x.x folder
  gROOT->ProcessLine(".include /afs/cern.ch/user/k/klawhorn/Delphes-3.0.10");
  // replace with location of your Delphes-3.x.x/external folder
  gROOT->ProcessLine(".include /afs/cern.ch/user/k/klawhorn/Delphes-3.0.10/external");
}
