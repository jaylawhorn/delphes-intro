//--------------------------------
// DEMO READ IN TEXT FROM TEXT FILE
//--------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__) 
#include <TROOT.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#endif

void read_txt(TString infile="read_txt.txt" // inputfile
	      ) {

  // loop variable
  TString filename;

  // define a file stream
  ifstream ifs;
  // open a file stream
  ifs.open(infile);
  // make sure the file is actually open, exit the program if not
  assert(ifs.is_open());
  // another loop variable
  string line;
  while (getline(ifs, line)) { // loop over lines in the file stream
    // make a string stream to read each line
    stringstream ss(line);
    // read the first chunk of the stream to the variable "filename"
    ss >> filename;
    // print to screen "root://eoscms.cern.ch/filename" followed by a new line
    cout << TString("root://eoscms.cern.ch/")+filename << endl;
  } // end line loop
  // close the file stream
  ifs.close();

}
