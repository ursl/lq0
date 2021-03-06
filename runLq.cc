#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "anaLq.hh"

using namespace std;


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: please see https://github.com/ursl/lq0
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "==> runLq: Running under process ID  " << processID << endl;

  string progName  = argv[0]; 
  string writeName, fileName, jsonName;
  int file(0), json(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  int randomSeed(processID);
  int verbose(-99); 
  int blind(1); 
  int isMC(0);
  int year(0); 

  // Change the MaxTreeSize to 100 GB (default since root v5.26)
  TTree::SetMaxTreeSize(100000000000ll); // 100 GB

  // -- Some defaults
  string dirBase("./");               // this could point to "/home/ursl/data/root/."
  string dirName("."); dirspec = 0;   // and this to, e.g. "bmm", "bee", "bem", ...
  string cutString("nada");

  string treeName("Delphes");

  string readerName("bmmReader"); 
  TString histfile("");

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i],"-h")) {
        cout << "List of arguments:" << endl;
        cout << "-b {0,1}      run blind?" << endl; 
        cout << "-c filename   chain definition file" << endl;
        cout << "-C cutstring  string with cuts" << endl;
        cout << "-D path       where to put the output" << endl;
        cout << "-f filename   single file instead of chain" << endl;
        cout << "-m            use MC information" << endl;
        cout << "-n integer    number of events to run on" << endl;
        cout << "-r class name which tree reader class to run" << endl;
        cout << "-s number     seed for random number generator" << endl;
        cout << "-S start     starting event number" << endl;
        cout << "-o filename   set output file" << endl;
        cout << "-v level      set verbosity level" << endl;
        cout << "-y year       set year" << endl;
        cout << "-h            prints this message and exits" << endl;
        return 0;
    }
    if (!strcmp(argv[i],"-b"))  {blind      = atoi(argv[++i]); }                 // run blind?
    if (!strcmp(argv[i],"-c"))  {fileName   = string(argv[++i]); file = 0; }     // file with chain definition
    if (!strcmp(argv[i],"-C"))  {cutString  = string(argv[++i]);           }     // file with cuts
    if (!strcmp(argv[i],"-D"))  {dirName    = string(argv[++i]);  dirspec = 1; } // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-j"))  {jsonName   = string(argv[++i]); json = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-m"))  {isMC       = 1; }                               // use MC information?
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run 
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-S"))  {start = atoi(argv[++i]); }                      // set start event number
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
    if (!strcmp(argv[i],"-v"))  {verbose    = atoi(argv[++i]); }                 // set verbosity level
    if (!strcmp(argv[i],"-y"))  {year       = atoi(argv[++i]); }                 // set year
  }


  // -- Determine filename for output histograms and 'final' small/reduced tree
  TString meta = fileName;
  if (histfile == "") {
    TString  barefile(fileName), chainFile, meta;
    if (file == 0) {
      // -- input from chain
      if (barefile.Contains("chains/")) {
        barefile.ReplaceAll("chains/", "");
        histfile = barefile + "." + ".root";
        if (dirspec) {
          if (dirName[0] == '/') {
            histfile = dirName + "/" + histfile;
          } else {
            histfile = dirBase + "/" + dirName + "/" + histfile;
          }
        }
      } else {
        histfile =  barefile + "." + ".root";
        if (dirspec) {
          if (dirName[0] == '/') {
            histfile = dirName + "/" + histfile;
          } else {
            histfile = dirBase + "/" + dirName + "/" + histfile;
          }
        }
      }
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla + "." + ".root";
      if (dirspec) {
        histfile = dirBase + "/" + dirName + "/" + histfile;
      }
    }  else if (file == 1) {
      // -- single file input
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla;
      histfile.ReplaceAll(".root", "");
      histfile +=  ".root";
      if (dirspec) {
        if (dirName[0] == '/') {
          histfile = dirName + "/" + histfile;
        } else {
          histfile = dirBase + "/" + dirName + "/" + histfile;
        }
      }
    }
  }

  cout << "==> runLq: Opening " << histfile.Data() << " for output histograms" << endl;
  cout << "==> runLq: Opening " << fileName.c_str() << " for input" << endl;


  // -- Set up chain
  TChain *chain = new TChain(TString(treeName));
  cout << "==> runLq: Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    // -- non-trivial chain input
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
      if (meta.Data()[0] == '#') continue; 
      sscanf(meta.Data(), "%s %d", pName, &nentries); 
      if (nentries > -1) {
        cout << "       " << pName << " -> " << nentries << " entries" << endl; 
        chain->Add(pName, nentries); 
      } else {
        cout << meta << endl;
        chain->Add(meta); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    // -- single file input
    cout << "           " << fileName << endl;
    chain->Add(TString(fileName));
  }

  anaLq *a = new anaLq(chain); 
  
  if (a) {
    string shistfile = histfile.Data(); 
    a->openHistFile(shistfile); 
    if (cutString.compare("nada")) {
      a->setCuts(cutString); 
    }

    a->startAnalysis(); 
    a->loop(nevents, start);
    a->endAnalysis();
    a->closeHistFile(); 
  } else
    cerr << "Readerclass '" << readerName << "' not found" << endl;

  delete a; // so we can dump some information in the destructor
  
  return 0;
}

