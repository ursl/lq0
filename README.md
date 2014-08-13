Code example for the LQ analysis. To start from scratch (assuming you
have github access and a tcsh):

* create a work directory: 
```
  mkdir lq
  cd lq
```

* setup a CMS release (I think you can ignore the warning about the non-production architecture)
```
  scramv1 project CMSSW_5_3_20
  cd CMSSW_5_3_20
  cmsenv
  cd - 
```

* install Delphes (I think you can ignore the many compilation warnings)
```
  wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.1.2.tar.gz
  tar zxvf Delphes-3.1.2.tar.gz
  cd Delphes-3.1.2
  perl -pi -e 's%#include "Pythia.h"%#include "Pythia8/Pythia.h"%g' modules/PileUpMergerPythia8.cc readers/DelphesPythia8.cpp
  make -j4
  setenv DELPHES `pwd`
  cd - 
```
  

* get user code
```
  git clone git@github.com:ursl/util  
  git clone git@github.com:ursl/lq0
  cd util && make && cd - 
  cd lq0 && make
```

* run user code
```
  bin/runLq -f /STORE/LQ/pair/Events/run_01/tag_1_delphes_events.root
```
