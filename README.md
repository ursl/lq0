Code example for the LQ analysis. To start from scratch (assuming you
have github access):

1. create a work directory: 
```
  cd
  mkdir lq
  cd lq
```

2. setup a CMS release
```
  scramv1 project CMSSW_5_3_20
  cd CMSSW_5_3_20
  cmsenv
  cd - 
```

3. install Delphes 
```
  wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.1.2.tar.gz
  tar zxvf Delphes-3.1.2.tar.gz
  cd Delphes-3.1.2
  perl -pi -e 's%#include "Pythia.h"%#include "Pythia8/Pythia.h"%g' modules/PileUpMergerPythia8.cc readers/DelphesPythia8.cpp
  make -j4
```

4. get user code
```
  git clone git@github.com:ursl/util  
  git clone git@github.com:ursl/lq0
  cd util && make && cd - 
  cd lq0 && make
```
