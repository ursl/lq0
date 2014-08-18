#ifndef REDTREEDATA_H
#define REDTREEDATA_H

struct redTreeData {
  int type, channel; 
  double w8; 

  double gpm;    
  double gpm2;    
  double gppt;   
       
  double gnm;    
  double gnm2;    
  double gnpt;   
        
  double glqpm;  
  double gljpm;  

  double glqnm;  
  double gljnm;  

  // -- reco 
  double ljpm;  
  double ljppt;  

  double ljnm;  
  double ljnpt;  

  double l0pt; 
  double l1pt; 

  double j0pt; 
  double j1pt; 

  double st; 
  double mll;
  double mljetmin;

};

#endif
