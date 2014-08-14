#ifndef REDTREEDATA_H
#define REDTREEDATA_H

struct redTreeData {
  int type; 
  double w8; 

  double gpm;    
  double gpm2;    
  double gppt;   
  double gpeta;  
       
  double gnm;    
  double gnm2;    
  double gnpt;   
  double gneta;  
        
  double glqpm;  
  double gljpm;  

  double glqnm;  
  double gljnm;  

  // -- reco 
  double ljpm;  
  double ljppt;  
  double ljpeta;  

  double ljnm;  
  double ljnpt;  
  double ljneta;  

  double st; 
  double mll;
  double mljetmin;

};

#endif
