#ifndef _HISTBOOKER_H
#define _HISTBOOKER_H

#include <map>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <string>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include <sstream>
/** 
\class HistBooker
  \author  Y.Haddad
  \date Nov 2012
  \version 0.1
  
  \brief ROOT histogram handling 
  
  <h2> Description</h2>
  Histogramer booker, and store it on the root file.
  
*/

class HistBooker
{
 public:
  HistBooker();
  /**
     @param TH1F,2F,3F parameters
     @return a TH1*,TH2*,TH3* ptr
     
  */
  TH1* book_h1(std::string name,int nbin,
	       double xmin,double xmax);
  TH2* book_h2(std::string name,int nbinx,int nbiny,
	       double xmin,double xmax,
	       double ymin,double ymax);
  TH3* book_h3(std::string name,int nbinx,int nbiny,int nbinz,
	       double xmin,double xmax,
	       double ymin,double ymax,
	       double zmin,double zmax);
  /**
     @param Tprofile1D,2D,3D parameters
     @return a TH1*,TH2*,TH3* ptr
  */
  
  TH1* book_profile1(std::string name,int nbinx,double xmin,double xmax,double ymin,double ymax);
  TH2* book_profile2(std::string name,int nbinx,double xmin,double xmax,int nbiny, double ymin,double ymax, double zmin, double zmax);
  
  /** 
      @param name Name of the histogram
      @return a TH1* ptr
  */
  
  TH1* get_h1(std::string name);
  TH2* get_h2(std::string name);
  TH3* get_h3(std::string name);
  
  /**
     @param Save the maps of Histograms on ROOT file
     @return hitograms.root
  */
  
  void write_hitograms(std::string file);
  std::vector<std::string> split(std::string string_list, char sep );
  
  static HistBooker* instance();
 private:
  std::map<std::string,TH1*> _map_h1;
  std::map<std::string,TH2*> _map_h2;
  std::map<std::string,TH3*> _map_h3;
  static HistBooker* _me;
  
};

#endif
