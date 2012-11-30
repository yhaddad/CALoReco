#ifndef _DHCAL_MIPS_H
#define _DHCAL_MIPS_H

#include <map>
#include <TMath.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <EVENT/CalorimeterHit.h>
#include "UTIL/CellIDDecoder.h"
#include "TRACK/MuonTrack.h"
#include "TCluster.h"
#include "HistBooker.h"


typedef EVENT::CalorimeterHit CaloHit;
typedef std::vector<CaloHit*> CaloHitVec;

using namespace std;

/** 
    \class DHCAL_Mips
    \author  Y.Haddad
    \date Nov 2012
    \version 0.1
    
    \brief DHCAL Mips studies 
    
    <h2> Description</h2>
    HDCAL Mips reconstruction and callibration 
    
*/

class DHCAL_Mips
{
 public:
  DHCAL_Mips();
  //DHCAL_Mips(std::vector<TCluster> cluster, int nlayer, double cut_region);
  /**
     @param 
     @return
  */
  void      set_nlayer(int nlayer);
  void      set_nhit_on_clus_cut(int cut);
  void      set_angle_cut(double cut);
  void      set_region_cut(double cut);
  void      set_clusters(std::vector<TCluster> clus);
  bool      do_calibration();
  void      get_charge_shape(CaloHitVec hits, double *p);
  
  void      asic_eff_map();
  void      asic_mul_map();
  
  double    get_incident_angle(double *parameter);
  void      set_geometry(std::map<int, double> geometry);
  double    get_track_lenght();
  double    get_euclidean_dist(double* point1, double* point2);
  void      save_plots();
  void      clear();

  static DHCAL_Mips* instance();
 private:
  int _nplate;
  double _cut_region;
  int _nhit_on_clus_cut;
  int _angle_cut;
  map<int, double > _plate_positions;
  
  std::vector<TCluster> _cluster;
  HistBooker* hist_book;
  static DHCAL_Mips* _me;
};

#endif
