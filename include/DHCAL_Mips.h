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

struct thr_t{
  int t_1;
  int t_2;
  int t_3;
};

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
  
  void      asic_map(int plate, double *p, TCluster clus);
  
  
  double    get_incident_angle(double *parameter);
  bool      set_geometry(std::string xml_geom_file);
  double    get_track_lenght();
  double    get_euclidean_dist(double* point1, double* point2);
  void      save_plots(std::string file="perfomance_plots.root");
  void      clear();
  void      set_threshold(int plate=-1,int t_1=170,int t_2=500, int t_3=345);
  void      set_xml_threshold(std::string filename);
  bool      read_xml_scan(std::string xmlfile);
  void      set_run_number(int run_id);
  int       get_nfired_layer();
  bool      do_characterization();
  bool      penetrated_muon_test();
  /**
     @param 3 DAQ threshold values 
     @return The values of threshold in pico-Coulomb
   */
  double    get_threshold_value(int threshold_id, int thr);
  double    set_id_decoder_string(std::string id_string = "I:7,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7"){_id_decoder_string = id_string;}
  
  static DHCAL_Mips* instance();
  
private:
  int _nplate;
  double _cut_region;
  int _nhit_on_clus_cut;
  int _angle_cut;
  int _run_id;
  int _nfired_layer;

  bool _nfired_layer_test;
  
  std::string _id_decoder_string;
  map<int, double > _plate_positions;
  map<int, double > _plate_threshold_1;
  map<int, double > _plate_threshold_2;
  map<int, double > _plate_threshold_3;
  
  map<int, map<int, thr_t > > scan_all;
  
  std::vector<TCluster> _cluster;
  HistBooker* hist_book;
  static DHCAL_Mips* _me;
};

#endif

