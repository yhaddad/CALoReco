#ifndef _HitProcess_hh_
#define _HitProcess_hh_


#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "IO/LCWriter.h"
#include "TCluster.h"
#include <map>
#include "DHCAL_Mips.h"

struct calibConst{
  int Ntrack;
  double m_Nhit;
  double m_Multiplicity;
  double m_Efficiency;
};


class HitProcess  : public marlin::Processor
{
public:
  Processor*  newProcessor() { return new HitProcess ; }
  HitProcess();
  ~HitProcess() {};
  void init();
  void processEvent( LCEvent * evtP );
  //void processRunHeader( LCRunHeader * runH);// added by me 
  bool fillCallibration(int K,double RegionCut,double *position,double *p);
  void effMap(EVENT::CalorimeterHitVec hits,double* p);
  double distance(double z,double x,double y,double *p);
  
  bool penetratedMuons(std::vector<TCluster> cluster, int minlayer, int minclus);
  void padEffMap(int Layer,double *position, double* p);
  void asicEffMap(int Layer,double multiplicity,double *position, double* p);
  int  getNhitOnClusterSet(std::vector<TCluster> cluster);
  //int* getNhitOnLayer();
  double getTrackLenght(std::vector<TCluster> cluster);
  // --------------
  // cleaning methods
  // ---------------
  double distEucl(double* point1,double* point2);
  /* remove the farest hits using hasdorf distance 
   * dist (pi, set of point ) = min(pi, {pj}_set )
   */
  std::vector<TCluster> getCleanHits(std::vector<TCluster> clusters);
 

  
  // --------------
  void initRootTree();
  //void saveRootTree();
  void end();
  
protected:
  // plots
  TProfile *Efficiency_layer;
  TProfile *Multiplicity_layer;
  
  TProfile2D *m_padScan_mul;
  TProfile2D *m_padScan_mul_1;
  TProfile2D *m_padScan_mul_2;
  TProfile2D *m_padScan_mul_3;
  
  TProfile2D *m_padScan_eff;
  TProfile2D *m_padScan_eff_1;
  TProfile2D *m_padScan_eff_2;
  TProfile2D *m_padScan_eff_3;
  
  TH1F * Nhit_dist;
  TH1F * Ncluster_dist;
  TH1F * Ncluster_cut1;
  TH1F * Ncluster_cut2;
  TH1F * Ncluster_cut3;
  TH1F * Nhit_dist_clean;
  TH1F * Ncluster_dist_clean;
  
  TH1F * Nhit_dist_thr[3];
  //TH1F * Nhit_dist_thr[3];
  TH1F * mul_dist;
  
  TH1F * theta_dist;
  TH1F * phi_dist;
  TH2F * Cor_hit_chi2;
  TH2F * Cor_Ncluster_chi2;
  TH2F * Cor_theta_phi;
  TH1F * Nhit_Layer_01;// fo efficiencient 
  TH1F * Nhit_Layer_00;// not efficienct
  TH1F * hitMinDist_dist;
  TH1F * chi2_dist;
  TH1F * Prob_chi2;
  TH1F * Traveling_lenght_gas;
  TH1F * h_effAsicDist;
  TH1F * h_TrackLenght;
  TH1F * h_mul_1;
  TH1F * h_mul_2;
  TH1F * h_mul_3;
  
  TH1F * h_res_x;
  TH1F * h_res_y;
  
  TH1F * h_noise_stm;
  
  TH1F * h_res_x_dist[50];
  TH1F * h_res_y_dist[50];
  
  TH2F * h2_TrackLenght_angle;


  TProfile2D * m_eff_asic_map[50];
  TProfile2D * m_mul_asic_map[50];
  TProfile2D * m_eff_pad_map[50];

  
  

  TProfile * m_angle_mul_theta;
  TProfile * m_angle_mul_phi;
  TProfile * m_angle_Nhit_phi;
  TProfile * m_inc_angle_mul;
  TProfile * m_inc_angle_mul_1;
  TProfile * m_inc_angle_mul_2;
  TProfile * m_inc_angle_mul_3;
  TProfile * m_inc_angle_mu_1;
  TProfile * m_inc_angle_mu_2;
  TProfile * m_inc_angle_mu_3;
  TProfile * m_eff_inc_angle;
  TProfile2D * m_eff_inc_angle_layer;
  TProfile * m_layer_inc_angle;
  TProfile * m_Nhit_inc_angle;
  TProfile * m_chi2_residual;
  TProfile * m_angle_residual;
  TProfile * m_trav_len_eff;
  TProfile * m_angle_track_lenght;
  TProfile * m_res_layer_x;
  TProfile * m_res_layer_y;
  
  TProfile * m_asic_eff_mu;
  
  
  TH2F * h2_res_x;
  TH2F * h2_res_y;
  
  TProfile2D * m_threshold_pad_scan;
  TProfile2D * m_threshold_res;
  
  
  TH2F * h_res_1;
  TH2F * h_res_2;
  TH2F * h_res_3;
  
  TProfile * m_thr_pos_x;
  TProfile * m_thr_pos_y;
  
  TH2F * mu_dist_layer;
  TH2F * cor_mu_thr;
  TH1F * inc_angle;
  TH2F * res;
  TH2F * Cor_inc_angle_Nhit;
  TH2F * Cor_inc_angle_Chi2;

  TH2F * Cor_n1_n2;
  TH2F * Cor_n1_n3;
  TH2F * Cor_n2_n3;
  TH2F * mul_scan_stat;
  // cuts 
  int _NhitClusterCut;
  int _NClusterCut;
  int _NLayer;
  double _Chi2Cut;
  double _RegionCut;
  double _rpc_Gap;
  double _minDistCut;
  enum {maxn = 50};
  
  int Nrejected[maxn];
  int Naccepted[maxn];
  
  // tree braible
  int _NCluster;
  int _Chi2;
  double _inc_angle;
  //other 
  unsigned int _eventNr;  
  Int_t _Nhit; 
  int evtnum;
  
  // Parameters
  DHCAL_Mips *calo;
  
  std::string _treeFileName;
  std::string _treeName;
  
  std::string _hist_file_calo;
  std::string _hist_file_current;
  std::string _colName; 
  std::string _fileName;
  std::string _outFileName;
  std::string _offset_file;

  std::string _xml_geom_file;
  std::string _xml_scan_file;
  



  std::ostream *_output;
  bool _usePenetratedMuons;
  std::vector<std::string> _hcalCollections;
  int _overwrite;
  TTree *_outputTree;  
  TClonesArray *_hitArray;
  TClonesArray *_rootArray;// added by me 
  std::string _RecoParticleColl; 
  std::map<int,offset> _offsets;
  //std::map<int,double> _offset_x;
  LCWriter* _lcWriter;
};
#endif


