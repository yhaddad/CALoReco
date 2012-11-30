#include "DHCAL_Mips.h"

#define DEBUG true
 
using namespace std;
using namespace TRACK;
using namespace TMath;

DHCAL_Mips* DHCAL_Mips::_me = 0;
DHCAL_Mips* DHCAL_Mips::instance()
{
  if( _me == 0)
    _me = new DHCAL_Mips();
}

void DHCAL_Mips::clear()
{
  _plate_positions.clear();
  _cluster.clear();
  _nplate =  50;
  _nhit_on_clus_cut = 4;
  _angle_cut = 3.14;
  _cut_region = 5.0;
}

DHCAL_Mips::DHCAL_Mips()
{
  hist_book = new HistBooker();
  _plate_positions.clear();

  for(int i = 0; i < 50; i++)
    _plate_positions[i] = (i * 2.8);
}

void DHCAL_Mips::set_clusters(std::vector<TCluster> clus)
{
  _cluster = clus;
}

void DHCAL_Mips::set_nlayer(int nlayer)
{
  _nplate = nlayer;
}

void DHCAL_Mips::set_nhit_on_clus_cut(int cut)
{
  _nhit_on_clus_cut = cut;
}

void DHCAL_Mips::set_angle_cut(double cut)
{
  _angle_cut = cut;
}

void DHCAL_Mips::set_region_cut(double cut)
{
  _cut_region = cut;
}

double DHCAL_Mips::get_incident_angle(double *parameter)
{
  return fabs(TMath::ACos(1/sqrt(parameter[1]*parameter[1]+parameter[3]*parameter[3]+1)));
}

void DHCAL_Mips::set_geometry(std::map<int, double> geometry){
  _plate_positions = geometry;
}

bool DHCAL_Mips::do_calibration()
{
  TRACK::MuonTrack *track = new TRACK::MuonTrack();
  double pos[3] = {-1.,-1.,-1.};
  for(int plate = 1; plate < 51 /*_nplate*/; plate++){
    double cluster_size = 0;
    track->reset();
    
    int NclusterCut = 0;
    for (std::vector<TCluster>::iterator cl= _cluster.begin();cl!=_cluster.end();cl++){
      if(cl->getLayer_id() != plate){
	NclusterCut++;
      }
    }
    
    CaloHitVec current_cluster;current_cluster.clear();
    int t_mul[3]= {0,0,0};
    for(std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++){
      if(c->getLayer_id() !=  plate){
	track->setX(c->getPosition()[0],c->getPos_Error()[0]);
	track->setY(c->getPosition()[1],c->getPos_Error()[1]);
	track->setZ(c->getPosition()[2],c->getPos_Error()[2]);
      } else {
	pos[0] = c->getPosition()[0];
	pos[1] = c->getPosition()[1];
	pos[2] = c->getPosition()[2];
	cluster_size = c->getN();
	current_cluster = c->getHits();
	// to be fixed
	t_mul[0] = c->getThrNhit()[0] ;
	t_mul[1] = c->getThrNhit()[1] ;
	t_mul[2] = c->getThrNhit()[2] ;
      }
    }
    track->trackReco();
    double chi2 = track->getChi2();
#if DEBUG
    track->positionStatus();
#endif 
    double *p; p = track->getParameter();
    double eta   = get_incident_angle(p);
    // ===== fill some histo =======
    hist_book->book_h1("prf/track/chi2",200,0,100)->Fill(track->getChi2());
    hist_book->book_h1("prf/track/prob_chi2",200,0.,1.)->Fill(track->getProb());
    hist_book->book_h1("prf/angle/eta",200,0.,Pi()/2.)->Fill(eta);
    if(cluster_size > 0){
      hist_book->book_h1("prf/mul/mu",20,0,20)->Fill(cluster_size);
      hist_book->book_h2("prf/mul/mu_plate",_nplate,-0.5,(_nplate-1)+0.5,20,0,20)->Fill(plate,cluster_size);
      hist_book->book_profile1("prf/mul/<mu>_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,4.)->Fill(plate,cluster_size);
      hist_book->book_profile1("prf/mul/<mu>_angle",200,0.,Pi()/2.,0.,4.)->Fill(eta,cluster_size);
      hist_book->book_profile1("prf/mul/<mu>_angle_1",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[0]+t_mul[1]+t_mul[2]);
      hist_book->book_profile1("prf/mul/<mu>_angle_2",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[0]+t_mul[2]);
      hist_book->book_profile1("prf/mul/<mu>_angle_3",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[2]);
    }
    // =============================
    
    if(( cluster_size <= 5)    && 
       ( chi2 > 0.0 )                          &&  
       ( NclusterCut > 15 )                    &&
       ( _cluster.size() > 5)//                  && 
       //( eta < _angle_cut)
       ){
      
      DHCAL_Mips::get_charge_shape(current_cluster,p);
      
      //double zz = _plate_positions.find(plate)->second;
      double zz = plate * 2.8;
      double dx = pos[0] - p[0] - p[1] * zz;
      double dy = pos[1] - p[2] - p[3] * zz;
      
      bool hit_on_plate = false;
      int found_thr[3] ={0,0,0};
      if(((p[0] +p[1] *zz) >= 0. && (p[0] +p[1] *zz) <= 100) &&
	 ((p[2] +p[3] *zz) >= 0. && (p[2] +p[3] *zz) <= 100) ){
	
	if( sqrt((dx*dx)+(dy*dy)) < 5.0 /*_cut_region*/ ){
	  hit_on_plate = true;
	  
	  if(t_mul[1] > 0 || t_mul[0] > 0 || t_mul[2] > 0){
	    found_thr[0] = 1;
	  }
	  if(t_mul[1] > 0 || t_mul[2] > 0){
	    found_thr[1] = 1;
	  }
	  if(t_mul[2] > 0){
	    found_thr[2] = 1;
	  }
	}
	
	
	hist_book->book_profile1("prf/eff/eff_plate_1",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[0]);
	hist_book->book_profile1("prf/eff/eff_plate_2",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[1]);
	hist_book->book_profile1("prf/eff/eff_plate_3",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[2]);
	
	
	if(hit_on_plate){
	  hist_book->book_profile1("prf/eff/eff_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,1.);
	  hist_book->book_profile1("prf/eff/eff_angle",100,0.,Pi()/2.,0.,1.)->Fill(eta,1.);
	  hist_book->book_profile1("prf/eff/eff_chi2",200,0.,100.,0.,1.)->Fill(chi2,1.);
	  //std::cout << "hit found ... "<< std::endl;
	}else{
	  hist_book->book_profile1("prf/eff/eff_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,0.);
	  hist_book->book_profile1("prf/eff/eff_angle",100,0.,Pi()/2.,0.,1.)->Fill(eta,0.);
	  hist_book->book_profile1("prf/eff/eff_chi2",200,0.,100.,0.,1.)->Fill(chi2,0.);
	}
      }
    }
  }
  return true;
}

void DHCAL_Mips::save_plots(){
  hist_book->write_hitograms("calo_perfermance.root");
}


void DHCAL_Mips::asic_eff_map(){
  
  
}

void DHCAL_Mips::asic_mul_map(){
  
  
}

void DHCAL_Mips::get_charge_shape(CaloHitVec hits, double *p){
  int mu = hits.size();
#if DEBUG
  std::cout << " charge distribution studies ..... > " << std::endl;
#endif
  CellIDDecoder<CalorimeterHit> idDecoder("I:7,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7");
  
  int mul_thr[3]={0,0,0};
  for (std::vector<EVENT::CalorimeterHit*>::iterator it= hits.begin();it!= hits.end();it++)
    {
      if((*it)->getEnergy() == 1.0){
	mul_thr[0]++;
      }else if((*it)->getEnergy() == 2.0){
	mul_thr[1]++;
      }else{
	mul_thr[2]++;
      }
    }
  
  for(CaloHitVec::iterator hit = hits.begin(); hit != hits.end(); hit++){
    double xx = (*hit)->getPosition()[0]/10.;
    double yy = (*hit)->getPosition()[1]/10.;
    double zz = (*hit)->getPosition()[2]/10.;
    
    int layer_id = idDecoder(*hit)["K"];
    
    
    int thr = int ((*hit)->getEnergy());
    double Dx =  (xx -p[0] -p[1] *zz);
    double Dy =  (yy -p[2] -p[3] *zz);
    
    
    if(fabs(Dx) <= 0.5 && fabs(Dy) <= 0.5){
      double rx = Dx+0.5 - int(Dx+0.5);
      double ry = Dy+0.5 - int(Dy+0.5);
      
      if(rx < 0.) rx = rx + 1;
      if(ry < 0.) ry = ry + 1;
      
      std::cout << "i'm here ... " << std::endl;
      hist_book->book_profile2("prf/mul/charge_shape",20,-0.5,0.5,20,-0.5,0.5,0.,4.)->Fill(Dx,Dy,mu);
      hist_book->book_profile1("prf/mul/charge_radius",20,0.0,1.0,0.,0.4)->Fill(sqrt(Dx*Dx + Dy*Dy) , mu);
      
      
      if( mul_thr[0] > 0 || mul_thr[1] > 0 || mul_thr[2] > 0){
	hist_book->book_profile2("prf/eff/pad_scan_eff_1",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("prf/eff/pad_scan_eff_1",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      if( mul_thr[0] > 0 || mul_thr[2] > 0){
	hist_book->book_profile2("prf/eff/pad_scan_eff_2",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("prf/eff/pad_scan_eff_2",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      if( mul_thr[2] > 0){
	hist_book->book_profile2("prf/eff/pad_scan_eff_3",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("prf/eff/pad_scan_eff_3",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      hist_book->book_profile2("prf/mul/pad_scan_mul",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mu);
      hist_book->book_profile2("prf/mul/pad_scan_mul_1",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[0]+mul_thr[1]+mul_thr[2]);
      hist_book->book_profile2("prf/mul/pad_scan_mul_2",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[0]+mul_thr[2]);
      hist_book->book_profile2("prf/mul/pad_scan_mul_3",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[2]);
    }
  }
  
}


double DHCAL_Mips::get_track_lenght()
{
  double p_a[3] = {0,0,0}; double p_b[3] = {0,0,0};
  int  k_min = 50, k_max = 0;
  for (std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++)
    {
      k_min = min(k_min, c->getLayer_id());
      k_max = min(k_max, c->getLayer_id());
    }
  
  for (std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++)
    {
      if(c->getLayer_id() == k_min ) { 
	p_a[0] = c->getPosition()[0];
	p_a[1] = c->getPosition()[1];
	p_a[2] = c->getPosition()[2];
      }
      if(c->getLayer_id() == k_max ) { 
	p_b[0] = c->getPosition()[0];
	p_b[1] = c->getPosition()[1];
	p_b[2] = c->getPosition()[2];
      }
    }
  return DHCAL_Mips::get_euclidean_dist(p_a,p_b);
}

double DHCAL_Mips::get_euclidean_dist(double* point1,double* point2){
  return sqrt(pow(point1[0]-point2[0],2)+
	      pow(point1[1]-point2[1],2)+
	      pow(point1[2]-point2[2],2));
}

