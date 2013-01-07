#include "DHCAL_Mips.h"
#include <marlin/tinyxml.h>

#define DEBUG false
 
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
  //_plate_positions.clear();
  _cluster.clear();
  _nplate =  50;
  _nhit_on_clus_cut = 4;
  _angle_cut = 3.14;
  _cut_region = 5.0;
  _nfired_layer = 0;
  _nfired_layer_test = false;
}

DHCAL_Mips::DHCAL_Mips()
{
  hist_book = new HistBooker();
  _plate_positions.clear();
  _nfired_layer = 0;
  for(int i = 0; i < 50; i++){
    //_plate_positions[i] = (i * 2.8);
    _plate_threshold_1[i] = get_threshold_value(1,160);
    _plate_threshold_2[i] = get_threshold_value(2,120);
    _plate_threshold_3[i] = get_threshold_value(3,120);
  }
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

bool DHCAL_Mips::set_geometry(std::string xml_geom_file){
  TiXmlDocument doc(xml_geom_file.c_str());
  if(doc.LoadFile())
    {
      TiXmlHandle hDoc(&doc);
      TiXmlElement *pRoot, *pParm;
      pRoot = doc.FirstChildElement("setup_geom");
      if(pRoot){
	pParm = pRoot->FirstChildElement("parameter");
	int i = 0; // for sorting the entries
	while(pParm)
	  {
	    string name = pParm->Attribute("name");
	    
	    if(name == "ChamberGeom") {
	      std::cout << " parameter name == " << name << std::endl;
	      
	      std::string value = pParm->GetText();
	      std::vector<std::string> lines;
	      istringstream iss(value);
	      copy(istream_iterator<string>(iss),
		   istream_iterator<string>(),
		   back_inserter<vector<string> >(lines));
	      
	      for(int iline = 0; iline < lines.size(); iline++){
		std::string line = lines.at(iline);
		//std::cout <<"\t line  == "<< line << std::endl;
		stringstream ss( line.c_str() );
		vector<string> result;
		int Layer_id;double position;
		while(ss.good()){
		  string substr;
		  getline( ss, substr, ',' );
		  result.push_back( substr );
		}
		istringstream ( result.at(0) ) >> Layer_id;
		istringstream ( result.at(3) ) >> position;
		//std::cout << " " << Layer_id <<" "<< position << std::endl;
		_plate_positions[Layer_id] = position;
		//std::cout << "maps set geometry "  <<  _plate_positions.size()<< std::endl;	      
}
	    }
	    pParm = pParm->NextSiblingElement("parameter");
	  }
	std::cout << "geom file : -----> " << std::endl;
	for(std::map<int, double >::iterator iter=_plate_positions.begin();iter!=_plate_positions.end();++iter)
	  {
	    std::cout << "plate == " << iter->first << "\t pos == " << iter->second << std::endl; 
	  }
	return true;
      }
    }
  else 
    {
      cerr << "Could not load XML File." << endl;
      return false;
    }  
}

void DHCAL_Mips::set_run_number(int run_id)
{
  _run_id = run_id;
}

bool DHCAL_Mips::read_xml_scan(std::string xmlfile){
  TiXmlDocument doc(xmlfile.c_str());
  //combo = (CComboBox*)GetDlgItem(IDC_EGC_CARD_TYPE);
  
  if(doc.LoadFile())
    {
      TiXmlHandle hDoc(&doc);
      TiXmlElement *pRoot, *pParm;
      pRoot = doc.FirstChildElement("scan_threshold");
      if(pRoot)
	{
	  pParm = pRoot->FirstChildElement("run");
	  int i = 0; // for sorting the entries
	  while(pParm)
	    {
	      int run_id; istringstream ( pParm->Attribute("value")) >> run_id ;
	      std::cout << " run value == " << run_id << std::endl;
	      // read the node content :
	      std::map<int, thr_t > scan_point;
	      
	      std::string value = pParm->GetText();
	      std::vector<std::string> lines;
	      istringstream iss(value);
	      copy(istream_iterator<string>(iss),
		   istream_iterator<string>(),
		   back_inserter<vector<string> >(lines));
	      
	      for(int iline = 0; iline < lines.size(); iline++){
		std::string line = lines.at(iline);
		std::cout <<"\t line  == "<< line << std::endl;
		stringstream ss( line.c_str() );
		vector<string> result;
		int Layer_id;thr_t threshold; 
		while(ss.good()){
		  string substr;
		  getline( ss, substr, ',' );
		  result.push_back( substr );
		}
		istringstream ( result.at(0) ) >> threshold.t_1;
		istringstream ( result.at(1) ) >> threshold.t_2;
		istringstream ( result.at(2) ) >> threshold.t_3;
		istringstream ( result.at(3) ) >> Layer_id;
		scan_point[Layer_id] = threshold;
	      }
	      
	      scan_all[run_id] = scan_point ;
	      pParm = pParm->NextSiblingElement("run");
	    }
	  return true;
	}
    }
  else 
    {
      cerr << "Could not load XML File." << endl;
      return false;
    }  
}


void DHCAL_Mips::set_xml_threshold(std::string filename)
{
  TiXmlDocument doc(filename.c_str());
  bool load_key = doc.LoadFile();
  
  if(load_key){
    std::cout << "file : " << filename.c_str() << std::endl;
    TiXmlHandle h_doc( &doc);
    TiXmlElement* p_elem;
    TiXmlHandle h_root(0);
    // load scan
    {
      p_elem = h_doc.FirstChildElement().Element();
      if(!p_elem) std::cout << "error ......" << std::endl;
      std::cout <<"scan node == " <<p_elem->Value() << std::endl;

      h_root = TiXmlHandle(p_elem);
    }
    // read run 
    {
      p_elem = h_root.FirstChild("run").Element();
      std::string key = p_elem->Attribute("value");
      std::cout <<"run node == " <<key.c_str() << std::endl;
      //p_elem = p_elem->NextSiblingElement(); // iteration 
    }
  }else{
    std::cout << "error .... " << std::endl;
  }

}


int DHCAL_Mips::get_nfired_layer(){
  std::vector<int> hits_plate(_nplate);//hits_plate.clear();
  for(std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++){
    hits_plate[c->getLayer_id()]++ ;
  }
  int n_tmp = 0;
  for(int pl = 1; pl < _nplate ; pl++){
    if( hits_plate[pl] > 0 ) n_tmp++ ;
  }
  //std::cout <<  "---> fired layer == " << n_tmp <<std::endl;
  return n_tmp;
}

bool DHCAL_Mips::penetrated_muon_test(){
  std::vector<int> hits_plate(_nplate);
  for(std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++){
    hits_plate[c->getLayer_id()]++ ;
  }
  int front_side = 0, back_side = 0,n_on_plate=0;
  for(int pl = 1; pl < 10 ; pl++){
    if( hits_plate[_nplate - pl ] > 0 ) back_side++ ;
    if( hits_plate[ 0      + pl ] > 0 ) front_side++ ;
  }
  
  for(int pl = 1; pl < _nplate ; pl++){
    if( hits_plate[pl] > 0 ) n_on_plate++ ;
  }
  
  hist_book->book_h1("/pntr_muon/n_on_layer",50,0,50)->Fill(n_on_plate);
  hist_book->book_h1("/pntr_muon/back_side",25,0,25)->Fill(back_side);
  hist_book->book_h1("/pntr_muon/front_side",25,0,25)->Fill(front_side);
  hist_book->book_h1("/pntr_muon/var_1",40,-20,20)->Fill(front_side - back_side );
  hist_book->book_h2("/pntr_muon/cor_front_back",25,25,0,25,0,25)->Fill(front_side,back_side);
  hist_book->get_h2("/pntr_muon/cor_front_back")->SetName(";front;back");
  hist_book->get_h1("/pntr_muon/var_1")->SetName(";(front+back)/nlayer;#");
  if(back_side > 0 && front_side > 0)
    return true;
  else 
    return false;
}

bool DHCAL_Mips::do_characterization(){
  
  std::vector<int> hits_plate(_nplate);//hits_plate.clear();
  for(std::vector<TCluster>::iterator c=_cluster.begin();c!=_cluster.end();c++){
    hits_plate[c->getLayer_id()]++ ;
  }
  int n_tmp = 0;
  for(int pl = 1; pl < _nplate ; pl++){
    if( hits_plate[pl] > 0 ) n_tmp++ ;
  }
  
  double density = _cluster.size()/(1.0*n_tmp);
  
  hist_book->book_h1("/char/nlayers",100,0,100)->Fill(n_tmp);
  hist_book->book_h1("/char/ncluster",100,0,100)->Fill(_cluster.size());
  hist_book->book_h1("/char/density",200,0,2)->Fill(density);
  
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
    TCluster tclus;
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
	tclus    = *c;
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
    std::map<int, thr_t> scan_threshold; scan_threshold.clear();
    if( scan_all.find(_run_id)->first == _run_id ){
      scan_threshold = scan_all.find(_run_id)->second;
    }
    double *p; p = track->getParameter();
    double eta   = get_incident_angle(p);
    // ===== fill some histo =======
    hist_book->book_h1("/prf/track/chi2",200,0,100)->Fill(track->getChi2());
    hist_book->book_h1("/prf/track/prob_chi2",200,0.,1.)->Fill(track->getProb());
    hist_book->book_h1("/prf/angle/eta",200,0.,Pi()/2.)->Fill(eta);
    if(cluster_size > 0 ){
      //scan part
      //if(scan_threshold.size() != 0){
	hist_book->book_profile1("/prf/scan/mu_scan_1",300,0.,0.6,0,4)->Fill(get_threshold_value(1,scan_threshold.find(6 )->second.t_1),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_1",300,0.,0.6,0,4)->Fill(get_threshold_value(1,scan_threshold.find(18)->second.t_1),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_1",300,0.,0.6,0,4)->Fill(get_threshold_value(1,scan_threshold.find(30)->second.t_1),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_2",300,0.,6. ,0,4)->Fill(get_threshold_value(2,scan_threshold.find(10)->second.t_2),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_2",300,0.,6. ,0,4)->Fill(get_threshold_value(2,scan_threshold.find(22)->second.t_2),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_2",300,0.,6. ,0,4)->Fill(get_threshold_value(2,scan_threshold.find(34)->second.t_2),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_3",500,0.,25.,0,4)->Fill(get_threshold_value(3,scan_threshold.find(10)->second.t_3),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_3",500,0.,25.,0,4)->Fill(get_threshold_value(3,scan_threshold.find(22)->second.t_3),cluster_size);
	hist_book->book_profile1("/prf/scan/mu_scan_3",500,0.,25.,0,4)->Fill(get_threshold_value(3,scan_threshold.find(34)->second.t_3),cluster_size);
	//}
      //
      
      hist_book->book_h1("/prf/mul/mu",20,0,20)->Fill(cluster_size);
      hist_book->book_h2("/prf/mul/mu_plate",_nplate,20,-0.5,(_nplate-1)+0.5,0,20)->Fill(plate,cluster_size);
      hist_book->book_profile1("/prf/mul/<mu>_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,4.)->Fill(plate,cluster_size);
      hist_book->book_profile1("/prf/mul/<mu>_angle",200,0.,Pi()/2.,0.,4.)->Fill(eta,cluster_size);
      hist_book->book_profile1("/prf/mul/<mu>_angle_1",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[0]+t_mul[1]+t_mul[2]);
      hist_book->book_profile1("/prf/mul/<mu>_angle_2",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[0]+t_mul[2]);
      hist_book->book_profile1("/prf/mul/<mu>_angle_3",200,0.,Pi()/2.,0.,4.)->Fill(eta,t_mul[2]);
    }
    // =============================
    
    if(( cluster_size <= 5)    && 
       ( chi2 > 0.0 )          &&  
       ( NclusterCut > 15 )    &&
       ( _cluster.size() > 5)
       ){
      
      //if(cluster_size > 0) 
      DHCAL_Mips::get_charge_shape(current_cluster,p);
      
      hist_book->book_h1("/prf/track/nhit_0",100,0,100)->Fill(t_mul[0]);
      hist_book->book_h1("/prf/track/nhit_1",100,0,100)->Fill(t_mul[1]);
      hist_book->book_h1("/prf/track/nhit_2",100,0,100)->Fill(t_mul[2]);
      hist_book->book_h1("/prf/track/nhit",100,0,100)->Fill(cluster_size);

      hist_book->book_h2("/prf/track/nhit_nnn",100,100,0,100,0,100)->Fill(cluster_size,t_mul[0]+t_mul[1]+t_mul[2] );
      
      //DHCAL_Mips::asic_map(plate,p,tclus);
      
      double zz = _plate_positions.find(plate)->second;
      //std::cout << "map size == " << _plate_positions.size() << std::endl;
      
      double xx =  p[0] + p[1] * zz ;
      double yy =  p[2] + p[3] * zz ;
      
      double dx = pos[0] - xx ;
      double dy = pos[1] - yy ;
      double dz = pos[2] - zz ;
#if DEBUG      
      std::cout <<"| "    <<  _plate_positions.find(plate)->first 
		<<"\t | " << xx <<" / " << pos[0] 
		<<"\t | " << yy <<" / " << pos[1] 
		<<"\t | " << zz <<" / " << pos[2] 
		<<"\t | " << std::endl;
#endif


      hist_book->book_h1("/prf/track/dx",100,-10,10)->Fill(dx);
      hist_book->book_h1("/prf/track/dy",100,-10,10)->Fill(dy);
      hist_book->book_h1("/prf/track/dz",100,-10,10)->Fill(dz);
      
      
      if( xx >= 0.  && xx <= 100 &&
	  yy >= 0.  && yy <= 100 ){
	double found_thr[3] ={0.,0.,0.};
	bool hit_on_plate = false;
	if( sqrt((dx*dx)+(dy*dy)) < 4.0 && fabs(dz) < 0.5 ){
	  hit_on_plate = true;
	
	  if(t_mul[1] != 0 || t_mul[0] != 0 || t_mul[2] != 0){
	    found_thr[0] = 1.;
	  }
	  if(t_mul[0] != 0 || t_mul[2] != 0){
	    found_thr[1] = 1.;
	  }
	  if(t_mul[2] != 0){
	    found_thr[2] = 1.;
	  }
	}
	//if( hit_on_plate ) 
	  // std::cout << "hit on plate ok ----> " << found_thr[0] << std::endl;
	
	
	
	hist_book->book_profile1("/prf/eff/eff_plate_1",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[0]);
	hist_book->book_profile1("/prf/eff/eff_plate_2",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[1]);
	hist_book->book_profile1("/prf/eff/eff_plate_3",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,found_thr[2]);
	
	hist_book->book_profile1("/prf/scan/eff_scan_1",300,0., 0.6,0,1)->Fill(get_threshold_value(1,scan_threshold.find(6 )->second.t_1),found_thr[0]);
	hist_book->book_profile1("/prf/scan/eff_scan_1",300,0., 0.6,0,1)->Fill(get_threshold_value(1,scan_threshold.find(18)->second.t_1),found_thr[0]);
	hist_book->book_profile1("/prf/scan/eff_scan_1",300,0., 0.6,0,1)->Fill(get_threshold_value(1,scan_threshold.find(30)->second.t_1),found_thr[0]);
	hist_book->book_profile1("/prf/scan/eff_scan_2",300,0., 6. ,0,1)->Fill(get_threshold_value(2,scan_threshold.find(10)->second.t_2),found_thr[1]);
	hist_book->book_profile1("/prf/scan/eff_scan_2",300,0., 6. ,0,1)->Fill(get_threshold_value(2,scan_threshold.find(22)->second.t_2),found_thr[1]);
	hist_book->book_profile1("/prf/scan/eff_scan_2",300,0., 6. ,0,1)->Fill(get_threshold_value(2,scan_threshold.find(34)->second.t_2),found_thr[1]);
	hist_book->book_profile1("/prf/scan/eff_scan_3",500,0., 25.,0,1)->Fill(get_threshold_value(3,scan_threshold.find(10)->second.t_3),found_thr[2]);
	hist_book->book_profile1("/prf/scan/eff_scan_3",500,0., 25.,0,1)->Fill(get_threshold_value(3,scan_threshold.find(22)->second.t_3),found_thr[2]);
	hist_book->book_profile1("/prf/scan/eff_scan_3",500,0., 25.,0,1)->Fill(get_threshold_value(3,scan_threshold.find(34)->second.t_3),found_thr[2]);
	
	if(hit_on_plate){
	  hist_book->book_profile1("/prf/eff/eff_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,1.);
	  hist_book->book_profile1("/prf/eff/eff_angle",100,0.,Pi()/2.,0.,1.)->Fill(eta,1.);
	  hist_book->book_profile1("/prf/eff/eff_chi2",200,0.,100.,0.,1.)->Fill(chi2,1.);
	}else{
	  hist_book->book_profile1("/prf/eff/eff_plate",_nplate,-0.5,(_nplate-1)+0.5,0.,1.)->Fill(plate,0.);
	  hist_book->book_profile1("/prf/eff/eff_angle",100,0.,Pi()/2.,0.,1.)->Fill(eta,0.);
	  hist_book->book_profile1("/prf/eff/eff_chi2",200,0.,100.,0.,1.)->Fill(chi2,0.);
	}
      }
    }
  }
  return true;
}

void DHCAL_Mips::save_plots(std::string file){
  hist_book->write_hitograms(file.c_str());
}

//void DHCAL_Mips::set_idDecoder_str(std::string str){
//  _idDecoderStr = str;
//}

void DHCAL_Mips::asic_map(int plate, double *p, TCluster clus){
  double *pos; pos = clus.getPosition();
  int I_0 = int(pos[0]/1.04125);
  int J_0 = int(pos[1]/1.04125);
  int I_r = int((p[0] +p[1] * _plate_positions.find(clus.getLayer_id())->second)/1.04125);
  int J_r = int((p[2] +p[3] * _plate_positions.find(clus.getLayer_id())->second)/1.04125);
  
  hist_book->book_profile2(Form("/scan/map/asic_mu_map_%i",plate),12,0,96,12,0,96,0.,4.)->Fill(I_0,J_0,clus.getN());
  
  if(fabs(I_0 - I_r) < 5 && fabs(J_0 - J_r) < 5){
    hist_book->book_profile2(Form("/scan/map/asic_eff_%i",plate),12,0,96,12,0,96,0.,1.)->Fill(I_0,J_0,1.);
  }else{
    hist_book->book_profile2(Form("/scan/map/asic_eff_%i",plate),12,0,96,12,0,96,0.,1.)->Fill(I_0,J_0,0.);
  }  
}

void DHCAL_Mips::get_charge_shape(CaloHitVec hits, double *p){
  int mu = hits.size();
#if DEBUG
  std::cout << "######################################" << std::endl;
  std::cout << "##### charge distribution studies ####" << std::endl;
  std::cout << "######################################" << std::endl;
#endif
  CellIDDecoder<EVENT::CalorimeterHit> idDecoder(_id_decoder_string.c_str());
  
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
      
      
      hist_book->book_profile2("/prf/mul/charge_shape",20,-0.5,0.5,20,-0.5,0.5,0.,4.)->Fill(Dx,Dy,mu);
      hist_book->book_profile1("/prf/mul/charge_radius",20,0.0,1.0,0.,0.4)->Fill(sqrt(Dx*Dx + Dy*Dy) , mu);
      
      
      if( mul_thr[0] > 0 || mul_thr[1] > 0 || mul_thr[2] > 0){
	hist_book->book_profile2("/prf/eff/pad_scan_eff_1",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("/prf/eff/pad_scan_eff_1",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      if( mul_thr[0] > 0 || mul_thr[2] > 0){
	hist_book->book_profile2("/prf/eff/pad_scan_eff_2",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("/prf/eff/pad_scan_eff_2",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      if( mul_thr[2] > 0){
	hist_book->book_profile2("/prf/eff/pad_scan_eff_3",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,1.);
      }else{
	hist_book->book_profile2("/prf/eff/pad_scan_eff_3",20,0.,1.,20,0.,1.,0.,1.)->Fill(rx,ry,0.);
      }
      
      hist_book->book_profile2("/prf/mul/pad_scan_mul",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mu);
      hist_book->book_profile2("/prf/mul/pad_scan_mul_1",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[0]+mul_thr[1]+mul_thr[2]);
      hist_book->book_profile2("/prf/mul/pad_scan_mul_2",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[0]+mul_thr[2]);
      hist_book->book_profile2("/prf/mul/pad_scan_mul_3",20,0.,1.,20,0.,1.,0.,4.)->Fill(rx,ry,mul_thr[2]);
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

double DHCAL_Mips::get_euclidean_dist(double* point1,double* point2)
{
  return sqrt(pow(point1[0]-point2[0],2)+
	      pow(point1[1]-point2[1],2)+
	      pow(point1[2]-point2[2],2));
}

void DHCAL_Mips::set_threshold(int plate, int  t_1, int t_2, int t_3)
{
  if(plate = -1){
    for(int iplate; iplate < _nplate; iplate++ ){
      _plate_threshold_1[iplate] = get_threshold_value(1,t_1);
      _plate_threshold_2[iplate] = get_threshold_value(2,t_2);
      _plate_threshold_3[iplate] = get_threshold_value(3,t_3);
      
    }
  }else if(plate != 0){
    _plate_threshold_1[plate] = get_threshold_value(1,t_1);
    _plate_threshold_2[plate] = get_threshold_value(2,t_2);
    _plate_threshold_3[plate] = get_threshold_value(3,t_3);
  }
  
}

double DHCAL_Mips::get_threshold_value(int threshold_id,int thr)
{
  if(threshold_id == 1) 
    return (1.e-3   )*(thr - 90.) / 0.7  ;
  else if(threshold_id == 2) 
    return (1.e-3   )*(thr - 98.) / 0.08 ;
  else if(threshold_id == 3) 
    return (1.      )*(thr - 90.) / 16.3 ;
  else 
    std::cout << "error [threshold doesn't exist !! ]"<< std::endl; 
}
