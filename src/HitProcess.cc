#include <HitProcess.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Vertex.h>
#include <EVENT/Cluster.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>// for display a track 
#include <IMPL/VertexImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include "UTIL/CellIDDecoder.h"
#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>  
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TObject.h"
#include "TCanvas.h"
#include <TMinuit.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <fstream>
#include <algorithm>
#include "TRACK/MuonTrack.h"
#include "HistBooker.h"
#include "DHCAL_Mips.h"

TH1F *Hit_dist;
TH2F *Hit_distance;

HistBooker* histo_book;
using namespace std;
using namespace TRACK;
using namespace TMath;
//=========================================================

HitProcess a_HitProcess_instance;
double HitProcess::distance(double z,double x,double y,double *p){
  double d2 = pow((x-p[0]-p[1]*z),2)+pow((y-p[2]-p[3]*z),2);
  return sqrt(d2);
}
double HitProcess::distEucl(double* point1,double* point2){
  return sqrt(pow(point1[0]-point2[0],2)+
	      pow(point1[1]-point2[1],2)+
	      pow(point1[2]-point2[2],2));
}
HitProcess::HitProcess():Processor("HitProcess"),
			 _output(0),
			 _outputTree(0),
			 _rootArray(0),
			 _NCluster(0){
  
  
  _description = "DHCAL Analysis" ;
 
  std::vector<std::string> hcalCollections;    
  hcalCollections.push_back(std::string("SDHCAL_HIT"));
  hcalCollections.push_back(std::string("SDHCAL_calohit_new"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "HCALCollections" ,
			    "HCAL Collection Names" ,
			    _hcalCollections ,
			    hcalCollections);
  
  _overwrite=0;
  registerProcessorParameter( "OverwriteFile" , 
			      "overwrit a file if zero." ,
			      _overwrite ,
			      _overwrite);
  
  _RecoParticleColl = "MuonRecoParticle";
  registerProcessorParameter( "RecoParticleCollName" ,
			      "The name of the Recostructed Particle Collection" ,
			      _RecoParticleColl ,
			      _RecoParticleColl);

  _NhitClusterCut = 4;
  registerProcessorParameter( "NhitCluster" , 
			      "Cut on number of the hits on cluster" ,
			      _NhitClusterCut ,
			      _NhitClusterCut);
  _Chi2Cut= 2.0;
  registerProcessorParameter( "Chi2Cut" , 
			      "Cut on chi2 " ,
			      _Chi2Cut ,
			      _Chi2Cut);
  _NClusterCut= 1 ;
  registerProcessorParameter( "NCluster" , 
			      "Cut on nbr of cluster" ,
			      _NClusterCut ,
			      _NClusterCut);
  _NLayer= 6 ;
  registerProcessorParameter( "NLayer" , 
			      "Number of Layer" ,
			      _NLayer ,
			      _NLayer);
  
  _RegionCut = 10.0;
  registerProcessorParameter( "RegionCut" , 
			      "Around the reco hit" ,
			      _RegionCut ,
			      _RegionCut);

  _treeFileName="plots/plots_file.root";
  registerProcessorParameter( "TreeOutputFile" ,
                              "The name of the file to which the ROOT tree will be written" ,
                              _treeFileName ,
                              _treeFileName);


  _xml_scan_file = "scan.xml";
  registerProcessorParameter( "scan_file" ,
                              "file contain the threshold scan setup" ,
                              _xml_scan_file ,
                              _xml_scan_file );

  
  
  _xml_geom_file = "scan.xml";
  registerProcessorParameter( "geom_file" ,
                              "file contain the stup geometry" ,
                              _xml_geom_file ,
                              _xml_geom_file );
  
  _offset_file="offsetfile.txt";
  registerProcessorParameter( "offsetfile" ,
                              "offset file for allignement" ,
                              _offset_file,
                              _offset_file);

  _treeName="recoco";
  registerProcessorParameter( "TreeName" ,
                              "The name of the ROOT tree" ,
                              _treeName ,
                              _treeName);

  _hist_file_calo = "hist_file_calo.root";
  registerProcessorParameter( "hist_file_calo" ,
                              "histogram book file for calibration" ,
                              _hist_file_calo ,
                              _hist_file_calo);
  
  _hist_file_current = "hist_file_calo.root";
  registerProcessorParameter( "hist_file_current" ,
                              "histogram book file for process" ,
                              _hist_file_current ,
                              _hist_file_current);

  _treeName="recoco";
  registerProcessorParameter( "TreeName" ,
                              "The name of the ROOT tree" ,
                              _treeName ,
                              _treeName);
  
  
  _rpc_Gap=2.8;
  registerProcessorParameter( "rpcGap" ,
                              "rpc Gap (cm)" ,
                              _rpc_Gap ,
                              _rpc_Gap);
  _minDistCut=12.;
  registerProcessorParameter( "minDistCut" ,
                              "noise remove by minimal distance (cm)" ,
                              _minDistCut ,
                              _minDistCut);

  
  _usePenetratedMuons=true;
  registerProcessorParameter( "usePenetratedMuons" ,
                              "Use Penetrated Muons" ,
			      _usePenetratedMuons ,
			      _usePenetratedMuons);

}

void HitProcess::initRootTree(){
  TFile *tree_file=new TFile(_treeFileName.c_str(),"RECREATE");
  
  if (!tree_file->IsOpen()) {
    delete tree_file;
    tree_file=new TFile(_treeFileName.c_str(),"NEW");
  }
  _outputTree = new TTree("recoco","recoco");
  _outputTree->SetAutoSave(32*1024*1024);
  
  // === tree of varable ===
  
  _outputTree->Branch("Nhit", &_Nhit,"Nhit/I");
  _outputTree->Branch("inc_angle", &_inc_angle,"inc_angle/D");
  _outputTree->Branch("chi2", &_Chi2,"chi2/D");
  _outputTree->Branch("NCluster",&_NCluster,"NCluster/I");
  
  
  
}

void HitProcess::init() {
  
  calo = new DHCAL_Mips();
  calo->set_geometry(_xml_geom_file.c_str());
  calo->read_xml_scan(_xml_scan_file.c_str());
  
  //calo->set_xml_threshold("scan_threshold.xml");
  
  histo_book = new HistBooker();
  printParameters();
  initRootTree();
  for(int i=0;i<maxn;i++){
    Nrejected[i] = 0;
    Naccepted[i] = 0;
  }
}

//=========================================================
std::vector<TCluster> HitProcess::getCleanHits(std::vector<TCluster> clusters){
  std::vector<TCluster> clean_hits;clean_hits.clear();
  for (std::vector<TCluster>::iterator cl1= clusters.begin();cl1!=clusters.end();cl1++){
    //cor_mu_thr->Fill(cl1->getThrNhit(1),1);
    //cor_mu_thr->Fill(cl1->getThrNhit(2),2);
    //cor_mu_thr->Fill(cl1->getThrNhit(3),3);
    double min = 100;
    for (std::vector<TCluster>::iterator cl2= clusters.begin();cl2!=clusters.end();cl2++){
      double pdist = distEucl(cl1->getPosition(),cl2->getPosition());
      
      if(min > pdist && pdist != 0){
	min = pdist;
      }
    }
    
    //hitMinDist_dist->Fill(min);
    
    if(min < _minDistCut){
      clean_hits.push_back(*cl1);
    }
  }
  return clean_hits;
}


double HitProcess::getTrackLenght(std::vector<TCluster> cluster)
{
  double p_a[3] = {0,0,0};
  double p_b[3] = {0,0,0};
  int  k_min = 50, k_max = 0;
  
  for (std::vector<TCluster>::iterator c=cluster.begin();c!=cluster.end();c++){
    k_min = min(k_min, c->getLayer_id());
    k_max = min(k_max, c->getLayer_id());
  }
  
  for (std::vector<TCluster>::iterator c=cluster.begin();c!=cluster.end();c++){
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
  
  return distEucl(p_a,p_b);
}

//=========================================================
void HitProcess::processEvent( LCEvent * evtP ){	
  if (evtP){	
    try{
      
      _eventNr=evtP->getEventNumber();
      _Nhit=0;
      
      int run_id  = evtP-> getRunNumber();

      
      LCCollection *event_coll = evtP ->getCollection(_hcalCollections[0].c_str());
      int Nhit = event_coll->getNumberOfElements();
      
      LCCollection *RecoParticle_coll = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      
      LCFlagImpl flag;
      flag.setBit(LCIO::CHBIT_LONG);
      RecoParticle_coll->setFlag(flag.getFlag());
      
      MuonTrack *track = new MuonTrack();
      std::cout<<"event no == "<<evtP->getEventNumber()<<std::endl;
      // 
      
      // reconstructed muon track 
      ReconstructedParticleImpl * reconstructedMuon     = new ReconstructedParticleImpl();   
      ReconstructedParticleImpl * reconstructedMuon_two = new ReconstructedParticleImpl(); 
      
      // ==== corrections for allignement ====
      double Xof[51];
      double Yof[51];
      for(int s=0; s < 51; s++){
	//Xof[s] = 0.437802 - 0.0559247*s + 0.00117728*s*s;
	Xof[s] = 0.0;
	Yof[s] = 0.0;
      }
      // ==== do clustering ====  
      std::vector<EVENT::CalorimeterHit*> event;
      std::vector<TCluster> clus;
      std::vector<TCluster> clus2;clus.clear();
      int nhit_thr[3] = {0,0,0};
      int nnn=0;
      
      for (int j(0); j < Nhit; ++j) { // fil the clustering 
        CalorimeterHit *hit_cluster = 
          dynamic_cast<CalorimeterHit*>(event_coll->getElementAt( j ));
	event.push_back(hit_cluster);
      }
      
      //for(int i =0; i < 3; i++){
      //	Nhit_dist_thr[i]->Fill(nhit_thr[i]);
      //}
      
      std::string initString = event_coll->getParameters().getStringVal(  LCIO::CellIDEncoding ) ;
      Clustering cc(event,initString);      //cc.readOffsetFile(_offset_file.c_str());
      cc.clearOffsets();
      cc.setOffsets(Xof,Yof);
      //cc.setOffsets();
      histo_book->book_h1("h1/raw/nhit",200,0,200)->Fill(Nhit);
      clus = cc.getClusters();
      clus2= cc.getClusters();
      //cc.Print();      
      
      // ------------------------------
      // --------  cleaning -----------
      // ------------------------------
      //Ncluster_cut1->Fill(clus.size());
      histo_book->book_h1("h1/raw/h1_NCluster",200,0,200)->Fill(clus.size());
      // ---  step 1
      std::vector<TCluster> clean_one;clean_one.clear();
      clean_one = getCleanHits(clus);
      histo_book->book_h1("h1/clean/h1_NCluster",200,0,200)->Fill(clean_one.size());
      // --- step 2
      
      //--- find the layers with NCluster > 1
      int nclus_layer[50];
      for(int i=0; i < 50 ; i++)
	nclus_layer[i] = 0;
      for (std::vector<TCluster>::iterator cl= clean_one.begin();cl!=clean_one.end();cl++){
	nclus_layer[cl->getLayer_id()]++;
      }
      track->reset();
      
      //--- fit without layers with NCluster > 1
      LCCollectionVec* col_clean_one = 
	new LCCollectionVec(LCIO::CALORIMETERHIT);
      col_clean_one->setFlag(col_clean_one->getFlag()|( 1 << LCIO::RCHBIT_LONG));
      float newcol_pos[3]={0,0,0};
      
      int test = 0;
      for(int k=0;k<_NLayer;k++){
	for (std::vector<TCluster>::iterator cl=clean_one.begin();cl!=clean_one.end();cl++){
	  //std::cout << "nclus_layer[k] == "  << nclus_layer[k]<< std::endl;
	  if((cl->getLayer_id() == k) && (nclus_layer[k] < 2)){
	    //std::cout << " que se passe-t-il ici == "<< cl->getPosition()[0] << std::endl;
	    track->setX(cl->getPosition()[0],cl->getPos_Error()[0]);
	    track->setY(cl->getPosition()[1],cl->getPos_Error()[1]);
	    track->setZ(cl->getPosition()[2],cl->getPos_Error()[2]);
	    test++;
	    
	    CalorimeterHitImpl* clusters = new CalorimeterHitImpl();
	    newcol_pos[0] = 10 * cl->getPosition()[0];
	    newcol_pos[1] = 10 * cl->getPosition()[1];
	    newcol_pos[2] = 10 * cl->getPosition()[2];
	    clusters->setPosition(newcol_pos);
	    col_clean_one->addElement(clusters);
	  }
	}
      }
      //std::cout << "nclus ==" << test << std::endl;
      //track->positionStatus();
      track->trackReco();
      // fisrt fit
      double chi2 = track->getChi2();
      float StartPoint[3];
      double direction[3];
      
      direction[0]=track->getParameter()[1];
      direction[1]=track->getParameter()[3];
      direction[2]=0.;
      
      StartPoint[0] = 10. * (track->getParameter()[0]);
      StartPoint[1] = 10. * (track->getParameter()[2]);
      StartPoint[2] = 0.;
      
      reconstructedMuon->setMomentum(direction);
      reconstructedMuon->setEnergy(chi2);// chi2
      reconstructedMuon->setMass(_NLayer*_rpc_Gap*10.); // longueur du det 
      reconstructedMuon->setCharge( 0 );
      reconstructedMuon->setReferencePoint(StartPoint);
      RecoParticle_coll->addElement( reconstructedMuon );
            
      // --- step 3
      std::vector<TCluster> clean_two;
      
      for (std::vector<TCluster>::iterator it = clean_one.begin();it!=clean_one.end();it++){
	double dmin = 10;
	TCluster tmp_cluster; 
	if(nclus_layer[it->getLayer_id()]<2){ // we find juste one cluster in the layer
	  clean_two.push_back((*it));
	}else{ // we find more than 1 cluster in the layer
	  for(std::vector<TCluster>::iterator cl = clean_one.begin();cl!=clean_one.end();cl++){
	    if(it->getLayer_id() == cl->getLayer_id()){
	      double d = distance(cl->getPosition()[2],cl->getPosition()[0],cl->getPosition()[1],track->getParameter());
	      if(d < dmin && d > 0.0 ) {
		dmin = d;
		tmp_cluster = (*cl);
	      }
	    }
	  }
	  if( tmp_cluster.getN() != 0 ){
	    clean_two.push_back(tmp_cluster);
	  }
	}
      }
      //Ncluster_cut3->Fill(clean_two.size());
      track->reset();
      
      LCCollectionVec* col_clean_two = 
	new LCCollectionVec(LCIO::CALORIMETERHIT);
      col_clean_two->setFlag(col_clean_two->getFlag()|( 1 << LCIO::RCHBIT_LONG));
      //std::cout << "selection two == "<< clean_two.size() <<std::endl;
      if(clean_two.size() > 0){
	for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	  track->setX(cl->getPosition()[0],cl->getPos_Error()[0]);
	  track->setY(cl->getPosition()[1],cl->getPos_Error()[1]);
	  track->setZ(cl->getPosition()[2],cl->getPos_Error()[2]);
	  
	  CalorimeterHitImpl* clusters2 = new CalorimeterHitImpl();
	  newcol_pos[0] = 10 * cl->getPosition()[0];
	  newcol_pos[1] = 10 * cl->getPosition()[1];
	  newcol_pos[2] = 10 * cl->getPosition()[2];
	  clusters2->setPosition(newcol_pos);
	  col_clean_two->addElement(clusters2);
	  
	}
      }
      //track->positionStatus();
      track->trackReco();
      chi2 = track->getChi2();
      


      direction[0]=track->getParameter()[1];
      direction[1]=track->getParameter()[3];
      direction[2]=0.;
      
      StartPoint[0] = 10. * (track->getParameter()[0]);
      StartPoint[1] = 10. * (track->getParameter()[2]);
      StartPoint[2] = 0.;
      
      reconstructedMuon_two->setMomentum(direction);
      reconstructedMuon_two->setEnergy(chi2);// chi2
      reconstructedMuon_two->setMass(_NLayer*_rpc_Gap*10.); // longueur du det 
      reconstructedMuon_two->setCharge( 0 );
      reconstructedMuon_two->setReferencePoint(StartPoint);
      RecoParticle_coll->addElement( reconstructedMuon_two );
      
      evtP->addCollection(RecoParticle_coll,_RecoParticleColl);
      evtP->addCollection(col_clean_one,"selection_1"); 
      evtP->addCollection(col_clean_two,"selection_2");
      // ----
      int n_hit = 0,n_1=0,n_2=0,n_3=0;
      for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	n_hit += cl->getN();
	//n_1   += cl->getThrNhit(1);
	//n_2   += cl->getThrNhit(2);
	//n_3   += cl->getThrNhit(3);
      }
      //Nhit_dist_clean->Fill(n_hit);
      histo_book->book_h1("h1/clean/Nhit",200,0,200)->Fill(n_hit);
      histo_book->book_h1("h1/clean/DAQ1",200,0,200)->Fill(n_1);
      histo_book->book_h1("h1/clean/DAQ2",200,0,200)->Fill(n_2);
      histo_book->book_h1("h1/clean/DAQ3",200,0,200)->Fill(n_3);


      histo_book->book_h2("h2/clean/N1_N2",200,200,0,200,0,200)->Fill(n_1,n_2);
      histo_book->book_h2("h2/clean/N1_N3",200,200,0,200,0,200)->Fill(n_1,n_3);
      histo_book->book_h2("h2/clean/N2_N3",200,200,0,200,0,200)->Fill(n_2,n_3);
            
      histo_book->book_h1("h1/clean2/NCluster",200,0,200)->Fill(clean_two.size());
      // Ncluster_dist_clean->Fill(clean_two.size());
      

      double ax = track->getParameter()[1];
      double ay = track->getParameter()[3];
      if(clean_two.size() > 0 ){
	_inc_angle = fabs(ACos(1/sqrt(ax*ax+ay*ay+1)));
	histo_book->book_profile1("m/angle/eta_nhit",400,0,Pi()/2.,0,200)->Fill(_inc_angle,Nhit);
	histo_book->book_profile1("m/angle/eta_chi2",400,0,Pi()/2.,0,200)->Fill(_inc_angle,chi2);
	
	histo_book->book_h2("h2/clean/chi2_nhit",200,200,0,200,0,200)->Fill(chi2,Nhit);
	histo_book->book_h2("h2/clean/chi2_clus",200,200,0,200,0,200)->Fill(chi2,Nhit);
	
	histo_book->book_h1("h1/track/chi2",200,0,200)->Fill(chi2);
	histo_book->book_h1("h1/track/prob_chi2",200,0,200)->Fill(track->getProb());
      } 
      
      // ALLIGEMENT BY FIXING TWO PLATE 9 and 39
      //MuonTrack *track_ali = new MuonTrack();
      //for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
      //	if(cl->getLayer_id() == 9 || cl->getLayer_id() == 39){
      //	  track_ali->setX(cl->getPosition()[0],cl->getPos_Error()[0]);
      //	  track_ali->setY(cl->getPosition()[1],cl->getPos_Error()[1]);
      //	  track_ali->setZ(cl->getPosition()[2],cl->getPos_Error()[2]);
      //	}
      //}
      //track_ali->trackReco();
      for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
      	double x_track = track->getParameter()[0] + track->getParameter()[1] * cl->getPosition()[2];
      	double y_track = track->getParameter()[2] + track->getParameter()[3] * cl->getPosition()[2];
      	double rx  = cl->getPosition()[0] - x_track;
      	double ry  = cl->getPosition()[1] - y_track;
      	
      	//h_res_x->Fill(rx);     
      	//h_res_y->Fill(ry);     
      	
      	//m_res_layer_x->Fill(cl->getLayer_id(),rx);
      	//m_res_layer_y->Fill(cl->getLayer_id(),ry);
      	//h_res_x_dist[cl->getLayer_id()]->Fill(rx);
      	//h_res_y_dist[cl->getLayer_id()]->Fill(ry);
      	//h2_res_x->Fill(cl->getLayer_id(),rx);
      	//h2_res_y->Fill(cl->getLayer_id(),ry);
      }
      
      //delete track_ali;
      
      //////////////////////////////////////////////////
      //////////////////////////////////////////////////
      // juste for residual biased 
      for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	double x_track = track->getParameter()[0] + track->getParameter()[1] * cl->getPosition()[2];
	double residu  = cl->getPosition()[0] - x_track;
	//m_chi2_residual->Fill(chi2,residu);
	//m_angle_residual->Fill(_inc_angle,residu);
      }
      
      
      

      
      
      /* noise estimation part 
	 the idea is to calculate the diffrance beetween 
	 the cleaned hit and the raw hits.
	 
	 Asic clock stamp is : 200e-09 s 
	 the time window is  : 600e-09 s
      */
      
      double noise_rate = (Nhit - n_hit)/600.0 ; // in MHz
      //h_noise_stm -> Fill(noise_rate);
      
      // efficiency part 
      //* - * /
      bool muons_pet_test;
      if( _usePenetratedMuons )
	muons_pet_test = penetratedMuons(clean_two,5,2);
      else
	muons_pet_test = true;
      
      //== DHCAL_Mips test ====
      
            

      
      // ======================
      if(chi2 > 0. 
	 && chi2 < _Chi2Cut &&
	 _inc_angle < 0.2){ // muon beam selection 
	double tl = getTrackLenght(clean_two);
	//h_TrackLenght->Fill(tl);
	//h2_TrackLenght_angle->Fill(tl,_inc_angle);
	//m_angle_track_lenght->Fill(_inc_angle,tl);
	// -- slecet a number of cluster
	//
	//	int NclusCut=0;
	//	for(int k=0;k<_NLayer;k++){
	//	  int cont = 0;
	//	  for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	//	    if(cl->getLayer_id() == k) cont++; 
	//	  }
	//	  if(cont++ > 0) NclusCut++;
	//	}
	//	if(NclusCut > _NlayerCut) 
	//--------------------------------------------------------------------------------
	// Fill the tree 
      
	_Nhit = Nhit;
	_NCluster = clean_two.size();
	_Chi2 = track->getChi2();
	double ax = track->getParameter()[1];
	double ay = track->getParameter()[3];
	_inc_angle = fabs(ACos(1/sqrt((ax*ax)+(ay*ay)+1)));
	//m_Nhit_inc_angle->Fill(_inc_angle,Nhit);
	_outputTree->Fill();
      
	//-----------------------------------ù ---------------------------------------------
	//Ncluster_dist->Fill(clean_two.size());
	//double ax = track->getParameter()[1];
	//double ay = track->getParameter()[3];
	//inc_angle->Fill(_inc_angle);
	//Traveling_lenght_gas->Fill(2./Cos(_inc_angle));
	//Cor_inc_angle_Nhit->Fill(Nhit,fabs(ACos(1/sqrt(ax*ax+ay*ay+1))));
      
	// -------------------------------------------------------------------------------
      
	calo->clear();
	calo->set_id_decoder_string( event_coll->getParameters().getStringVal(  LCIO::CellIDEncoding )) ;
	calo->set_run_number(run_id);
	calo->set_nlayer(50);
	calo->set_clusters(clean_two);
	calo->set_nhit_on_clus_cut(5);
	calo->set_angle_cut(3.14);
	calo->set_region_cut(_RegionCut);
	if(calo->penetrated_muon_test()){
	  calo->do_calibration();
	  calo->do_characterization();
	}

	//add the condition of number of the hits here 
	//for(int ilayer=0; ilayer < _NLayer; ilayer++){ // loop over the layer 
	//  // cut on number of fired layer
	//  int NclusterCut = 0;
	//  for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	//    if(cl->getLayer_id() != ilayer){
	//      NclusterCut++;
	//    }
	//  }
	//  if(NclusterCut > 15 && clean_two.size() > 5){
	//    double pos[3] = {-1.,-1.,-1.};// the target cluster position il layer i
	//    int Multiplicity = 0;
	//    int mul_thr[3]={0,0,0};
	//    chi2 = -1;
	//    track->reset();
	//    std::cout<<"*********"<<std::endl;
	//    std::vector<EVENT::CalorimeterHit*> current_cluster;
	//    //for (std::vector<TCluster>::iterator cl= clus.begin();cl!=clus.end();cl++){
	//    for (std::vector<TCluster>::iterator cl= clean_two.begin();cl!=clean_two.end();cl++){
	//      //if(((_NLayer) - cl->getLayer_id()) != ilayer){cleanned_events
	//      if(cl->getLayer_id() != ilayer){
	//	//std::cout<<" layer ---- "<< cl->getLayer_id() << " ---> "<< ilayer <<std::endl;
	//	track->setX(cl->getPosition()[0],cl->getPos_Error()[0]);
	//	track->setY(cl->getPosition()[1],cl->getPos_Error()[1]);
	//	track->setZ(cl->getPosition()[2],cl->getPos_Error()[2]);
	//	//track->setZ((_Nlayer*_rpc_Gap) - cl->getPosition()[2],cl->getPos_Error()[2]);
	//      }else{
	//	pos[0]=cl->getPosition()[0];
	//	pos[1]=cl->getPosition()[1];
	//	//std::cout<< "getPosition()[2] == " << cl->getPosition()[2] << std::endl;
	//	//pos[2]=(_Nlayer*_rpc_Gap) - cl->getPosition()[2];
	//	pos[2]=cl->getPosition()[2];
	//	//std::cout<< "pos[2] == " << pos[2]<< std::endl;
	//	Multiplicity = cl->getN();
	//	mul_thr[0] = cl->getThrNhit(1);
	//	mul_thr[1] = cl->getThrNhit(2);
	//	mul_thr[2] = cl->getThrNhit(3);
	//	current_cluster = cl->getHits();
	//      }
	//    }
	//  
	//    //track->positionStatus();
	//    track->trackReco();
	//    chi2 = track->getChi2();
	//    // fill 
	//    Cor_theta_phi->Fill(TMath::ATan(track->getParameter()[1]),
	//			TMath::ATan(track->getParameter()[3]));
	//  
	//    theta_dist->Fill(TMath::ATan(track->getParameter()[1]));
	//    phi_dist->Fill(TMath::ATan(track->getParameter()[3]));
	//  
	//
	//    // ======================
	//    double *p;
	//    p = track -> getParameter(); 
	//  
	//    // ======================
	//
	//    if((Multiplicity <= _NhitClusterCut ) && (chi2 > 0)){
	//      //Cor_hit_chi2->Fill(chi2,Nhit);
	//      //chi2_dist->Fill(chi2);
	//      mu_dist_layer->Fill(ilayer,Multiplicity);
	//      double aax = track->getParameter()[1];
	//      double aay = track->getParameter()[3];
	//    
	//      if(fabs(ACos(1/sqrt(aax*aax+aay*aay+1))) < 1.5){ // cut on the angle
	//	fillCallibration(ilayer,_RegionCut,pos,track->getParameter());
	//	asicEffMap(ilayer,Multiplicity,pos,track->getParameter());
	//	effMap(current_cluster,track->getParameter());
	//	padEffMap(ilayer,pos,track->getParameter());
	//	//padEffScan();
	//      }
	//      if(Multiplicity>0){
	//	//double aax = track->getParameter()[1];
	//	//double aay = track->getParameter()[3];
	//	m_inc_angle_mul_1->Fill(fabs(ACos(1/sqrt(aax*aax+aay*aay+1))),mul_thr[0]);	
	//	m_inc_angle_mul_2->Fill(fabs(ACos(1/sqrt(aax*aax+aay*aay+1))),mul_thr[1]);	
	//	m_inc_angle_mul_3->Fill(fabs(ACos(1/sqrt(aax*aax+aay*aay+1))),mul_thr[2]);	
	//	m_inc_angle_mul->Fill(fabs(ACos(1/sqrt(aax*aax+aay*aay+1))),Multiplicity);	
	//	m_angle_mul_theta->Fill(fabs(TMath::ATan(track->getParameter()[1])),Multiplicity);
	//	m_angle_mul_phi->Fill(fabs(TMath::ATan(track->getParameter()[3])),Multiplicity);
	//
	//	Multiplicity_layer->Fill(ilayer,Multiplicity); 
	//	mul_dist->Fill(Multiplicity);
	//
	//	h_mul_1->Fill(mul_thr[0]);
	//	h_mul_2->Fill(mul_thr[1]);
	//	h_mul_3->Fill(mul_thr[2]);
	//      }
	//    }
	//    track->reset();
	//  }
      }
      //// loop over layer 
      //setReturnValue(true);  
      //}else {
      //setReturnValue(false);
      //
      // * - */
      //std::cout << " test 1"<<std::endl;
      delete track;
      //std::cout << " test 2"<<std::endl;
    }	
    catch (lcio::DataNotAvailableException err) { }
    //std::cout << "hello 99"<< std::endl;
  }
}	
//=========================================================
void HitProcess::end(){
  calo->save_plots(_hist_file_calo.c_str());
  histo_book->write_hitograms( _hist_file_current.c_str());
  
  if (_outputTree) {
    TFile *tree_file = _outputTree->GetCurrentFile(); 
    tree_file->Write();
    delete tree_file;
  }
}
 
//void  HitProcess::padEffMap(int Layer,double *position, double* p){
//  int I0 = int(position[0]/1.04125);
//  int J0 = int(position[1]/1.04125);
//
//  int It = int((p[0] +p[1] *Layer *_rpc_Gap)/1.04125);
//  int Jt = int((p[2] +p[3] *Layer *_rpc_Gap)/1.04125);
//  
//  if(fabs(I0-It)< 2 && fabs(J0-Jt)< 2){
//    m_eff_pad_map[Layer]->Fill(It,Jt,1.);
//  }else{
//    m_eff_pad_map[Layer]->Fill(It,Jt,0.);
//  }
//    
//}
//
//void  HitProcess::effMap(std::vector<EVENT::CalorimeterHit*> hits
//			 ,double* p){
//  int mul=hits.size();
//  CellIDDecoder<CalorimeterHit> idDecoder("I:7,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7");
//  int mul_thr[3]={0,0,0};
//  for (std::vector<EVENT::CalorimeterHit*>::iterator it= hits.begin();it!= hits.end();it++)
//    {
//      if((*it)->getEnergy() == 1.0){
//	mul_thr[0]++;
//      }else if((*it)->getEnergy() == 2.0){
//	mul_thr[1]++;
//      }else{
//	mul_thr[2]++;
//      }
//    }
//  
//  for (std::vector<EVENT::CalorimeterHit*>::iterator it= hits.begin();it!= hits.end();it++)
//    {
//      double xx = (*it)->getPosition()[0]/10.;
//      double yy = (*it)->getPosition()[1]/10.;
//      double zz = (*it)->getPosition()[2]/10.;
//      int layer_id = idDecoder(*it)["K"];
//      
//      int thr = int ((*it)->getEnergy());
//      double Dx =  (xx -p[0] -p[1] *zz);
//      double Dy =  (yy -p[2] -p[3] *zz);
//      
//      
//      //if(fabs(Dx) < 1. && fabs(Dy) < 1.){
//      //h_res_x->Fill(Dx);     
//      //h_res_y->Fill(Dy);     
//	
//      //m_res_layer_x->Fill(layer_id,Dx);
//      //m_res_layer_y->Fill(layer_id,Dy);
//      //h_res_x_dist[layer_id]->Fill(Dx);
//      //h_res_y_dist[layer_id]->Fill(Dy);
//      //h2_res_x->Fill(layer_id,Dx);
//      //h2_res_y->Fill(layer_id,Dy);
//	
//      //m_thr_pos_x->Fill(Dx,thr);
//      //m_thr_pos_y->Fill(Dy,thr);
//      //}
//      m_threshold_res->Fill(Dx,Dy,thr);
//      if(thr == 1)
//	h_res_1->Fill(Dx,Dy);
//      else if(thr == 2)
//	h_res_2->Fill(Dx,Dy);
//      else
//	h_res_3->Fill(Dx,Dy);
//      
//      res->Fill(Dx,Dy);
//      
//      if(fabs(Dx) <= 0.5 && fabs(Dy) <= 0.5){
//	double rx = Dx+0.5 - int(Dx+0.5);
//	double ry = Dy+0.5 - int(Dy+0.5);
//	
//	if(rx < 0.) rx = rx + 1;
//	if(ry < 0.) ry = ry + 1;
//	
//	m_padScan_mul->Fill(rx,ry,mul);
//	
//	if( mul_thr[0] > 0 || mul_thr[1] > 0 || mul_thr[2] > 0){
//	  m_padScan_eff_1->Fill(rx,ry,1.); 
//	}else{
//	  m_padScan_eff_1->Fill(rx,ry,0.);
//	}
//	
//	if( mul_thr[0] > 0 || mul_thr[2] > 0){
//	  m_padScan_eff_2->Fill(rx,ry,1.); 
//	}else{
//	  m_padScan_eff_2->Fill(rx,ry,0.);
//	}
//	
//	if( mul_thr[2] > 0){
//	  m_padScan_eff_3->Fill(rx,ry,1.); 
//	}else{
//	  m_padScan_eff_3->Fill(rx,ry,0.);
//	}
//	
//	m_padScan_mul_1->Fill(rx,ry,mul_thr[0]+mul_thr[1]+mul_thr[2]);
//	m_padScan_mul_2->Fill(rx,ry,mul_thr[0]+mul_thr[2]);
//	m_padScan_mul_3->Fill(rx,ry,mul_thr[2]);
//	mul_scan_stat->Fill(rx,ry);
//	
//	m_threshold_pad_scan->Fill(rx,ry,thr);
//      }
//    }
//}
////=========================================================
//
//bool HitProcess::penetratedMuons(std::vector<TCluster> cluster, int minlayer, int minclus){
//  int clus_on_layer[50]={0};
//  int nclus_in = 0;
//  int nclus_out= 0;
//  for (std::vector<TCluster>::iterator cl =cluster.begin();cl!=cluster.end();cl++){
//    clus_on_layer[cl->getLayer_id()]++;
//    if(cl->getLayer_id() < minlayer )
//      nclus_in++;
//    if(cl->getLayer_id() > 50-minlayer )
//      nclus_out++;
//  }
//  if( nclus_in >= minclus &&  nclus_out >= minclus)
//    return true;
//  else 
//    return false;
//}
//
////=========================================================
