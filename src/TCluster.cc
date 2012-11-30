/*
 * DHCAL neigbor clustering 
 * Haddad Yacine 
 * 01/2012 
 * LLR - Ecole Polytechnique 
 */

#include "TCluster.h"
#include <math.h>
#include <TMath.h>
#include "UTIL/CellIDDecoder.h"
#include <string>
#include <algorithm>

using namespace std;

double const LCell = 1.;
double const SigmaGap = 0.010421 ;
double const sigma = 1. / sqrt(12.);    
//namespace TRACK{
//____________________________________________________
void TCluster::setHits(EVENT::CalorimeterHit *hit){
  _Event.push_back(hit);
}

//____________________________________________________
int TCluster::distance(const EVENT::CalorimeterHit *h1,const EVENT::CalorimeterHit *h2){
  CellIDDecoder<CalorimeterHit> idDecoder(_initString);
  if(idDecoder(h1)["K"] != idDecoder(h2)["K"]) return -1;
  return fabs((idDecoder(h1)["I"]-idDecoder(h2)["I"])+(idDecoder(h1)["J"]-idDecoder(h2)["J"]))/10.;
}
//____________________________________________________
void TCluster::setOffset(double Delta_x,double Delta_y){
  _delta_x = Delta_x;
  _delta_x = Delta_x;
}
//____________________________________________________
bool TCluster::Append(EVENT::CalorimeterHit *hit){
  CellIDDecoder<CalorimeterHit> idDecoder(_initString);
  if( _hits.size()==0 ){
    _hits.push_back(hit);
    return true;
  }
  bool append = false;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it= _hits.begin();it!=_hits.end();it++){
    //if (idDecoder(hit)["K"] != idDecoder(*it)["K"] ) return false;
    //if (hit == (*it)) break; 
    if(TCluster::distance(hit,(*it)) < 2 && TCluster::distance(hit,(*it)) >= 0){
      append= true;
      break;
    }
  }
  if(append) _hits.push_back(hit);
  return append;
}
//____________________________________
void TCluster::Print()
{
  for (std::vector<EVENT::CalorimeterHit*>::iterator it= _hits.begin();it!=_hits.end();it++)
    {
      std::cout<<"\t"<< (*it)->getPosition()[0]<<"\t"<< (*it)->getPosition()[0]<<"\t"<< (*it)->getPosition()[0] <<std::endl; 
    }
}
//____________________________________________________
double* TCluster::getPosition(){
  double  min_x   = 105.;double  min_y  = 105.;
  double  max_x   = 0.  ;double  max_y  = 0.  ;
  double  sum_x   = 0.  ;double  sum_y  = 0.  ;
  double  sum2_x  = 0.  ;double  sum2_y = 0.  ;
  
  int  min_I   = 100 ;int  min_J  = 100;
  int  max_I   = 0   ;int  max_J  = 0  ;
  //double _errorx  = 0.  ;double _errory = 0.  ; 
  //int Nhit = _hits.size();
  int Nhit = 0;
  CellIDDecoder<CalorimeterHit> idDecoder(_initString);
  for (std::vector<EVENT::CalorimeterHit*>::iterator it= _hits.begin();it!=_hits.end();it++)
    { 
      int I = idDecoder(*it)["I"];
      int J = idDecoder(*it)["J"];
      
      if(min_x > (*it)->getPosition ()[0]/10.) min_x = (*it)->getPosition ()[0]/10.;
      if(min_y > (*it)->getPosition ()[1]/10.) min_y = (*it)->getPosition ()[1]/10.;
      if(max_x < (*it)->getPosition ()[0]/10.) max_x = (*it)->getPosition ()[0]/10.;
      if(max_y < (*it)->getPosition ()[1]/10.) max_y = (*it)->getPosition ()[1]/10.;
      
      if(min_I > I) min_I = I;
      if(min_J > J) min_J = J;
      if(max_I < I) max_I = I;
      if(max_J < J) max_J = J;
      
      if ( _useTreshold ){
	Nhit     += (*it)->getEnergy();
	
	sum_x    += (*it)->getEnergy() * (*it)->getPosition ()[0]/10.;
	sum_y    += (*it)->getEnergy() * (*it)->getPosition ()[1]/10.;
	sum2_x   +=  ((*it)->getPosition ()[0]/10.)*((*it)->getPosition ()[0]/10.);
	sum2_y   +=  ((*it)->getPosition ()[1]/10.)*((*it)->getPosition ()[1]/10.);
	
      }else{
	Nhit     += 1;
	
	sum_x    +=  (*it)->getPosition ()[0]/10.;
	sum_y    +=  (*it)->getPosition ()[1]/10.;
	sum2_x   +=  ((*it)->getPosition ()[0]/10.)*((*it)->getPosition ()[0]/10.);
	sum2_y   +=  ((*it)->getPosition ()[1]/10.)*((*it)->getPosition ()[1]/10.);
	
      }
      _position[2] =(*it)->getPosition ()[2]/10.;
    }
  
  double span_x =  fabs( max_x - min_x );
  double span_y =  fabs( max_y - min_y );
  
  int nhit_x = fabs( max_I - min_I ) +1;
  int nhit_y = fabs( max_J - min_J ) +1;
  
  double spread_x = sqrt((sum2_x*Nhit - sum_x*sum_x)/Nhit/Nhit);
  double spread_y = sqrt((sum2_y*Nhit - sum_y*sum_y)/Nhit/Nhit);
  
  if ( Nhit  > 0 ) {
    _position[0] =sum_x/double(Nhit) - _delta_x;
    _position[1] =sum_y/double(Nhit) - _delta_y;
  }
      
  _posError[0] =(LCell - 2 * SigmaGap) * sqrt(nhit_x/12.);
  _posError[1] =(LCell - 2 * SigmaGap) * sqrt(nhit_y/12.);
  //_posError[0] = 1.;
  //_posError[1] = 1.;
  _posError[2] = 0.0;
  
  return _position;
}

//____________________________________
double* TCluster::getPos_Error(){
  return _posError;
}
//____________________________________
void TCluster::setUseTreshold(bool useTreshold){
  _useTreshold = useTreshold;
}
//____________________________________
EVENT::CalorimeterHitVec TCluster::getHits(){ return _hits;}

//____________________________________
int* TCluster::getThrNhit(){ 
  int nhit[3] = {0,0,0};
  for (std::vector<EVENT::CalorimeterHit*>::iterator it= _hits.begin();it!=_hits.end();it++)
    {
      nhit[int((*it)->getEnergy()-1)]++;
    }
  return nhit;
}
//____________________________________
int TCluster:: getLayer_id(){
  //std::cout << "dDecoder(_initString) == "<< _initString  << std::endl;
  CellIDDecoder<CalorimeterHit> idDecoder(_initString);
  for (std::vector<EVENT::CalorimeterHit*>::iterator it= _hits.begin();it!=_hits.end();it++){
    _Layer_id = idDecoder(*it)["K"];
  }
  //std::cout << " Layer id == " <<  _Layer_id << std::endl;
  return _Layer_id;
}
//____________________________________
void TCluster::setCollString(const std::string initString){
  //_initString = initString ;
  //std::cout << "cellID decoder == "<<initString << std::endl;
  //_initString = "I:7,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7";
  _initString = initString ;
}
//____________________________________
Clustering::Clustering(std::vector<EVENT::CalorimeterHit*> hits, std::string initString){
  TCluster cl;
  cl.setCollString(initString);
  std::vector<EVENT::CalorimeterHit*> other_hit;
  
  for (std::vector<EVENT::CalorimeterHit*>::iterator ihit= hits.begin();ihit!=hits.end();ihit++){
    bool append_hit = cl.Append(*ihit);
    if(append_hit == false){
      other_hit.push_back(*ihit);
    }
  }
  clusters.push_back(cl);
  
  bool toutou =  false;
  if(other_hit.size() != 0) toutou = true;
  while(toutou){
    
    cl.clear();
    bool find_hit = true;
    
    while(find_hit){
      std::vector<EVENT::CalorimeterHit*> hit_mem = other_hit;
      other_hit.clear();
      int nb_hit = cl.getN();
      for (std::vector<EVENT::CalorimeterHit*>::iterator it= hit_mem.begin();it!=hit_mem.end();it++){
	bool test = cl.Append(*it);
	if(test == false) 
	  other_hit.push_back(*it);
      }
      if(cl.getN() > nb_hit) find_hit = true;
      else find_hit = false;
    }
    if(other_hit.size() != 0) 
      toutou = true;
    else
      toutou = false; 
    clusters.push_back(cl);
  }
}

//____________________________________
std::vector<TCluster> Clustering::getClusters(){
  return clusters;
}
//____________________________________________________
void Clustering::setOffsets(){
  for (std::vector<TCluster>::iterator cl= clusters.begin();cl!=clusters.end();cl++){
    for(std::map<int,offset>::iterator of = _offsets.begin();of!=_offsets.end();of++)     {
      if(of->first == cl->getLayer_id()){
	cl->setOffset(of->second.delta_x,of->second.delta_y);
      }
    }
  }
}

void Clustering::setOffsets(double* Xof, double* Yof){
  for(int ilayer=0;ilayer < 51; ilayer++){
    offset contenu;
    contenu.delta_x = Xof[ilayer];
    contenu.delta_y = Yof[ilayer];
    _offsets[ilayer] = contenu;
  }
}

void Clustering::clearOffsets(){
  for (std::vector<TCluster>::iterator cl= clusters.begin();cl!=clusters.end();cl++){
    cl->setOffset(0.0,0.0);
  }
}
//____________________________________________________
void Clustering::readOffsetFile(std::string offsetfile){
  cout << "read the offset file .."<< endl;
  
  offset contenu;
  ifstream file(offsetfile.c_str(), ios::in);
  if(file){
    while(!file.eof()){
      int Layer;
      char co;
      file >> Layer >> co
           >> contenu.delta_x >> co
           >> contenu.delta_y;
      _offsets[Layer] = contenu;
    }
    file.close();
  }
  else
    cerr << "ERROR ... offset file not correct !" << endl;
}
//____________________________________
void Clustering::Print(){
  int i = 0;
  for (std::vector<TCluster>::iterator cl= clusters.begin();cl!=clusters.end();cl++){
    std::cout<<"cl :: "<< i <<" Nhit=="<< cl->getN()
	     <<"\t"<< cl->getLayer_id()
	     <<"\t"<< cl->getPosition()[0] <<"+-"<<cl->getPos_Error()[0]
	     <<"\t"<< cl->getPosition()[1] <<"+-"<<cl->getPos_Error()[1]
	     <<"\t"<< cl->getPosition()[2] <<"+-"<<cl->getPos_Error()[2]
	     <<std::endl;
    i++;
  }
} 
//}
