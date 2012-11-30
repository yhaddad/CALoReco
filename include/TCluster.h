// -*- C++ -*-

#ifndef TCLUSTER_H
#define TCLUSTER_H

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <EVENT/CalorimeterHit.h>
#include <map>
//namespace TRACK {
/** Clusterisation Algorithm for SDHCAL
 * 
 * @author Yacine Haddad
 * @version $Id: v0.1 $
 */
struct offset{
  double delta_x;
  double delta_y;
};

class TCluster{
public:
  TCluster(){
    _hits.clear();
    _Layer_id= 0;
    _position[0]=_posError[0] = 0;
    _position[1]=_posError[1] = 0;
      _position[2]=_posError[2] = 0;
      _x = _y = _z = 0;
  }
  virtual ~TCluster(){
    _hits.clear();
    _Layer_id=-1;
    _position[0]=_posError[0] = 0;
    _position[1]=_posError[1] = 0;
    _position[2]=_posError[2] = 0;
    _x = _y = _z = 0;
  }
  virtual void clear(){
    _hits.clear();
    _Layer_id=-1;
    _position[0]=_posError[0] = 0;
    _position[1]=_posError[1] = 0;
    _position[2]=_posError[2] = 0;
    _x = _y = _z = 0;
  }
  
  virtual bool       Append(EVENT::CalorimeterHit *hit);
  virtual int        distance(const EVENT::CalorimeterHit *h1,
			      const EVENT::CalorimeterHit *h2);
  virtual void       setHits(EVENT::CalorimeterHit *hit);
  virtual void       setOffset(double Delta_x=0.0,double Delta_y=0.0);
  virtual EVENT::CalorimeterHitVec getHits();
  virtual int        getN(){return _hits.size();};
  virtual double*    getPosition();
  virtual void       setUseTreshold(bool useTreshold = false);
  virtual void       Print();
  virtual double*    getPos_Error();
  virtual void       setCollString(const std::string initString);
  virtual int        getLayer_id();
  virtual int*       getThrNhit();
private:
  
  EVENT::CalorimeterHitVec _hits;
  double _x,_y,_z;
  std::string _initString;
  bool _useTreshold;
  int _Layer_id;
  double _position[3];
  double _posError[3];
  double _delta_x;
  double _delta_y;
  EVENT::CalorimeterHitVec _Event;
};


class Clustering{
public:
  Clustering(){
    clusters.clear();
  }
  
  Clustering(std::vector<EVENT::CalorimeterHit*> hit,std::string initString);
  virtual ~Clustering(){}
  void setOffsets();
  void clearOffsets();
  void setOffsets(double* Xof, double* Yof);
  void readOffsetFile(std::string offsetfile);
  std::vector<TCluster> getClusters();
  void Print();
private:
  std::map<int,offset> _offsets;
  std::vector<TCluster> clusters;
  
};

//}
#endif
