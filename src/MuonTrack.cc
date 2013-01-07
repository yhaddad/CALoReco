
#include "TRACK/MuonTrack.h"

#include <iostream>
#include <assert.h>
#include <vector>
#include <math.h> 
#include <algorithm>
#include <TSystem.h>
#include <TMath.h>
TGraph2DErrors *gr;
using std::cout;
using std::endl;
namespace TRACK{
  //__________________________________________________________________
  MuonTrack::MuonTrack()
  {
    
    _Parameter[0]=50;
    _Parameter[1]=0;
    _Parameter[2]=50;
    _Parameter[3]=0;
    _Parameter[4]=0;
    _x.clear();
    _y.clear();
    _z.clear();
    _ex.clear();
    _ey.clear();
    _ez.clear();
  }
  
  //__________________________________________________________________
  double MuonTrack::distance2(double z,double x,double y,double ex,double ey,double *p) { 
    double d2 = pow((x-p[0]-p[1]*z)/ex,2)+pow((y-p[2]-p[3]*z)/ey,2); 
    return d2; 
  }
  
  //__________________________________________________________________
  void MuonTrack::linearMinFunction(int &, double *, double & sum, double * par, int) { 
    assert(gr != 0);
    double * x = gr->GetX();
    double * y = gr->GetY();
    double * z = gr->GetZ();
    
    double * ex = gr->GetEX();
    double * ey = gr->GetEY();
    double * ez = gr->GetEZ();
    
    int npoints = gr->GetN();
    sum = 0;
    //double sigma = 10./sqrt(12);
    for (int i  = 0; i < npoints; ++i) { 
      double d = distance2(z[i],x[i],y[i],ex[i],ey[i],par);
      sum += d;
    }
    sum = sum/double(npoints - 4);
  }
  
  //__________________________________________________________________
  void MuonTrack::trackFitter(double *p){
    
    //! this fuction is called to fit the line and return somme parameter of line fit
    //! the minimisation of sumfunction
    Double_t chi2;
    int ier;
    Double_t arglist[10];
    double * test;
    _Parameter[4]=0;
    arglist[0] = 3; 
    gSystem->Load("libMinuit");
    
    
    TMinuit *min = new TMinuit(4);
    
    min->SetPrintLevel(-1);
    
    min->SetFCN(linearMinFunction);
    min->mnparm(0,"p0",p[0],0.5,0,0, ier);
    min->mnparm(1,"p1",p[1],0.5,0,0, ier);
    min->mnparm(2,"p2",p[2],0.5,0,0, ier);
    min->mnparm(3,"p3",p[3],0.5,0,0, ier);
  
    arglist[0] = 5000; // number of function calls 
    arglist[1] = 0.001; // tolerance 
  
    min->mnexcm("MIGRAD",arglist,2, ier);
  
    double parFit[4];
    double errFit[4];
  
    for (int i = 0; i <4; ++i){
      min->GetParameter(i, parFit[i], errFit[i]);
    }
    min->Eval(4, test, chi2, parFit, ier);
  
    for(int k=0;k < 4;k++){
      p[k] = parFit[k];
#if DEBUG
      std::cout<<"Prameter fit p["<<k<<"] == "<< parFit[k] << std::endl;
#endif
    }
    p[4] = chi2;
    delete min;
  }
  //__________________________________________________________________
  void MuonTrack::trackReco(){
    
    int Nhit=  _z.size();
    gr=new TGraph2DErrors();
    for(int k=0; k<Nhit; k++){
      gr->SetPoint(k,_x[k],_y[k],_z[k]);
      gr->SetPointError(k,_ex[k],_ey[k],_ez[k]);
    }
    
    trackFitter(_Parameter);
#if DEBUG    
    std::cout<<"Chi2 == "<<_Parameter[4]<<std::endl;
#endif 
    gr->Delete();
  }
  
  //__________________________________________________________________
  void MuonTrack::setX(double x,double ex){_x.push_back(x);_ex.push_back(ex);}
  void MuonTrack::setY(double y,double ey){_y.push_back(y);_ey.push_back(ey);}
  void MuonTrack::setZ(double z,double ez){_z.push_back(z);_ez.push_back(ez);}
  
  //__________________________________________________________________
  double MuonTrack::getChi2(){return _Parameter[4];}  
  //__________________________________________________________________
  void MuonTrack::reset(){
    _x.clear();_y.clear();_z.clear();
    _ex.clear();_ey.clear();_ez.clear();
  }
  
  //__________________________________________________________________
  double* MuonTrack::getParameter(){
    return _Parameter;
  }
  
  //__________________________________________________________________
  void MuonTrack::positionStatus(){
    std::cout<<"position status ...."<<std::endl;
    for(int i=0;i<_z.size();i++){
      std::cout<<_x[i]<<"+-"<<_ex[i]<<" .. "
	       <<_y[i]<<"+-"<<_ey[i]<<" .. "
	       <<_z[i]<<"+-"<<_ez[i]<<std::endl;
    }
  }
  //__________________________________________________________________
  int MuonTrack::getNhit(){return _z.size();}
  //__________________________________________________________________
  void MuonTrack::trackFinder(){
    // chi2distribution == [0]*1/(TMath::Gamma([1]/2.)*2^([1]/2))*x^(([1]/2)-1)*exp(-x/2.) 
    
    std::cout<< "track finder"<<std::endl;
  }
  double MuonTrack::getProb(){
    return TMath::Prob(_Parameter[4],(_z.size())-4);
  }
} /// namespace TRACK
