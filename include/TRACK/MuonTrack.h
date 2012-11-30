// -*- C++ -*-
#ifndef TRACK_MUONTRACK_H
#define TRACK_MUONTRACK_H 1
#include <vector>
#include <TROOT.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TMinuit.h>
#include <map>


namespace TRACK {
  
  /** Implementation of MCParticle.
   * 
   * @author Yacine Haddad
   * @version $Id: v0.0 $
   */
  
  class MuonTrack{
    
  public: 
    MuonTrack();
    virtual ~MuonTrack(){};
    
    /** Set The Hits position of muons track.
     */
    //virtual void setHits(EVENT::CalorimeterHit *hit);
    virtual void setX(double x,double ex);
    virtual void setY(double y,double ey);
    virtual void setZ(double z,double ez);
    /** Set The Hits position of muons track.
     */
    virtual int getNhit();
    
    /** Reset all intern variable 
     */
    virtual void reset();
    
    /** 
	@return chi2 value.
     */
    virtual double getChi2();
    
    /** 
	@return chi2 value.
     */
    virtual double *getParameter();
    
    /** Distance of one hit in .
     */
    static double distance2(double z,
			    double x,
			    double y,
			    double ex,
			    double ey,
			    double *p);
    
    /** Distance between hit_reco and cluster
     */
    //static double distance(double z,
    //			   double x,
    //			   double y,
    //			   double *p);
    //
    /** Minimisation function.
     */
    static   void     linearMinFunction(int &, 
					double *,
					double & sum, 
					double * par, 
					int);
    
    /** fitting function.
     */
    virtual  void     trackFitter(double *p);
    
    /** get the reconstructed track function.
     */
    virtual  void     trackReco();
    //virtual  void     getCallibrationConst();
    virtual  void     setNlayer(int Nlayer){_Nlayer = Nlayer;};
    
    /** print the the position.
     */
    virtual  void     positionStatus();
    
    /** print the the position.
     */
    //    virtual void      getSigma();  
    

    /** use the hough transform to find the track
     */
    virtual void     trackFinder();


    virtual double   getProb();
    
  protected:
    // the poisitions of the hit
    std::vector<double> _x;
    std::vector<double> _y;
    std::vector<double> _z;

    std::vector<double> _ex;
    std::vector<double> _ey;
    std::vector<double> _ez;
    
    std::map<int,double*> found_track;
    
    int _Nlayer;
    //EVENT::CalorimeterHitVec _Event;
    /** Fit parameter 
     * p[0..3] = the line parameter 
     * p[4] = the chi2
     */
    
    double _Parameter[5];
    
  }; // class
} // namespace IMPL
#endif /* ifndef IMPL_MCPARTICLEIMPL_H */
