#include "src/RecoObject.cpp"

#ifndef MUON_HIT
#define MUON_HIT


namespace Muon {

  class Hit : public RecoObject {
  public:
    Hit();
    Hit(double tdc_time, double adc_time, double drift_time, 
           double corr_time, int _tdc, int _chan);
    ~Hit() {};

    double     TDCTime  () const;
    double     ADCTime  () const;
    double     DriftTime() const;
    double     CorrTime () const;
    int        TDC      () const;
    int        Channel  () const;
    void       Draw     () override;     
    static double     RadiusError(double radius);     
  private:
    double     tdctime;
    double     adctime;
    double     driftTime;
    double     corrTime;
    int        tdc;
    int        chan;

  };

  Hit::Hit() {
    tdctime    = 0;
    adctime    = 0;
    driftTime  = 0;
    corrTime   = 0;
  }
  
  Hit::Hit(double tdc_time, double adc_time, double drift_time, 
           double corr_time, int _tdc, int _chan) {
    tdctime    = tdc_time;
    adctime    = adc_time;
    driftTime  = drift_time;
    corrTime   = corr_time;
    tdc        = _tdc;
    chan       = _chan;
  }

  double Hit::TDCTime()   const{return tdctime;}
  double Hit::ADCTime()   const{return adctime;}
  double Hit::DriftTime() const{return driftTime;}
  double Hit::CorrTime()  const{return corrTime;}
  int    Hit::TDC()       const{return tdc;}
  int    Hit::Channel()   const{return chan;}  
  static double RadiusError(double radius);  

  void   Hit::Draw()           {}
  
  // double Hit::TimeError() {
  //   if (corrTime < 90)
  //     return 0.11;
  //   else
  //     return 0.11-corrTime/3.0;
  // }

  double Hit::RadiusError(double radius) {

  //   /*
  //   if (radius <3.5) 
  //     return .25-.028*radius;//(.2-.0400*radius+.00475*radius*radius);
  //   else
  //     return .15;
  //   */

  //   //return .2-.034*radius + .00384*radius*radius;
  //   //return 1.3*(.172-.0377*radius+.003815*radius*radius);
  //   // return .217-0.053*radius+0.0053*radius*radius
  //   //return .237-.068*radius+.007*radius*radius;
  //   //return .256-.081*radius+0.00856*radius*radius;
    
  //   // nominal resolution:
     return .226-0.025*radius+0.001248*radius*radius;
  }

  // void Hit::SetRadius(double r) {
  //   radius = r;
  // }

  
}

#endif
