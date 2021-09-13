//////////////////////////////////////
////header for varibles and functions
////
//////////////////////////////////////


#ifndef variable_function_h
#define variable_function_h



#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
//#include <assert.h>
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <sys/mman.h>
//#include <fcntl.h>
#include <unistd.h>
#include <TTree.h>
#include <iomanip>
#include <stdio.h>
#include <TFile.h>
#include "TF1.h"
#include <TLeaf.h>


int m_inoise ; //npoints/4

TBranch *b_w4f;
TBranch *b_w0f;
TBranch *b_w3f; ///!
TBranch *b_w2f; ///!
TBranch *b_w1f;
TBranch *b_w5f;


TBranch *b_t4f;
TBranch *b_t3f; ///!
TBranch *b_t2f; ///!
TBranch *b_t1f;
TBranch *b_boardf;

//double  w11[2000]; //pulses holder
//std::vector<double> t0, t1, t2, t3, t4, t5, w0, w1, w2, w3, w4, w5;
std::vector<int> board0;
std::vector<double> b0, b1, b2, b3, b4, b5;
std::vector<double> x1, x2, y11, y2;
std::pair<double,double> noise2, noise3;
std::pair<double,double> baseline2, baseline3;

//new Tfile and tree generation.

//////////////////////


//////variables
std::vector<double>  pmax0, tmax0, rise0_1040, rise0_1090, pulse_area0, rms0, dvdt0,  t01;
std::vector<double>  pmax1, tmax1, rise1_1040, rise1_1090, pulse_area1, rms1, dvdt1,  t02;
std::vector<double>  pmax2, tmax2, rise2_1040, rise2_1090, pulse_area2, rms2, dvdt2,  t03;
std::vector<double>  pmax3, tmax3, rise3_1040, rise3_1090, pulse_area3, rms3, dvdt3,  t04;
std::vector<double>  pmax4, tmax4, rise4_1040, rise4_1090, pulse_area4, rms4, dvdt4,  t05;
std::vector<double>  pmax5, tmax5, rise5_1040, rise5_1090, pulse_area5, rms5, dvdt5,  t06;
std::vector<double>  w1, w2, w3, w4, t1, t2, t3,t4, cfd1, cfd2,cfd3, cfd4, wbase0, wbase1, wbase2, wbase3, wbase4, wbase5;
/////////////////

void base_line(int npoints, std::vector<double>& w, int m_inoise)
{
	double mean =0;

	for(int j = 0; j < m_inoise; j++){mean += w[j];}
	
	mean = mean/m_inoise;

	for(int j = 0; j < npoints; j++){w[j] = w[j]- mean;}
} 

bool pulse_baseline(std::vector<double> w, std::pair<double, double>& baseline, std::pair<double, double>& noise_rms)

{

  baseline = std::make_pair(-1, -1);

  noise_rms = std::make_pair(-1, -1);

  double mean = 0.0;

  double strd = 0.0;

  for (int j = 0; j < m_inoise; j++) mean += w[j];

  mean = mean / m_inoise;

  for (int i = 0; i < m_inoise; ++i) strd += pow(w[i] - mean, 2);

  strd = sqrt(strd / m_inoise);

  double limit1 = ceil((mean - 8*strd)*1000)/1000;

  double limit2 = ceil((mean + 8*strd)*1000)/1000;

  int bins = ceil(16*strd*1000);

  TH1D* noise = new TH1D("PulseNoise", "PulseNoise", bins, limit1, limit2);

  for (int i = 0; i<m_inoise; i++) noise->Fill(w[i]);

  double m = 0.0;

  double dm = 0.0;

  double sigma = 0.0;

  double dsigma = 0.0;

  if (noise->Integral()>0)

     {

      double mx = noise->GetMean();

      double rms = noise->GetRMS();

      double rmin = mx - 3 * rms;

      double rmax = mx + 3 * rms;

      TF1* mygau = new TF1("mygau", "gaus", rmin, rmax);

      noise->Fit(mygau,"Q0");

      m = mygau->GetParameter(1);

      dm = mygau->GetParError(1);

      sigma = sqrt(mygau->GetParameter(2));

      dsigma = 0.5*pow(mygau->GetParameter(2), -0.5)*mygau->GetParError(2);

      mygau->Delete();

     }

  noise->Delete();

  baseline = std::make_pair(m, dm);

  noise_rms = std::make_pair(sigma, dsigma);

  return true;

}

double pulse_max( int npoints, std::vector<double> w, int inoise){ // number of data points in event, voltage vector,
    // function to calculate pulse maximum
    double pmax =0.0;
    for( int j = 0; j < npoints; j++){
        if (j == 0){pmax = 0.0;}
        if(j != 0 && w[j] > pmax){pmax = w[j];}// update pulse max
    }
    return pmax;
}

double pulse_max_fit( int npoints, std::vector<double> w,std::vector<double> t , double bottom, double top, int inoise){ // number of data points in event, voltage vector,
    // function to calculate pulse maximum
    double rise{0}, pmax =0, tmax{0}, lowerval = 0, upperval = 0, m1=0, m2=0, tbottom =0, ttop=0, peak{0};
    int imax=0;
    bool ten=true, ninety=true;

    for( int j = inoise; j < npoints; j++)
    {
        if(w[j] > pmax)
        {
            imax = j;
            pmax = w[j];
            tmax = t[j];
            lowerval = w[j]*bottom;
            upperval = w[j]*top;
        }
    } // find index of max, pulse max, and the amplitudes at 10%,90% of pmax
    for( int j = imax; j > -1; j--)
    {
        if( w[j]<lowerval)
        {
            tbottom=t[j];      //find the index right below 10%
            break;
        }
    }
    for( int j = imax; j < npoints; j++)
      {
          if(w[j]<upperval)
          {
              ttop=t[j];     //find the index right below 90% (past peak value)
              break;
          }
      }
    //fit with gaussian
    TF1* f = new TF1("f","gaus",tbottom,ttop,"r");
    f->SetParLimits(0,pmax*0.5,pmax*1.05);
    f->SetParameter(0,10);
    f->SetParLimits(1,tbottom,ttop);
    f->SetParameter(1,tmax);
    f->SetParLimits(2,0,100);
    f->SetParameter(2,15);
    TGraphErrors* g = new TGraphErrors();
    for( int j = 0; j < npoints; j++){g->SetPoint(j,t[j],w[j]);}
    g->Fit("f","","0",tbottom,ttop);
  
    peak = f->GetParameter(0);
    delete f;
    delete g;
    return peak;
}




double rise_fit( int npoints, std::vector<double> w,std::vector<double> t , double bottom, double top, int inoise){ // number of data points in event, voltage vector,
    // function to calculate pulse maximum
    double rise{0}, pmax =0, tmax{0}, lowerval = 0, upperval = 0, m1=0, m2=0, tbottom =0, ttop=0, sigma{0};
    int imax=0;
    bool ten=true, ninety=true;

    for( int j = inoise; j < npoints; j++)
    {
        if(w[j] > pmax)
        {
            imax = j;
            pmax = w[j];
            tmax = t[j];
            lowerval = w[j]*bottom;
            upperval = w[j]*top;
        }
    } // find index of max, pulse max, and the amplitudes at 10%,90% of pmax
    for( int j = imax; j > -1; j--)
    {
        if( w[j]<lowerval)
        {
            tbottom=t[j];      //find the index right below 10%
            break;
        }
    }
    for( int j = imax; j < npoints; j++)
      {
          if(w[j]<upperval)
          {
              ttop=t[j];     //find the index right below 90% (past peak value)
              break;
          }
      }
    //fit with gaussian
    TF1* f = new TF1("f","gaus",tbottom,ttop,"r");
    f->SetParLimits(0,0,100);
    f->SetParameter(0,10);
    f->SetParLimits(1,tbottom,ttop);
    f->SetParameter(1,tmax);
    f->SetParLimits(2,0,100);
    f->SetParameter(2,15);
    TGraphErrors* g = new TGraphErrors();
    for( int j = 0; j < npoints; j++){g->SetPoint(j,t[j],w[j]);}
    g->Fit("f","","0",tbottom,ttop);
  
    rise = f->GetParameter(2);
    delete f;
    delete g;
    rise = rise *1.7;
    return rise;
}


double time_max( int npoints, std::vector<double> t, std::vector<double> w, int inoise){   //number of data points in event, voltage vector, time vector
    // function to calculate time at pulse maximum
    double pmax = 0.0, tmax = 0.0;
    for (int j = 0; j < npoints; j++){
        
        if (j == 0){
            pmax = 0.0;
            tmax = 0.0;
        }
        
        if( j != 0 && w[j] > pmax){
            pmax = w[j];
            tmax = t[j];
        }
    }
    return tmax;
}






double pulse_area( int npoints, std::vector<double> t, std::vector<double> w, int inoise){// number of data points in event, time vector, voltage vector, termination resisitace in ohms, effective gain from amplifier(s)

    double pulse_area = 0.0, time_difference, pmax =0; int imax =  0, istart = 0, iend = 0; // collected charge, size of time bins; pulse max, index of pulse max, start and end. Start and and end defined as index of first points before and after max that are less than or equal to zero.
    //NEED to add better baseline correction!///////////////
    time_difference = t[1] - t[0];
    for( int j = 0; j < npoints; j++){      if(w[j] > pmax){imax = j; pmax = w[j];}    } // find index of pulse maximum
    for( int j = imax; j>-1; j--){     if(w[j] <= 0){ istart = j+1; break;}  } // find index of start of pulse
    for( int j = imax; j< npoints; j++){     if(w[j] <= 0){ iend = j-1; break;}  } // find index of end of pulse
    if (istart ==0){ return 0;}
    
    for( int j = istart; j < iend + 1; j++){
        pulse_area = pulse_area + w[j]/1000 * (t[j+1] - t[j-1])/2 ;
    }
    pulse_area = pulse_area /1E12;
    return pulse_area; // collected charge in Coulombs, assuming termination is in ohms and voltage is in volts, time is in seconds
}


double rise_int( double x1, double y1, double x2, double y2, double y3){
    double m, x3;
    
    m = (y2 - y1) / (x2 - x1);
    
    x3 = ( y3 - y1 ) / m + x1;
    
    return x3;
}

double cfd_index(int npoints, int fraction, std::vector<double> t, std::vector<double> w){
    // function to calculate index of constant fraction - not truly a constant fraction discriminator
    
    double pmax = 0, time_fraction = 0;
    int imax = 0, ifraction =0;
    
    // find index of pulse maxima
    for( int j = 0; j < npoints; j++){
        if(w[j] > pmax){
            imax = j;
            pmax = w[j];
        }
    }
    
    for( int j = imax; j>-1; j--){
        if(w[j] <= pmax*double(fraction)/100){
            ifraction = j;              //find index of first point before constant fraction of pulse
            time_fraction = t[j];
            break;
        }
        
    }
    time_fraction = time_fraction + (t[ifraction+1] - t[ifraction])* (pmax*double(fraction)/100 - w[ifraction]) /(w[ifraction+1] - w[ifraction]);
    return time_fraction;
}

double cfd_index_discrete(int npoints, int fraction, std::vector<double> t, std::vector<double> w){
    // function to calculate index of constant fraction - not truly a constant fraction discriminator
    
    double pmax = 0, time_fraction = 0;
    int imax = 0, ifraction =0;
    
    // find index of pulse maxima
    for( int j = 0; j < npoints; j++){
        if(w[j] > pmax){
            imax = j;
            pmax = w[j];
        }
    }
    
    for( int j = imax; j>-1; j--){
        if(w[j] <= pmax*double(fraction)/100){
            ifraction = j;              //find index of first point before constant fraction of pulse
            time_fraction = t[j];
            break;
        }
        
    }
    time_fraction = time_fraction + (t[ifraction+1] - t[ifraction])* (pmax*double(fraction)/100 - w[ifraction]) /(w[ifraction+1] - w[ifraction]);
    
    if (time_fraction  - t[ifraction] <= (t[ifraction+1] - t[ifraction]/2.0) ){ time_fraction = t[ifraction];
    }else {time_fraction =  t[ifraction + 1];}
    
    return time_fraction;
}



double pulse_min( int npoints, std::vector<double> w){  // number of data points in event, voltage vector
    // function to calculate pulse minimum
    double pmin = 0;
    for (int j = 0; j < npoints; j++){        if (j == 0){ pmin = 0;}
        if (w[j] < pmin){ pmin = w[j];}
    }
    return pmin;
}



double rise_time( int npoints, std::vector<double> t, std::vector<double> w, double bottom, double top){ // number of data points in event, time vector, voltage vector
    // function to calculate 10-90 rise time
    // would be nice to add linear (or better) interpolation
    double rise{0}, pmax =0, lowerval = 0, upperval = 0, m1=0, m2=0, tbottom =0, ttop=0;
    int imax=0, itop = 0, ibottom = 0;
    bool ten=true, ninety=true;

    
    
    for( int j = 0; j < npoints; j++)
    {
        if(w[j] > pmax)
        {
            imax = j;
            pmax = w[j];
            lowerval = w[j]*bottom;
            upperval = w[j]*top;
        }
    } // find index of max, pulse max, and the amplitudes at 10%,90% of pmax
    
    for( int j = imax; j > -1; j--)
    {
        if(ninety && w[j]<upperval)
        {
            itop=j;     //find the index right below 90%
            ninety = false;
        }
        
        if(ten && w[j]<lowerval)
        {
            ibottom=j;      //find the index right below 10%
            ten = false;
        }
        
        if( !ten && !ninety ){
            break;
        }
    }
    
    tbottom = rise_int( t[ibottom], w[ibottom], t[ibottom+1], w[ibottom + 1], lowerval);
    ttop = rise_int( t[itop], w[itop], t[itop + 1], w[itop + 1], upperval);
    rise = ttop - tbottom; // rise
    return rise;
}




double th_t0( int npoints, std::vector<double> t, std::vector<double> w, double thresh){ // number of data points in event, time vector, voltage vector
    // function to calculate time of threshold crossing
    // would be nice to add linear (or better) interpolation
    double t0, pmax =0;
    int imax=0, itop = 0;
    
    
    for( int j = 0; j < npoints; j++)
    {
        if(w[j] > pmax)
        {
            imax = j;
            pmax = w[j];
        }
    } // find index of max, pulse max, and the amplitudes at 10%,90% of pmax
    
    for( int j = imax; j > -1; j--)
    {
        if(w[j]<thresh)
        {
            itop=j;     //find the index of point just below threhsold
            break;
        }
    }
    
    t0 = t[itop]; // rise
    return t0;
}




double noise_rms( int npoints, std::vector<double> w, int inoise){// number of data points in event, time vector, voltage vector, index at end of noise
    // function to calculate pulse maximum
    double rms = 0, mean =0, var;
    for(int j = 0; j < inoise; j++){rms += w[j]*w[j]; mean += w[j];}
    mean = mean / inoise;
    rms = rms/inoise;
    var = rms - mean * mean;
    rms = pow(var, 0.5); // stupid convention where we call the standard deviation the rms
    return rms;
}





double pulse_dvdt_cfd(int npoints, int fraction,   int ndif, std::vector<double> t, std::vector<double> w){
    // function to dv/dt at a given constant fraction value.
    double pmax = 0,  time_difference, dvdt;
    int imax = 0,ifraction =0;
    //time_difference = t[1] -t[0];
    // find index of pulse maxima
    for( int j = 0; j < npoints; j++){  if(w[j] > pmax){imax = j; pmax = w[j]; }    }
    for( int j = imax; j>-1; j--){     if(w[j] <= pmax*double(fraction)/100){ ifraction = j; break;}  }//find index of first point before constant fraction of pulse
    if(ndif ==0){         dvdt = (w[ifraction+1] - w[ifraction])/(t[ifraction + 1] - t[ifraction]);}
    else {dvdt = (w[ifraction+ndif] - w[ifraction-ndif])/(t[ifraction+ndif] - t[ifraction-ndif]);}
    return dvdt;
}










#endif // #ifdef anal_cxx
