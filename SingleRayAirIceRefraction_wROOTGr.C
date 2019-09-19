#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>

const double pi=4.0*atan(1.0); /**< Gives back value of Pi */
const double spedc=299792458.0; /**< Speed of Light in m/s */
const double A_ice=1.78;/**< Value of Parameter A for SP model of refractive index */
const double B_ice=-0.43; /**< Value of Parameter B for SP model of refractive index */
const double C_ice=0.0132;/**< Value of Parameter A for SP model of refractive index. There will be a negative sign for postive depth depth. */

////Define vectors to store data from the file
vector <vector <double>> nh_data;////n(h) refractive index profile of the atmosphere as a function of height
vector <vector <double>> lognh_data;////log(n(h)-1) log of the refractive index profile of the atmosphere as a function of height subtracted by 1
vector <vector <double>> h_data;////height data

////Define Arrays for storing values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file
double ATMLAY[5];
double abc[5][3];

////This Function reads in the values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file. The a,b and c values are mass overburden values and are not required in this code.
int readATMpar(){

  ////Open the file
  ifstream ain("Atmosphere.dat");
  
  int n1=0;////variable for counting total number of data points
  string line;
  double dummya[5]={0,0,0,0,0};////temporary variable for storing data values from the file
  
  //Check if file is open and store data
  if(ain.is_open()){
    
    while (getline(ain,line)){

      if(n1<4){////only read in the lines which contain the ATMLAY and a,b and c values in the file
	ain>>dummya[0]>>dummya[1]>>dummya[2]>>dummya[3]>>dummya[4];
	//cout<<n1<<" "<<dummya[0]<<" , "<<dummya[1]<<" , "<<dummya[2]<<" , "<<dummya[3]<<" , "<<dummya[4]<<endl;
      }

      ////Store the values in their respective arrays
      if(n1==0){
	for (int i=0; i<5; i++){ ATMLAY[i]=dummya[i]; }
      }    
      if(n1==1){
	for (int i=0; i<5; i++){ abc[i][0]=dummya[i]; }
      }
      if(n1==2){
	for (int i=0; i<5; i++){ abc[i][1]=dummya[i]; }
      }
      if(n1==3){
	for (int i=0; i<5; i++){ abc[i][2]=dummya[i]; }
      }
      n1++;
    }////end the while loop
    
    ain.close();
  }////if condition to check if file is open

  return 0;
}

int readnhFromFile(){

  ////Open the file
  ifstream ain("Atmosphere.dat");
  ain.precision(10); 

  int n1=0;////variable for counting total number of data points
  int layer=0;
  string line;

  ////Ignore the lines containing ATMLAY and a,b and c values.
  for(int i=0; i<5; i++){ ain.ignore(256,'\n'); }
  
  ////Check if file is open and store data
  if(ain.is_open()){
    ////define dummy/temporary variables for storing data
    double dummy1,dummy2;
    ////define dummy/temporary vectors for storing data.
    vector <double> temp1,temp2,temp3;
    
    while (getline(ain,line)){
      ain>>dummy1>>dummy2;
      
      if(dummy1>=0){////start storing height at above and equal to 0 m
	////push in the height values for a single layer in the temporary vector
	temp1.push_back(dummy1);
	temp2.push_back(dummy2);
	temp3.push_back(log(dummy2-1));
	
	if(dummy1*100>=ATMLAY[layer]){////change the layer once the data of all the heights of that layer has been read in
	  if(layer>0){////now since the layer has finished and the temporary vectors have been filled in. Now we push the vectors in the main 2d height and refractice index vectors
	    h_data.push_back(temp1);
	    nh_data.push_back(temp2);
	    lognh_data.push_back(temp3);
	    ////clear the vectors now for storing the next layer
	    temp1.clear();
	    temp2.clear();
	    temp3.clear();
	  } 
	  layer++;
	}
	n1++;
      } 
    }////end the while loop
    
    if(layer>0){////For storing the last layer
      h_data.push_back(temp1);
      nh_data.push_back(temp2);
      lognh_data.push_back(temp3);
      ////clear the vectors now for storing the next layer
      temp1.clear();
      temp2.clear();
      temp3.clear();
    }
    layer++;
    
    ain.close();
  }////if condition to check if file is open

  ////The file reading condition "while (getline(ain,line))" reads the last the datapoint of the file twice. This is to to remove the last repeat data point in all the data arrays
  h_data[h_data.size()-1].erase(h_data[h_data.size()-1].end() - 1);
  nh_data[nh_data.size()-1].erase(nh_data[nh_data.size()-1].end() - 1);
  lognh_data[lognh_data.size()-1].erase(lognh_data[lognh_data.size()-1].end() - 1);
  
  return 0;
}

/////Functions used for Raytracing in Ice using the analytical solution

////Analytical solution describing the ray path in ice
struct fDnfR_params { double a, b, c, l; };
double fD(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
struct fdxdz_params { double a, b, c, lang, z0,z1; };
double fdxdz(double x,void *params){
  
  struct fdxdz_params *p= (struct fdxdz_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Lang = p->lang;
  double Z0 = p->z0;
  double Z1 = p->z1;
  
  return tan(asin( ((A+B*exp(C*Z0))*sin(x))/(A+B*exp(C*Z1)) ) ) - tan(Lang);
}

////The function used to calculate ray propogation time in ice
struct ftimeD_params { double a, b, c, speedc,l; };
double ftimeD(double x,void *params){

  struct ftimeD_params *p= (struct ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;
  
  return (1.0/(Speedc*C*sqrt(pow(A+B*exp(x*C),2)-L*L)))*(pow(A+B*exp(x*C),2)-L*L+(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)))*(A*A*sqrt(pow(A+B*exp(x*C),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(A+B*exp(x*C),2)-L*L)*log(A+B*exp(x*C)+sqrt(pow(A+B*exp(x*C),2)-L*L)) );
}

////Use GSL minimiser which uses Brent's Method to find root for a given function
double FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  //printf ("using %s method\n", gsl_root_fsolver_name (s));
  //printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);

      if (status == GSL_SUCCESS){
	//printf ("Converged:");
	//printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

////Get the distance on the 2ndLayer at which the ray should hit given an incident angle such that it hits an target at depth of z0 m in the second layer.
//// n_layer1 is the refractive index value of the previous layer at the boundary of the two mediums
//// A,B and C values are the values required for n(z)=A+B*exp(Cz) for the second layer
//// TxDepth is the starting height or depth
//// RxDepth is the final height or depth
//// AirOrIce variable is used to determine whether we are working in air or ice as that sets the range for the GSL root finder.
double *GetLayerHitPointPar(double n_layer1,double A, double B, double C, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){
  double *output=new double[3];
  
  double x0=0;////Starting horizontal point of the ray. Always set at zero
  double x1=0;////Variable to store the horizontal distance that will be traveled by the ray
  
  double ReceiveAngle=0;////Angle from the vertical at which the target will recieve the ray
  double Lvalue=0;//// L parameter of the ray for that layer
  double RayTimeIn2ndLayer=0;////Time of propagation in 2ndLayer 
  double AngleOfEntryIn2ndLayer=0;////Angle at which the ray enters the layer

  double SurfaceRayIncidentAngle=IncidentAng*(pi/180.0);////Angle at which the ray is incident on the second layer
  double RayAngleInside2ndLayer=asin((n_layer1/(A+B*exp(-C*TxDepth)))*sin(SurfaceRayIncidentAngle));////Use Snell's Law to find the angle of transmission in the 2ndlayer

  double GSLFnLimit=90*(pi/180);
  if(AirOrIce==0){
    //cout<<"we are in ice"<<endl;
    GSLFnLimit=50*(pi/180);
  }
  if(AirOrIce==1){
    //cout<<"we are in air"<<endl;
    GSLFnLimit=88.5*(pi/180);
  }
  
  ////calculate the angle at which the target receives the ray
  gsl_function F1;
  struct fdxdz_params params1 = {A, B, -C, RayAngleInside2ndLayer, RxDepth,TxDepth};
  F1.function = &fdxdz;
  F1.params = &params1;
  ReceiveAngle=FindFunctionRoot(F1,0*(pi/180),GSLFnLimit);
  //cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/pi)<<" deg"<<endl;
  
  ////calculate the distance of the point of incidence on the 2ndLayer surface and also the value of the L parameter of the solution
  Lvalue=(A+B*exp(-C*RxDepth))*sin(ReceiveAngle);
  struct fDnfR_params params2 = {A, B, -C, Lvalue};
  x1=+fD(RxDepth,&params2)-fD(TxDepth,&params2);
  if(AirOrIce==1){
    x1*=-1;
  }
  //cout<<"The hit point horizontal distance is from the Rx target "<<x1<<" m  on the surface"<<endl; 

  ////calculate the propagation time in 2ndLayer
  struct ftimeD_params params3 = {A, B, -C, spedc,Lvalue};
  RayTimeIn2ndLayer=+ftimeD(RxDepth,&params3)-ftimeD(TxDepth,&params3);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  //cout<<"The propagation time in 2ndLayer is: "<<RayTimeIn2ndLayer<<" s"<<endl;

  ///////calculate the initial angle when the ray enters the 2ndLayer. This should be the same as RayAngleInside2ndLayer. This provides a good sanity check to make sure things have worked out.
  // gsl_function F4;
  // double result, abserr;
  // F4.function = &fD;
  // F4.params = &params2;
  // gsl_deriv_central (&F4, TxDepth, 1e-8, &result, &abserr);
  // AngleOfEntryIn2ndLayer=atan(result)*(180.0/pi);
  // if(TxDepth==RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=180-ReceiveAngle;
  // }
  // if(TxDepth!=RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=90;
  // }
  // cout<<" ,AngleOfEntryIn2ndLayer= "<<AngleOfEntryIn2ndLayer<<" ,RayAngleInside2ndLayer="<<RayAngleInside2ndLayer*(180/pi)<<endl;

  output[0]=x1;
  output[1]=ReceiveAngle*(180/pi);
  output[2]=Lvalue;
  
  return output;
}

////This function flattens out 2d vectors into 1d vectors
vector<double> flatten(const vector<vector<double>>& v) {
    size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    vector<double> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

////This function is used to measure the amount of time the code takes
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp (){
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void SingleRayAirIceRefraction_wROOTGr(){

  ////For recording how much time the process took
  timestamp_t t0 = get_timestamp();
  
  double AntennaDepth=200;////Depth of antenna in the ice
  double RayLaunchAngle=170;////Initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 is vertically down
  double TxHeight=20000;////Height of the source
  double IceLayerHeight=3000;////Height where the ice layer starts off

  if(TxHeight>h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]){
    ////Maximum height available with the refractive index data
    cout<<"Tx Height is set higher than maximum available height for atmospheric refractive index which is "<<h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]<<endl;
    TxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];
    cout<<"Setting Tx Height to be the maximum available height"<<endl;
  }

  if(RayLaunchAngle<=90){
    cout<<"RayLaunchAngle has been set at "<<RayLaunchAngle<<" deg which is outside of the allowed range of 90 deg <RayLaunchAngle< 180 deg"<<endl;
    RayLaunchAngle=135;
    cout<<"Setting RayLaunchAngle at"<<RayLaunchAngle<<endl;
  }
  
  ////Print out the ray path x and y values in a file
  ofstream aout("RayPathinAirnIce.txt");
  
  ////Fill in the n(h) and h arrays and ATMLAY and a,b and c (these 3 are the mass overburden parameters) from the data file
  readATMpar();
  readnhFromFile();
  int MaxLayers=h_data.size();////store the total number of layers present in the data
  
  ////Flatten out the height and the refractive index vectors to be used for setting the up the spline interpolation.
  vector <double> flattened_h_data=flatten(h_data);
  vector <double> flattened_nh_data=flatten(nh_data);

  ////Set up the GSL cubic spline interpolation. This is used for interpolating values of refractive index at different heights.
  gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,flattened_h_data.size());
  gsl_spline_init(spline, flattened_h_data.data(), flattened_nh_data.data(), flattened_h_data.size());

  ////Plot and Fit each of the three layers to the function n(h)-1=ln(B)+C*h
  TCanvas *cc=new TCanvas("cc","cc");
  cc->cd();
  TGraph * gr=new TGraph();
  int isample=0;
  for(int i=0;i<h_data.size();i++){
    for(int j=0;j<h_data[i].size();j++){
      gr->SetPoint(isample,h_data[i][j],log(nh_data[i][j]-1));
      isample++;
    }
  }
  gr->SetTitle("Refractive index profile of the atmosphere; Height (m); ln(n(h)-1)");
  gr->Draw("ALP");

  TF1 *nz=new TF1("nz","[0]+[1]*x",0,50000);
  nz->SetParameters(-7.1,-0.0001);
  ////Now fit the layers witht the predefined function
  for(int i=0;i<h_data.size();i++){
    ////Start from Lowest Layer in altitude-Layer 1
    //gr->Fit(nz,"","",ATMLAY[i]/100,ATMLAY[i+1]/100);
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Section for propogating the ray through the atmosphere
  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
    //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<endl;
    if(TxHeight<ATMLAY[ilayer]/100 && TxHeight>ATMLAY[ilayer-1]/100){
      cout<<"Tx Height is in this layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer-1]/100<<" m and is at a height of "<<TxHeight<<" m"<<endl;
      ilayer=-100;
    }
    if(ilayer>-1){
      skiplayer++;
    }
  }
  int SkipLayersAbove=skiplayer;
  //cout<<"The tota number of layers that need to be skipped from above is "<<skiplayer<<endl;

  ////Find out how many atmosphere layers are below the ice height which we do not need
  skiplayer=0;
  for(int ilayer=0;ilayer<MaxLayers;ilayer++){
    //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer+1]/100<<endl;
    if(IceLayerHeight>ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
      cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
      ilayer=100;
    }
    if(ilayer<MaxLayers){
      skiplayer++;
    }
  }
  int SkipLayersBelow=skiplayer;
  //cout<<"The total number of layers that need to be skipped from below is "<<skiplayer<<endl;
  
  ////Define variables for ray propogation through mutliple layers in the atmosphere
  double Start_nh=0;
  double StartHeight=0;
  double StopHeight=0;
  double StartAngle=0;
  double A=0,B=0,C=0;
  double TotalHorizontalDistance=0;
  vector <double> layerAs,layerBs,layerCs,layerLs;////vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them

  double c0, c1;////variables to store the fit results from the GSL linear fitter for y=c0+c1*x. Here we are trying to fit the function: log(n(h)-1)=log(B)+C*h which basically comes from n(h)=A+B*exp(C*h)
  double cov00, cov01, cov11, chisq;////variables to store the covariance matrix elements outputted the GSL linear fitting function

  ////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
  cout<<"Fitting the atmosphere refrative index profile with multiple layers and propogate the ray"<<endl;
  for(int ilayer=MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){

    ////Run the GSL fitter to get the refractive index profile for the layer
    gsl_fit_linear (h_data[ilayer].data(), 1, lognh_data[ilayer].data() ,1, h_data[ilayer].size(),
  		    &c0, &c1, &cov00, &cov01, &cov11,
  		    &chisq);
    ////Store the results from the GSL fitter
    ////A, B and C parameters for the refractive index profile
    A=1.0;////This value is fixed at 1 as that is what the refractive index is for vacuum
    B=exp(c0);
    C=-c1;
    
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=TxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=ATMLAY[ilayer+1]/100;
    }

    ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
    Start_nh=gsl_spline_eval(spline, StartHeight, accelerator);

    ////Set the stopping height of the ray for propogation for that layer
    if(ilayer==(SkipLayersBelow-1)+1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      StopHeight=IceLayerHeight;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      StopHeight=ATMLAY[ilayer]/100;
    }

    ////If this is the first layer then set the initial launch angle of the ray through the layers
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      StartAngle=180-RayLaunchAngle;
    }
    //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;

    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    double* GetHitPar=GetLayerHitPointPar(Start_nh, A, B, C, StopHeight, StartHeight, StartAngle, 1);
    TotalHorizontalDistance+=GetHitPar[0];
    StartAngle=GetHitPar[1];

    ////Store in the values of A,B,C and L for tha layer
    layerLs.push_back(GetHitPar[2]);
    layerAs.push_back(A);
    layerBs.push_back(B);
    layerCs.push_back(C);
    
    // printf ("# best fit: Y = %g + %g X\n", c0, c1);
    // // printf ("# covariance matrix:\n");
    // // printf ("# [ %g, %g\n#   %g, %g]\n",cov00, cov01, cov01, cov11);
    // printf ("# chisq = %g\n", chisq);

    ////dont forget to delete the pointer!
    delete []GetHitPar;
  }

  double IncidentAngleonIce=StartAngle;
  cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<endl;
  cout<<"Now treating the atmosphere refrative index profile as a single layer and fitting it and propogating the ray"<<endl;

  ////Run the GSL fitter to get the refractive index profile for the whole atmosphere
  gsl_fit_linear (flatten(h_data).data(), 1, flatten(lognh_data).data() ,1, flatten(h_data).size() ,
                   &c0, &c1, &cov00, &cov01, &cov11,
                   &chisq);
  ////A, B and C parameters for the refractive index profile
  A=1.0;////This value is fixed at 1 as that is what the refractive index is for vacuum
  B=exp(c0);
  C=-c1;
  
  ////Set the starting height of the ray for propogation to be the height of the transmitter
  StartHeight=TxHeight;
  ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
  Start_nh=gsl_spline_eval(spline, StartHeight, accelerator);
  ////Set the stopping height of the ray for propogation to be the height of the ice layer
  StopHeight=IceLayerHeight;
  ////Set the initial launch angle of the ray
  StartAngle=180-RayLaunchAngle;
  //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;
  
  ////Get the hit parameters from the function. The output is:
  //// How much horizontal distance did the ray travel through the whole atmosphere
  //// The angle of reciept/incidence at the end 
  //// The value of the L parameter for whole atmosphere fit
  double* GetHitPar=GetLayerHitPointPar(Start_nh, A, B, C, StopHeight, StartHeight, StartAngle, 1);

  ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
  double TotalHorizontalDistanceSLF=GetHitPar[0];
  double StartAngleSLF=GetHitPar[1];
  double LvalueSLF=GetHitPar[2];

  //printf ("# best fit: Y = %g + %g X\n", c0, c1);
  ////printf ("# covariance matrix:\n");
  ////printf ("# [ %g, %g\n#   %g, %g]\n",cov00, cov01, cov01, cov11);
  //printf ("# chisq = %g\n", chisq);
  
  cout<<"Total horizontal distance travelled by the ray using Single Layer fitting is  "<<TotalHorizontalDistanceSLF<<endl;
  cout<<"Now propagating the ray through the ice towards the antenna"<<endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Section for propogating the ray through the ice
  
  ////A, B and C parameters for the refractive index profile of ice
  A=A_ice;
  B=B_ice;
  C=C_ice;
  
  ////Set the starting depth of the ray for propogation to at the ice surface
  double StartDepth=0.0;
  ////Since we have the starting height of the ice layer we can find out the refactive index of air at that height from data using spline interpolation
  //double Start_nh=gsl_spline_eval(spline, IceLayerHeight, accelerator);
  Start_nh=gsl_spline_eval(spline, IceLayerHeight, accelerator);
  ////Set the stopping depth of the ray for propogation to be the depth of the antenna
  StopHeight=AntennaDepth;
  ////Set the initial launch angle or the angle of incidence of the ray
  StartAngle=IncidentAngleonIce;
  //cout<<"Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingDepth="<<StartDepth<<" ,StoppingDepth="<<AntennaDepth<<" ,RayLaunchAngle="<<StartAngle<<endl;
  
  ////Get the hit parameters from the function. The output is:
  //// How much horizontal distance did the ray travel through ice to hit the antenna
  //// The angle of reciept/incidence at the end at the antenna
  //// The value of the L parameter for whole atmosphere fit
  GetHitPar=GetLayerHitPointPar(Start_nh, A, B, C, AntennaDepth,StartDepth, StartAngle, 0);

  ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
  double TotalHorizontalDistanceinIce=GetHitPar[0];
  double RecievdAngleinIce=GetHitPar[1];
  double LvalueIce=GetHitPar[2];

  delete[] GetHitPar;
  
  cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceinIce<<endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice
  
  TMultiGraph *mg=new TMultiGraph();
  
  TGraph *grRefracted=new TGraph();
  TGraph *grStraight=new TGraph();
  TGraph *grResidual=new TGraph();
  TGraph *grAirIce=new TGraph();

  ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
  double StraightLine_slope=tan(pi/2-RayLaunchAngle*(pi/180));
  double StraightLine_y_intercept=TxHeight;
  //cout<<"StraightLine slope "<<StraightLine_slope<<" ,StraightLine_y_intercept "<<StraightLine_y_intercept<<endl;
  
  ////Define ray variables for plotting and/or storing ray path as it comes down from the atmosphere
  double Refracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray
  double LastRefracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray in the previous iteration
  double LastHeight=0;////Y coordinate variable which stores the Height of the ray in the previous iteration
  double LayerHoriOffset=0;////X axes layer offset that needs to be calculated to align the layers with each other
  double LayerStartHeight=0;////The starting height for the propagation in the layer
  double LayerStopHeight=0;////The stopping height for the propagation in the layer
  int ipoints=0;////variable for counting the total number of samples that make up the ray path

  ////Start looping over the layers to trace out the ray
  for(int il=0;il<MaxLayers-SkipLayersAbove-SkipLayersBelow;il++){

    ////Get and Set the A,B,C and L parameters for the layer
    struct fDnfR_params params2 = {layerAs[il], layerBs[il], layerCs[il], layerLs[il]};

    if(il==0){      
      ////If this is the first layer then set the start height to be the height of the source
      LayerStartHeight=TxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the next layer or the end height of the previous layer
      LayerStartHeight=LastHeight;
    }

    if(il==MaxLayers-SkipLayersAbove-SkipLayersBelow-1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      LayerStopHeight=IceLayerHeight-1;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      LayerStopHeight=(ATMLAY[MaxLayers-SkipLayersBelow-il-1]/100)-1;
    }
    //cout<<il<<" A="<<layerAs[il]<<" ,B="<<layerBs[il]<<" ,C="<<layerCs[il]<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<endl;

    ////Start tracing out the ray as it propagates through the layer
    for(double i=LayerStartHeight;i>LayerStopHeight;i--){

      ////Calculate the x (distance) value of the ray for given y (height) value
      Refracted_x=fD(-i,&params2)-fD(-TxHeight,&params2);

      ////If the ray just started off in a new layer we need to offset the x values of the new layer so that the ray aligns with the previous layer.
      if(ipoints>0 && fabs(i-LayerStartHeight)<0.0001){
  	LayerHoriOffset=(Refracted_x-LastRefracted_x);
  	//cout<<il<<" layer offset is "<<LayerHoriOffset<<endl;
      }
      Refracted_x=Refracted_x-LayerHoriOffset;

      ////Caclulate the y value of the straight line
      double StraightLine_y=StraightLine_slope*Refracted_x+StraightLine_y_intercept;

      grAirIce->SetPoint(ipoints,Refracted_x,i);
      grRefracted->SetPoint(ipoints,Refracted_x,i);
      grStraight->SetPoint(ipoints,Refracted_x,StraightLine_y);

      ////To make sure that the residual in air is only caculated above the ice surface
      if(StraightLine_y>=IceLayerHeight){
  	grResidual->SetPoint(ipoints,Refracted_x,i-StraightLine_y);
      }

      //cout<<ipoints<<" "<<Refracted_x<<" "<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<i-StraightLine_y<<endl;
      aout<<ipoints<<" "<<Refracted_x<<" "<<i<<endl;

      ////If you want to check the the transition between different layers uncomment these lines
      //if(fabs(i-LayerStopHeight)<10){
  	//cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<endl;;
      //}
      
      //if(fabs(i-LayerStartHeight)<10){
  	//cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<endl;;
      //}
      
      ipoints++;
      LastHeight=i;
    }
    LastRefracted_x=Refracted_x;
  }
  
  ////Print out the ray path in ice too  
  struct fDnfR_params params2 = {A_ice, B_ice, C_ice, LvalueIce};
  for(int i=0;i>-(AntennaDepth+1);i--){
    grAirIce->SetPoint(ipoints,TotalHorizontalDistance-fD((double)i,&params2)+fD(0,&params2),(double)i+IceLayerHeight);
    aout<<ipoints<<" "<<TotalHorizontalDistance-fD((double)i,&params2)+fD(0,&params2)<<" "<<(double)i+IceLayerHeight<<endl;
    ipoints++;
  }

  grRefracted->SetLineColor(2);//red
  grStraight->SetLineColor(4);//blue
  grRefracted->SetLineWidth(5);
  grStraight->SetLineWidth(2);
  
  grStraight->SetTitle("Straight Line with same launch angle");
  grRefracted->SetTitle("Refracted Ray");

  TString mgtitle="Launch Angle=";
  mgtitle+=180-RayLaunchAngle;
  mgtitle+=" deg where 0 deg is straight down;Distance (m); Height (m)";
  mg->SetTitle(mgtitle);
  grResidual->SetTitle("Difference of Blue with Red;Distance (m); Height (m)");
  mg->Add(grRefracted);
  mg->Add(grStraight);

  TCanvas *ca=new TCanvas("ca","ca");
  ca->Divide(1,2);
  ca->cd(1);
  mg->Draw("AL");
  ca->cd(1)->BuildLegend();

  ca->cd(2);
  grResidual->Draw("AL");
  //grResidual->Fit("pol1","","");

  grAirIce->SetTitle("RayPath through Air and Ice;Distance (m); Height (m)");
  grAirIce->SetMarkerStyle(20);
  TCanvas *cb=new TCanvas("cb","cb");
  cb->cd();
  grAirIce->Draw("ALP");

  delete accelerator;
  delete spline;
  flattened_h_data.clear();
  flattened_nh_data.clear();

  timestamp_t t1 = get_timestamp();
  
  double secs = (t1 - t0) / 1000000.0L;
  cout<<"total time taken by the script: "<<secs<<" s"<<endl;

}
