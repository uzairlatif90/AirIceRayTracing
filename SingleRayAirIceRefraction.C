#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstdlib>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>

#include <sys/time.h>

const double pi=4.0*atan(1.0); /**< Gives back value of Pi */
const double spedc=299792458.0; /**< Speed of Light in m/s */

////Define std::vectors to store data from the file
std::vector <std::vector <double>> nh_data;////n(h) refractive index profile of the atmosphere as a function of height
std::vector <std::vector <double>> lognh_data;////log(n(h)-1) log of the refractive index profile of the atmosphere as a function of height subtracted by 1
std::vector <std::vector <double>> h_data;////height data

////Define Arrays for storing values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file
double ATMLAY[5];
double abc[5][3];

////define dummy variables which will be filled in later after fitting
double C_air[5];
double B_air[5];

////define variables which are going to be used by GSL for linear interpolation
gsl_interp_accel * accelerator;
gsl_spline *spline;

////The variable which will store the max layers available in an atmosphere model
int MaxLayers=0;

////This Function reads in the values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file. The a,b and c values are mass overburden values and are not required in this code.
int readATMpar(){

  ////Open the file
  std::ifstream ain("Atmosphere.dat");
  
  int n1=0;////variable for counting total number of data points
  std::string line;
  double dummya[5]={0,0,0,0,0};////temporary variable for storing data values from the file
  
  //Check if file is open and store data
  if(ain.is_open()){
    
    while (getline(ain,line)){

      if(n1<4){////only read in the lines which contain the ATMLAY and a,b and c values in the file
	ain>>dummya[0]>>dummya[1]>>dummya[2]>>dummya[3]>>dummya[4];
	//std::cout<<n1<<" "<<dummya[0]<<" , "<<dummya[1]<<" , "<<dummya[2]<<" , "<<dummya[3]<<" , "<<dummya[4]<<std::endl;
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
  std::ifstream ain("Atmosphere.dat");
  ain.precision(10); 

  int n1=0;////variable for counting total number of data points
  int layer=0;
  std::string line;

  ////Ignore the lines containing ATMLAY and a,b and c values.
  for(int i=0; i<5; i++){ ain.ignore(256,'\n'); }
  
  ////Check if file is open and store data
  if(ain.is_open()){
    ////define dummy/temporary variables for storing data
    double dummy1,dummy2;
    ////define dummy/temporary std::vectors for storing data.
    std::vector <double> temp1,temp2,temp3;
    
    while (getline(ain,line)){
      ain>>dummy1>>dummy2;
      
      if(dummy1>-1){////start storing height at above and equal to 0 m
	////push in the height values for a single layer in the temporary std::vector
	temp1.push_back(dummy1);
	temp2.push_back(dummy2);
	temp3.push_back(log(dummy2-1));
	
	if(dummy1*100>=ATMLAY[layer]){////change the layer once the data of all the heights of that layer has been read in
	  if(layer>0){////now since the layer has finished and the temporary std::vectors have been filled in. Now we push the std::vectors in the main 2d height and refractice index std::vectors
	    h_data.push_back(temp1);
	    nh_data.push_back(temp2);
	    lognh_data.push_back(temp3);
	    ////clear the std::vectors now for storing the next layer
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
      ////clear the std::vectors now for storing the next layer
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

  MaxLayers=h_data.size();////store the total number of layers present in the data
  
  return 0;
}

////Set the value of the asymptotic parameter of the ice refractive index model
const double A_ice=1.78;

////Get the value of the B parameter for the ice refractive index model
double GetB_ice(double z){
  double zabs=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

////Get the value of the C parameter for the ice refractive index model
double GetC_ice(double z){
  double zabs=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth in ice
double Getnz_ice(double z){
  z=fabs(z);
  return A_ice+GetB_ice(z)*exp(-GetC_ice(z)*z);
}

////Set the value of the asymptotic parameter of the air refractive index model
const double A_air=1.00;


int FillInAirRefractiveIndex(){
  
  double N0=0;
  for(int ilayer=0;ilayer<MaxLayers;ilayer++){
    double hlow=ATMLAY[ilayer]/100;
    C_air[ilayer]=1.0/(abc[ilayer][2]/100);
    if(ilayer>0){
      N0=A_air+B_air[ilayer-1]*exp(-hlow*C_air[ilayer-1]);
    }
    if(ilayer==0){
      N0=gsl_spline_eval(spline, 0, accelerator);
    }
    B_air[ilayer]=((N0-1)/exp(-hlow*C_air[ilayer]));
  }

  return 0;   
}

////Get the value of the B parameter for the air refractive index model
double GetB_air(double z){
  double zabs=fabs(z);
  double B=0;
  int whichlayer=0;
 
  for(int ilayer=0;ilayer<MaxLayers-1;ilayer++){

    if(zabs<ATMLAY[ilayer+1]/100 && zabs>=ATMLAY[ilayer]/100){
      whichlayer=ilayer;
      ilayer=100;
    }  
  }
  if(zabs>=ATMLAY[MaxLayers-1]/100){
    whichlayer=MaxLayers-1;
  }
 
  B=B_air[whichlayer];
  return B;
}

////Get the value of the C parameter for the air refractive index model
double GetC_air(double z){
  double zabs=fabs(z);
  double C=0;
  int whichlayer=0;
  
  for(int ilayer=0;ilayer<MaxLayers-1;ilayer++){
    if(zabs<ATMLAY[ilayer+1]/100 && zabs>=ATMLAY[ilayer]/100){
      whichlayer=ilayer;
      ilayer=100;
    }
  }
  
  if(zabs>=ATMLAY[MaxLayers-1]/100){
    whichlayer=MaxLayers-1;
  }
  C=C_air[whichlayer];
  return C;
}

////Get the value of refractive index model for a given depth in air
double Getnz_air(double z){
  z=fabs(z);
  return A_air+GetB_air(z)*exp(-GetC_air(z)*z);
}

////Use GSL minimiser which uses Brent's Method to find root for a given function
double FindFunctionRoot(gsl_function F,double x_lo, double x_hi,const gsl_root_fsolver_type *T)
{
  int status;
  int iter = 0, max_iter = 100;
  gsl_root_fsolver *s;
  double r = 0;
  
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
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
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.000001);

      if (status == GSL_SUCCESS){
	//printf ("Converged:");
	//printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

/////Functions used for Raytracing in Ice using the analytical solution

////Analytical solution describing the ray path in ice
struct fDnfR_params { double a, b, c, l; };
double fDnfR(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
struct fdxdz_params { double lang, z0,z1; int airorice;};
double fdxdz(double x,void *params){
  
  struct fdxdz_params *p= (struct fdxdz_params *) params;
  double Lang = p->lang;
  double Z0 = p->z0;
  double Z1 = p->z1;
  int AirOrIce = p->airorice;
  
   double output=0;
  if(AirOrIce==0){
    output=tan(asin( (Getnz_ice(Z0)*sin(x))/Getnz_ice(Z1) ) ) - tan(Lang);
  }
  
  if(AirOrIce==1){
    output=tan(asin( (Getnz_air(Z0)*sin(x))/Getnz_air(Z1) ) ) - tan(Lang);
  }
  return output;
}

////The function used to calculate ray propogation time in ice
struct ftimeD_params { double a, b, c, speedc,l; int airorice; };
double ftimeD(double x,void *params){
  
  struct ftimeD_params *p= (struct ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;
  int AirOrIce=p->airorice;

  double result=0;
  if(AirOrIce==0){//in ice
    result=(1.0/(Speedc*C*sqrt(pow(Getnz_ice(x),2)-L*L)))*(pow(Getnz_ice(x),2)-L*L+(C*x-log(A*Getnz_ice(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz_ice(x),2)-L*L)))*(A*A*sqrt(pow(Getnz_ice(x),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(Getnz_ice(x),2)-L*L)*log(Getnz_ice(x)+sqrt(pow(Getnz_ice(x),2)-L*L)) );
  }
  if(AirOrIce==1){//in air
    result=(1.0/(Speedc*C*sqrt(pow(Getnz_air(x),2)-L*L)))*(pow(Getnz_air(x),2)-L*L+(C*x-log(A*Getnz_air(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz_air(x),2)-L*L)))*(A*A*sqrt(pow(Getnz_air(x),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(Getnz_air(x),2)-L*L)*log(Getnz_air(x)+sqrt(pow(Getnz_air(x),2)-L*L)) );
  }
  
  return result;
}

////Get the distance on the 2ndLayer at which the ray should hit given an incident angle such that it hits an target at depth of z0 m in the second layer.
//// n_layer1 is the refractive index value of the previous layer at the boundary of the two mediums
//// A,B and C values are the values required for n(z)=A+B*exp(Cz) for the second layer
//// TxDepth is the starting height or depth
//// RxDepth is the final height or depth
//// AirOrIce variable is used to determine whether we are working in air or ice as that sets the range for the GSL root finder.
double *GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){
  double *output=new double[4];
  
  double x0=0;////Starting horizontal point of the ray. Always set at zero
  double x1=0;////Variable to store the horizontal distance that will be traveled by the ray
  
  double ReceiveAngle=0;////Angle from the vertical at which the target will recieve the ray
  double Lvalue=0;//// L parameter of the ray for that layer
  double RayTimeIn2ndLayer=0;////Time of propagation in 2ndLayer 
  double AngleOfEntryIn2ndLayer=0;////Angle at which the ray enters the layer

  double SurfaceRayIncidentAngle=IncidentAng*(pi/180.0);////Angle at which the ray is incident on the second layer
  double RayAngleInside2ndLayer=0;////Use Snell's Law to find the angle of transmission in the 2ndlayer

  double A=0;
  double nzRx=0;
  double nzTx=0;
  double GSLFnLimit=0;

  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    A=A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    A=A_air;
    nzRx=Getnz_air(RxDepth);
    nzTx=Getnz_air(TxDepth);
  }

  ////LimitAngle sets a limit on the range to which the GSL minimisation will work. This limit comes from the fact that in fdxdx() you have tan(asin(x)) which goes to infinity at x=1. In our case x=(nz(Z0)*sin(Angle))/nz(Z1) . Solving for Angle gives us our limit.
  double LimitAngle=asin(nzTx/nzRx);
  
  GSLFnLimit=LimitAngle;
  RayAngleInside2ndLayer=asin((n_layer1/nzTx)*sin(SurfaceRayIncidentAngle));////Use Snell's Law to find the angle of transmission in the 2ndlayer
  
  ////calculate the angle at which the target receives the ray
  gsl_function F1;
  struct fdxdz_params params1 = {RayAngleInside2ndLayer, RxDepth, TxDepth, AirOrIce};
  F1.function = &fdxdz;
  F1.params = &params1;
  ReceiveAngle=FindFunctionRoot(F1,0*(pi/180),GSLFnLimit,gsl_root_fsolver_bisection);
  //std::cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/pi)<<" deg"<<std::endl;
  
  ////calculate the distance of the point of incidence on the 2ndLayer surface and also the value of the L parameter of the solution
  Lvalue=nzRx*sin(ReceiveAngle);
  struct fDnfR_params params2a;
  struct fDnfR_params params2b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params2a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), Lvalue};
    params2b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), Lvalue};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params2a = {A, GetB_air(RxDepth), -GetC_air(RxDepth), Lvalue};
    params2b = {A, GetB_air(TxDepth), -GetC_air(TxDepth), Lvalue};
  }
  x1=+fDnfR(RxDepth,&params2a)-fDnfR(TxDepth,&params2b);
  if(AirOrIce==1){
    x1*=-1;
  }
  
  ////calculate the propagation time in 2ndLayer
  struct ftimeD_params params3a;
  struct ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), spedc, Lvalue,0};
    params3b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, GetB_air(RxDepth), -GetC_air(RxDepth), spedc, Lvalue,1};
    params3b = {A, GetB_air(TxDepth), -GetC_air(TxDepth), spedc, Lvalue,1};
  }
  RayTimeIn2ndLayer=+ftimeD(RxDepth,&params3a)-ftimeD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  //std::cout<<"The propagation time in 2ndLayer is: "<<RayTimeIn2ndLayer<<" s"<<std::endl;

  ///////calculate the initial angle when the ray enters the 2ndLayer. This should be the same as RayAngleInside2ndLayer. This provides a good sanity check to make sure things have worked out.
  // gsl_function F4;
  // double result, abserr;
  // F4.function = &fDnfR;
  // F4.params = &params2b;
  // gsl_deriv_central (&F4, TxDepth, 1e-8, &result, &abserr);
  // AngleOfEntryIn2ndLayer=atan(result)*(180.0/pi);
  // if(TxDepth==RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=180-ReceiveAngle;
  // }
  // if(TxDepth!=RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=90;
  // }
  //std::cout<<"AngleOfEntryIn2ndLayer= "<<AngleOfEntryIn2ndLayer<<" ,RayAngleInside2ndLayer="<<RayAngleInside2ndLayer*(180/pi)<<std::endl;

  output[0]=x1;
  output[1]=ReceiveAngle*(180/pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  
  return output;
}

////This function flattens out 2d std::vectors into 1d std::vectors
std::vector<double> flatten(const std::vector<std::vector<double>>& v) {
    size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<double> result;
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

int MakeAtmosphere(){
   
  ////Fill in the n(h) and h arrays and ATMLAY and a,b and c (these 3 are the mass overburden parameters) from the data file
  readATMpar();
  readnhFromFile();
  
  ////Flatten out the height and the refractive index std::vectors to be used for setting the up the spline interpolation.
  std::vector <double> flattened_h_data=flatten(h_data);
  std::vector <double> flattened_nh_data=flatten(nh_data);
  
  ////Set up the GSL cubic spline interpolation. This used for interpolating values of refractive index at different heights.
  accelerator =  gsl_interp_accel_alloc();
  spline = gsl_spline_alloc (gsl_interp_cspline,flattened_h_data.size());
  gsl_spline_init(spline, flattened_h_data.data(), flattened_nh_data.data(), flattened_h_data.size());
  
  FillInAirRefractiveIndex();
  // delete accelerator;
  // delete spline;

  return 0;
}

int main(int argc, char **argv){
  
  ////Name of the main function if you want to run it with ROOT
  //void SingleRayAirIceRefraction(){

  ////For recording how much time the process took
  timestamp_t t0 = get_timestamp();
  
  if(argc==1){
    std::cout<<"No Extra Command Line Argument Passed Other Than Program Name"<<std::endl;
    std::cout<<"Example run command: ./SingleRayAirIceRefraction 200 170 20000 3000"<<std::endl;
    std::cout<<"Here 200 is Antenna Depth in m, 170 is the ray launch angle (where is 0 vertically up) in deg, 20000 is the Tx Height in m and 3000 is Ice Layer Height in m"<<std::endl;
    return 0;
  }
  if(argc<5){
    std::cout<<"More Arguments needed!"<<std::endl;
    std::cout<<"Example run command: ./SingleRayAirIceRefraction 200 170 20000 3000"<<std::endl;
    std::cout<<"Here 200 is Antenna Depth in m, 170 is the ray launch angle (where is 0 vertically up) in deg, 20000 is the Tx Height in m and 3000 is Ice Layer Height in m"<<std::endl;
    return 0;
  }
  if(argc==5){
    std::cout<<"Antenna Depth is set at "<<atof(argv[1])<<" m, The Ray Launch Angle is set at "<<atof(argv[2])<<" deg, Tx Height is set at "<<atof(argv[3]) <<" m, Ice Layer Height is set as "<<atof(argv[4])<<" m"<<std::endl;
  } 
  if(argc>5){
    std::cout<<"More Arguments than needed!"<<std::endl;
    std::cout<<"Example run command: ./SingleRayAirIceRefraction 200 170 20000 3000"<<std::endl;
    std::cout<<"Here 200 is Antenna Depth in m, 170 is the ray launch angle (where is 0 vertically up) in deg, 20000 is the Tx Height in m and 3000 is Ice Layer Height in m"<<std::endl;
    return 0;
  }

  MakeAtmosphere();

  double AntennaDepth=atof(argv[1]);//200;////Depth of antenna in the ice
  double RayLaunchAngle=atof(argv[2]);//170;////Initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 is vertically down
  double AirTxHeight=atof(argv[3]);//20000;////Height of the source
  double IceLayerHeight=atof(argv[4]);//3000;////Height where the ice layer starts off

  if(AirTxHeight>h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]){
    ////Maximum height available with the refractive index data
    std::cout<<"Tx Height is set higher than maximum available height for atmospheric refractive index which is "<<h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]<<std::endl;
    AirTxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];
    std::cout<<"Setting Tx Height to be the maximum available height"<<std::endl;
  }

  if(RayLaunchAngle<=90){
    std::cout<<"RayLaunchAngle has been set at "<<RayLaunchAngle<<" deg which is outside of the allowed range of 90 deg <RayLaunchAngle< 180 deg"<<std::endl;
    RayLaunchAngle=135;
    std::cout<<"Setting RayLaunchAngle at"<<RayLaunchAngle<<std::endl;
  }

  ////Print out the ray path x and y values in a file
  std::ofstream aout("RayPathinAirnIce.txt");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Section for propogating the ray through the atmosphere
  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
    //std::cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<std::endl;
    if(AirTxHeight<ATMLAY[ilayer]/100 && AirTxHeight>=ATMLAY[ilayer-1]/100){
      std::cout<<"Tx Height is in this layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<std::endl;
      ilayer=-100;
    }
    if(ilayer>-1){
      skiplayer++;
    }
  }
  int SkipLayersAbove=skiplayer;
  //std::cout<<"The total number of layers that need to be skipped from above is "<<skiplayer<<std::endl;

  ////Find out how many atmosphere layers are below the ice height which we do not need
  skiplayer=0;
  for(int ilayer=0;ilayer<MaxLayers;ilayer++){
    //std::cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer+1]/100<<std::endl;
    if(IceLayerHeight>=ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
      std::cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<std::endl;
      ilayer=100;
    }
    if(ilayer<MaxLayers){
      skiplayer++;
    }
  }
  int SkipLayersBelow=skiplayer;
  //std::cout<<"The tota number of layers that need to be skipped from below is "<<skiplayer<<std::endl;
  
  ////Define variables for ray propogation through mutliple layers in the atmosphere
  double Start_nh=0;
  double StartHeight=0;
  double StopHeight=0;
  double StartAngle=0;
  double TotalHorizontalDistance=0;
  std::vector <double> layerLs;////std::vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them

  ////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
  for(int ilayer=MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
    
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=ATMLAY[ilayer+1]/100-0.0001;
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
    //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<std::endl;

    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    double* GetHitPar=GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
    TotalHorizontalDistance+=GetHitPar[0];
    StartAngle=GetHitPar[1];

    ////Store in the values of A,B,C and L for tha layer
    layerLs.push_back(GetHitPar[2]);
    
    ////dont forget to delete the pointer!
    delete []GetHitPar;
  }

  double IncidentAngleonIce=StartAngle;
  std::cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<std::endl;
  std::cout<<"Now treating the atmosphere refrative index profile as a single layer and fitting it and propogating the ray"<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Section for propogating the ray through the ice 
  
  ////Set the starting depth of the ray for propogation to at the ice surface
  double StartDepth=0.0;

  ////Since we have the starting height of the ice layer we can find out the refactive index of air at that height from data using spline interpolation
  Start_nh=gsl_spline_eval(spline, IceLayerHeight, accelerator);
  
  ////Set the stopping depth of the ray for propogation to be the depth of the antenna
  StopHeight=AntennaDepth;
  ////Set the initial launch angle or the angle of incidence of the ray
  StartAngle=IncidentAngleonIce;
  //std::cout<<"Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingDepth="<<StartDepth<<" ,StoppingDepth="<<AntennaDepth<<" ,RayLaunchAngle="<<StartAngle<<std::endl;
  
  ////Get the hit parameters from the function. The output is:
  //// How much horizontal distance did the ray travel through ice to hit the antenna
  //// The angle of reciept/incidence at the end at the antenna
  //// The value of the L parameter for whole atmosphere fit
  double *GetHitPar=GetLayerHitPointPar(Start_nh, AntennaDepth,StartDepth, StartAngle, 0);

  ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
  double TotalHorizontalDistanceinIce=GetHitPar[0];
  double RecievdAngleinIce=GetHitPar[1];
  double LvalueIce=GetHitPar[2];

  delete[] GetHitPar;
  
  std::cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceinIce<<std::endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice
  
  ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
  double StraightLine_slope=tan(pi/2-RayLaunchAngle*(pi/180));
  double StraightLine_y_intercept=AirTxHeight;
  //std::cout<<"StraightLine slope "<<StraightLine_slope<<" ,StraightLine_y_intercept "<<StraightLine_y_intercept<<std::endl;
  
  ////Define ray variables for plotting and/or storing ray path as it comes down from the atmosphere
  double Refracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray
  double LastRefracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray in the previous iteration
  double LastHeight=0;////Y coordinate variable which stores the Height of the ray in the previous iteration
  double LayerHoriOffset=0;////X axes layer offset that needs to be calculated to align the layers with each other
  double LayerStartHeight=0;////The starting height for the propagation in the layer
  double LayerStopHeight=0;////The stopping height for the propagation in the layer
  int ipoints=0;////variable for counting the total number of samples that make up the ray path

  ////Get and Set the A,B,C and L parameters for the layer
  struct fDnfR_params params2a;
  struct fDnfR_params params2b;
  
  ////Start looping over the layers to trace out the ray
  for(int il=0;il<MaxLayers-SkipLayersAbove-SkipLayersBelow;il++){
    
    if(il==0){
      ////If this is the first layer then set the start height to be the height of the source
      LayerStartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the next layer or the end height of the previous layer
      LayerStartHeight=LastHeight-0.00001;
    }
    
    if(il==MaxLayers-SkipLayersAbove-SkipLayersBelow-1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      LayerStopHeight=IceLayerHeight-1;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      LayerStopHeight=(ATMLAY[MaxLayers-SkipLayersAbove-SkipLayersBelow-il-1]/100)-1;
    }
    
    //std::cout<<il<<" A="<<layerAs[il]<<" ,B="<<layerBs[il]<<" ,C="<<layerCs[il]<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<" ,LayerStartHeight="<<LayerStartHeight<<" ,LayerStopHeight="<<LayerStopHeight<<std::endl;
    //std::cout<<" new layer "<<std::endl;
    ////Start tracing out the ray as it propagates through the layer
    for(double i=LayerStartHeight-1;i>LayerStopHeight+1;i--){
      
      ////Get and Set the A,B,C and L parameters for the layer
      params2a = {A_air, GetB_air(-i), GetC_air(-i), layerLs[il]};
      params2b = {A_air, GetB_air(-(i)), GetC_air(-(i)), layerLs[il]};
      
      ////Calculate the x (distance) value of the ray for given y (height) value
      Refracted_x=fDnfR(-i,&params2a)-fDnfR(-(LayerStartHeight),&params2b)+LastRefracted_x;
      
      ////If the ray just started off in a new layer we might need to offset the x values of the new layer so that the ray aligns with the previous layer.
      // if(ipoints>0 && fabs(i-LayerStartHeight-1)<2){
      //   LayerHoriOffset=(Refracted_x-LastRefracted_x);
      //   std::cout<<il<<" layer offset is "<<LayerHoriOffset<<std::endl;
      // }
      // Refracted_x=Refracted_x-LayerHoriOffset;

      ////Caclulate the y value of the straight line
      double StraightLine_y=StraightLine_slope*Refracted_x+StraightLine_y_intercept;      
      
      //std::cout<<ipoints<<" "<<Refracted_x<<" "<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<i-StraightLine_y<<std::endl;
      aout<<ipoints<<" "<<Refracted_x<<" "<<i<<std::endl;

      // ////If you want to check the the transition between different layers uncomment these lines
      // if(fabs(i-LayerStopHeight)<10){
      //   std::cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<" "<<GetB_air(-i)<<" "<<GetC_air(-i)<<" "<<layerLs[il]<<" "<<fDnfR(-i,&params2a)<<" "<<-fDnfR(-(LayerStartHeight),&params2b)<<" "<<LayerStartHeight<<" "<<LastRefracted_x<<" "<<GetB_air(-LayerStartHeight)<<" "<<GetC_air(-LayerStartHeight)<<std::endl;
      // }
      
      // if(fabs(i-LayerStartHeight)<10){
      //   std::cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<" "<<GetB_air(-i)<<" "<<GetC_air(-i)<<" "<<layerLs[il]<<" "<<fDnfR(-i,&params2a)<<" "<<-fDnfR(-(LayerStartHeight),&params2b)<<" "<<LayerStartHeight<<" "<<LastRefracted_x<<" "<<GetB_air(-LayerStartHeight)<<" "<<GetC_air(-LayerStartHeight)<<std::endl;
      // }
      
      ipoints++;
      LastHeight=i;
    }
    LastRefracted_x=Refracted_x;
  }      
  
  ////Print out the ray path in ice too  
  struct fDnfR_params params3a;
  struct fDnfR_params params3b;
  
  for(int i=0;i>-(AntennaDepth+1);i--){
    params3a = {A_ice, GetB_ice(i), GetC_ice(i), LvalueIce};
    params3b = {A_ice, GetB_ice(0), GetC_ice(0), LvalueIce};
    
    double refractedpath=LastRefracted_x-fDnfR((double)i,&params3a)+fDnfR(0,&params3b);
    aout<<ipoints<<" "<<refractedpath<<" "<<(double)i+IceLayerHeight<<std::endl;
    ipoints++;
  }
  
  timestamp_t t1 = get_timestamp();
  
  double secs = (t1 - t0) / 1000000.0L;
  std::cout<<"total time taken by the script: "<<secs<<" s"<<std::endl;
  
  return 0;
}
