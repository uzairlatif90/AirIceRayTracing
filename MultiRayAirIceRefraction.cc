#include "MultiRayAirIceRefraction.h"

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

////Store the maximum possible height allowed by GDAS tables
double MaxAirTxHeight=0;
////Store the minimum possible height allowed antenna tables
double MinAirTxHeight=0;

////Multidimensional vector to store data read in from the antenna tables for interpolation puprposes
std::vector<std::vector<std::vector <double>>> AllTableAllAntData;

////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
double AngleStepSize=0.5;
double LoopStartAngle=92;
double LoopStopAngle=180;
int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize)+1;

////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
double HeightStepSize=20;
double LoopStartHeight=0;
double LoopStopHeight=0;
int TotalHeightSteps=0;

////This Function reads in the values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file. The a,b and c values are mass overburden values and are not required in this code.
int MultiRayAirIceRefraction::readATMpar(){

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

  abc[4][0]=abc[3][0];
  abc[4][1]=abc[3][1];
  abc[4][2]=abc[3][2];
  
  return 0;
}

int MultiRayAirIceRefraction::readnhFromFile(){

  nh_data.clear();
  lognh_data.clear();
  h_data.clear();
  
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

  MaxLayers=h_data.size()+1;////store the total number of layers present in the data
  
  return 0;
}

////Get the value of the B parameter for the ice refractive index model
double MultiRayAirIceRefraction::GetB_ice(double z){
  //double zabs=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

////Get the value of the C parameter for the ice refractive index model
double MultiRayAirIceRefraction::GetC_ice(double z){
  //double zabs=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth in ice
double MultiRayAirIceRefraction::Getnz_ice(double z){
  double zabs=fabs(z);
  return MultiRayAirIceRefraction::A_ice+GetB_ice(zabs)*exp(-GetC_ice(zabs)*zabs);
}

int MultiRayAirIceRefraction::FillInAirRefractiveIndex(){
  
  double N0=0;
  for(int ilayer=0;ilayer<5;ilayer++){
    double hlow=ATMLAY[ilayer]/100;
    C_air[ilayer]=1.0/(abc[ilayer][2]/100);
    if(ilayer>0){
      N0=A_air+B_air[ilayer-1]*exp(-hlow*C_air[ilayer-1]);
    }
    if(ilayer==0){
      N0=gsl_spline_eval(spline, 0, accelerator);
    }
    B_air[ilayer]=((N0-1)/exp(-hlow*C_air[ilayer]));
    //std::cout<<ilayer<<" "<<B_air[ilayer]<<" "<<C_air[ilayer]<<std::endl;
  }
  
  return 0;   
}

////Get the value of the B parameter for the air refractive index model
double MultiRayAirIceRefraction::GetB_air(double z){
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
  //B=0.000333;
  return B;
}

////Get the value of the C parameter for the air refractive index model
double MultiRayAirIceRefraction::GetC_air(double z){
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
  //C=0.0000000000000001;
  return C;
}

////Get the value of refractive index model for a given depth in air
double MultiRayAirIceRefraction::Getnz_air(double z){
  double zabs=fabs(z);
  return MultiRayAirIceRefraction::A_air+GetB_air(zabs)*exp(-GetC_air(zabs)*zabs);
}



////E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R
double MultiRayAirIceRefraction::Refl_S(double thetai, double IceLayerHeight){
  double Nair=MultiRayAirIceRefraction::Getnz_air(IceLayerHeight);
  double Nice=MultiRayAirIceRefraction::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;
  
  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double RS=(num*num)/(den*den);
  if(isnan(RS)){
    RS=1;
  }
  return (RS);
}

////E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R
double MultiRayAirIceRefraction::Refl_P(double thetai, double IceLayerHeight){
  double Nair=MultiRayAirIceRefraction::Getnz_air(IceLayerHeight);
  double Nice=MultiRayAirIceRefraction::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;

  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double RP=(num*num)/(den*den);
  if(isnan(RP)){
    RP=1;
  }
  return (RP);
}


////Use GSL minimiser which uses Brent's Method to find root for a given function
double MultiRayAirIceRefraction::FindFunctionRoot(gsl_function F,double x_lo, double x_hi,const gsl_root_fsolver_type *T,double tolerance)
{
  int status;
  int iter = 0, max_iter = 100;
  //const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  //double tolerance=0.000000001;
  
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
      status = gsl_root_test_interval (x_lo, x_hi,0, tolerance);

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
double MultiRayAirIceRefraction::fDnfR(double x,void *params){
  
  struct MultiRayAirIceRefraction::fDnfR_params *p= (struct MultiRayAirIceRefraction::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
double MultiRayAirIceRefraction::fdxdz(double x,void *params){
  
  struct MultiRayAirIceRefraction::fdxdz_params *p= (struct MultiRayAirIceRefraction::fdxdz_params *) params;
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
double MultiRayAirIceRefraction::ftimeD(double x,void *params){
  
  struct MultiRayAirIceRefraction::ftimeD_params *p= (struct MultiRayAirIceRefraction::ftimeD_params *) params;
  double A = p->a;
  //double B = p->b;
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

double MultiRayAirIceRefraction::GetRayOpticalPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct MultiRayAirIceRefraction::fDnfR_params params2a;
  struct MultiRayAirIceRefraction::fDnfR_params params2b;
  ////just initialise the parameters
  params2a = {0,0,0,0};
  params2b = {0,0,0,0};
  
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
  double x1=+fDnfR(RxDepth,&params2a)-fDnfR(TxDepth,&params2b);
  if(AirOrIce==1){
    x1*=-1;
  }
  
  return x1;
}

double MultiRayAirIceRefraction::GetRayPropagationTime(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct MultiRayAirIceRefraction::ftimeD_params params3a;
  struct MultiRayAirIceRefraction::ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
    params3b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, GetB_air(RxDepth), -GetC_air(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
    params3b = {A, GetB_air(TxDepth), -GetC_air(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
  }
  double RayTimeIn2ndLayer=+ftimeD(RxDepth,&params3a)-ftimeD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  
  return RayTimeIn2ndLayer;
}

////Get the distance on the 2ndLayer at which the ray should hit given an incident angle such that it hits an target at depth of z0 m in the second layer.
//// n_layer1 is the refractive index value of the previous layer at the boundary of the two mediums
//// A,B and C values are the values required for n(z)=A+B*exp(Cz) for the second layer
//// TxDepth is the starting height or depth
//// RxDepth is the final height or depth
//// AirOrIce variable is used to determine whether we are working in air or ice as that sets the range for the GSL root finder.
double *MultiRayAirIceRefraction::GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){

  double *output=new double[4];
  
  //double x0=0;////Starting horizontal point of the ray. Always set at zero
  double x1=0;////Variable to store the horizontal distance that will be traveled by the ray
  
  double ReceiveAngle=0;////Angle from the vertical at which the target will recieve the ray
  double Lvalue=0;//// L parameter of the ray for that layer
  double RayTimeIn2ndLayer=0;////Time of propagation in 2ndLayer 
  //double AngleOfEntryIn2ndLayer=0;////Angle at which the ray enters the layer

  double SurfaceRayIncidentAngle=IncidentAng*(MultiRayAirIceRefraction::pi/180.0);////Angle at which the ray is incident on the second layer
  double RayAngleInside2ndLayer=0;////Use Snell's Law to find the angle of transmission in the 2ndlayer

  double A=0;
  double nzRx=0;
  double nzTx=0;
  double GSLFnLimit=0;

  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    A=MultiRayAirIceRefraction::A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    A=MultiRayAirIceRefraction::A_air;
    nzRx=Getnz_air(RxDepth);
    nzTx=Getnz_air(TxDepth);
  }

  ////LimitAngle sets a limit on the range to which the GSL minimisation will work. This limit comes from the fact that in fdxdx() you have tan(asin(x)) which goes to infinity at x=1. In our case x=(nz(Z0)*sin(Angle))/nz(Z1) . Solving for Angle gives us our limit.
  double LimitAngle=asin(nzTx/nzRx);
  
  GSLFnLimit=LimitAngle;
  RayAngleInside2ndLayer=asin((n_layer1/nzTx)*sin(SurfaceRayIncidentAngle));////Use Snell's Law to find the angle of transmission in the 2ndlayer
  
  ////calculate the angle at which the target receives the ray
  gsl_function F1;
  struct MultiRayAirIceRefraction::fdxdz_params params1 = {RayAngleInside2ndLayer, RxDepth, TxDepth, AirOrIce};
  F1.function = &fdxdz;
  F1.params = &params1;
  ReceiveAngle=FindFunctionRoot(F1,0*(MultiRayAirIceRefraction::pi/180),GSLFnLimit, gsl_root_fsolver_falsepos,0.000000001);
  //std::cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/MultiRayAirIceRefraction::pi)<<" deg"<<std::endl;
  
  ////calculate the distance of the point of incidence on the 2ndLayer surface and also the value of the L parameter of the solution
  Lvalue=nzRx*sin(ReceiveAngle);

  x1=GetRayOpticalPath(A, RxDepth, TxDepth, Lvalue, AirOrIce);
  //std::cout<<"The hit point horizontal distance is from the Rx target "<<x1<<" m  on the surface"<<std::endl;
  
  ////calculate the propagation time in 2ndLayer 
  RayTimeIn2ndLayer=GetRayPropagationTime(A, RxDepth, TxDepth, Lvalue, AirOrIce);
  //std::cout<<"The propagation time in 2ndLayer is: "<<RayTimeIn2ndLayer<<" s"<<std::endl;

  ///////calculate the initial angle when the ray enters the 2ndLayer. This should be the same as RayAngleInside2ndLayer. This provides a good sanity check to make sure things have worked out.
  // gsl_function F4;
  // double result, abserr;
  // F4.function = &MultiRayAirIceRefraction::fDnfR;
  // F4.params = &params2b;
  // gsl_deriv_central (&F4, TxDepth, 1e-8, &result, &abserr);
  // AngleOfEntryIn2ndLayer=atan(result)*(180.0/MultiRayAirIceRefraction::pi);
  // if(TxDepth==RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=180-ReceiveAngle;
  // }
  // if(TxDepth!=RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=90;
  // }
  //std::cout<<"AngleOfEntryIn2ndLayer= "<<AngleOfEntryIn2ndLayer<<" ,RayAngleInside2ndLayer="<<RayAngleInside2ndLayer*(180/MultiRayAirIceRefraction::pi)<<std::endl;

  output[0]=x1;
  output[1]=ReceiveAngle*(180/MultiRayAirIceRefraction::pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  
  return output;
}

////This function flattens out 2d std::vectors into 1d std::vectors
std::vector<double> MultiRayAirIceRefraction::flatten(const std::vector<std::vector<double>>& v) {
    size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<double> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

////Get Propogation parameters for ray propagating in air
double * MultiRayAirIceRefraction::GetAirPropagationPar(double LaunchAngleAir, double AirTxHeight, double IceLayerHeight){
  double *output=new double[4*MaxLayers+1];
  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
    if(AirTxHeight<ATMLAY[ilayer]/100 && AirTxHeight>ATMLAY[ilayer-1]/100){
      //cout<<"Tx Height is in this layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<endl;
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
    if(IceLayerHeight>ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
      //cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
      ilayer=100;
    }
    if(ilayer<MaxLayers){
      skiplayer++;
    }
  }
  int SkipLayersBelow=skiplayer;
  
  double StartAngle=0;
  double StartHeight=0;
  double Start_nh=0;
  double StopHeight=0;

  std::vector <double> TotalHorizontalDistance;
  std::vector <double> ReceiveAngle;
  std::vector <double> Lvalue;
  std::vector <double> PropagationTime;

  int ipoints=0;
  for(int ilayer=MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
    
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=ATMLAY[ilayer+1]/100-0.00001;
    }
    
    ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
    Start_nh=Getnz_air(StartHeight);
    
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
      StartAngle=180-LaunchAngleAir;
    }
    //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;
    
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      double* GetHitPar=GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
      TotalHorizontalDistance.push_back(GetHitPar[0]);
      ReceiveAngle.push_back(GetHitPar[1]);
      Lvalue.push_back(GetHitPar[2]);
      PropagationTime.push_back(GetHitPar[3]);
      StartAngle=GetHitPar[1];
      delete []GetHitPar;  
    }
    if(ilayer<MaxLayers-SkipLayersAbove-1){
      Lvalue.push_back(Lvalue[0]);
      double nzStopHeight=Getnz_air(StopHeight);
      double RecAng=asin(Lvalue[0]/nzStopHeight);
      RecAng=RecAng*(180/MultiRayAirIceRefraction::pi);
      ReceiveAngle.push_back(RecAng);
      double THD=GetRayOpticalPath(A_air, StopHeight, StartHeight, Lvalue[0], 1);
      TotalHorizontalDistance.push_back(THD);
      double PropTime=GetRayPropagationTime(A_air, StopHeight, StartHeight, Lvalue[0], 1);
      PropagationTime.push_back(PropTime);
      StartAngle=RecAng;
    }
    //cout<<ilayer<<" "<<TotalHorizontalDistance[ipoints]<<" "<<ReceiveAngle[ipoints]<<" "<<Lvalue[ipoints]<<" "<<PropagationTime[ipoints]<<endl;
    
    ipoints++;
    ////dont forget to delete the pointer!
    
  }
  
  for(int i=0;i<Lvalue.size();i++){
    output[4*i+0]=TotalHorizontalDistance[i];
    output[4*i+1]=ReceiveAngle[i];
    output[4*i+2]=Lvalue[i];
    output[4*i+3]=PropagationTime[i];
  }
  output[4*MaxLayers]=Lvalue.size();
  return output;
}

////Get Propogation parameters for ray propagating in ice
double * MultiRayAirIceRefraction::GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth, double Lvalue){
  double *output=new double[4];

  //double StartAngle=IncidentAngleonIce;
  double StartDepth=0.0;
  double StopDepth=AntennaDepth;
  double nzStopDepth=Getnz_ice(StopDepth);
  
  double TotalHorizontalDistance=GetRayOpticalPath(A_ice, StopDepth, StartDepth, Lvalue, 0);
  double ReceiveAngle=asin(Lvalue/nzStopDepth)*(180/MultiRayAirIceRefraction::pi);
  double PropagationTime=GetRayPropagationTime(A_ice, StopDepth, StartDepth, Lvalue, 0);

  output[0]=TotalHorizontalDistance;
  output[1]=ReceiveAngle;
  output[2]=Lvalue;
  output[3]=PropagationTime;

  return output;
}


////This function will be minimised to get the launch angle in air. It looks at the total horizontal distance and tries to find a ray in air and ice which can cover that horizontal distance perfectly. It takes into account the discontinuity in refractive index at the air ice boundary
double MultiRayAirIceRefraction::MinimizeforLaunchAngle(double x, void *params){
  
  struct MultiRayAirIceRefraction::MinforLAng_params *p= (struct MultiRayAirIceRefraction::MinforLAng_params *) params;
  double AirTxHeight = p->airtxheight;
  double IceLayerHeight = p->icelayerheight;
  double AntennaDepth = p->antennadepth;
  double HorizontalDistance = p->horizontaldistance;
  //std::cout<<AirTxHeight<<" "<<IceLayerHeight<<" "<<AntennaDepth<<" "<<HorizontalDistance<<std::endl;

  double * GetResultsAir=GetAirPropagationPar(x,AirTxHeight,IceLayerHeight);
  double TotalHorizontalDistanceinAir=0;
  int FilledLayers=GetResultsAir[4*MaxLayers];
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[i*4];
  }
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*4];
  double Lvalue=GetResultsAir[2];
  delete [] GetResultsAir;

  double * GetResultsIce=GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth, Lvalue);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  delete [] GetResultsIce;

  double checkmin=((TotalHorizontalDistanceinIce + TotalHorizontalDistanceinAir) - HorizontalDistance);
  
  return checkmin;
  
}

///This function loads in the GDAS atmosphere file. It calls the other functions to load in the tabulated refractive index values and the sea level refractive index value from the file. It also reads the mass overburden A,B and C values from the file
int MultiRayAirIceRefraction::MakeAtmosphere(){
   
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

  // flattened_h_data.clear();
  // flattened_nh_data.clear();
  
  return 0;
}

////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){
  
  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice
  // double AirTxHeight=5000;////Height of the source
  // double HorizontalDistance=1000;////Horizontal distance
  // double IceLayerHeight=3000;////Height where the ice layer starts off
  // double AntennaDepth=200;////Depth of antenna in the ice  

  //cout<<"Minimser values "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<endl;
  
  gsl_function F1;
  struct MinforLAng_params params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  F1.function = &MinimizeforLaunchAngle;
  F1.params = &params1;

  ////Set the initial angle limits for the minimisation
  double startanglelim=91.5;
  double endanglelim=178.5;

  ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
  bool checknan=false;
  double TotalHorizontalDistanceinAirt=0;
  int FilledLayerst=0;
  while(checknan==false){
    double *GetResultsAirTest1=GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
    TotalHorizontalDistanceinAirt=0;
    FilledLayerst=GetResultsAirTest1[4*MaxLayers];
    for(int i=0;i<FilledLayerst;i++){
      TotalHorizontalDistanceinAirt+=GetResultsAirTest1[i*4];
    }
    delete []GetResultsAirTest1;
    
    startanglelim=startanglelim-0.1;
    if(isnan(TotalHorizontalDistanceinAirt)==true){
      checknan=true;
      startanglelim=startanglelim+0.1*2;
      //cout<<"startangle is "<<startanglelim<<endl;
    }
  }

  checknan=false;
  while(checknan==false){
    double *GetResultsAirTest2=GetAirPropagationPar(endanglelim,AirTxHeight,IceLayerHeight);
    TotalHorizontalDistanceinAirt=0;
    FilledLayerst=GetResultsAirTest2[4*MaxLayers];
    for(int i=0;i<FilledLayerst;i++){
      TotalHorizontalDistanceinAirt+=GetResultsAirTest2[i*4];
    }
    delete []GetResultsAirTest2;
    
    endanglelim=endanglelim+0.1;
    if(isnan(TotalHorizontalDistanceinAirt)==true){
      checknan=true;
      endanglelim=endanglelim-0.1*2;
    }  
  }
  //std::cout<<"startangle "<<startanglelim<<" endangle "<<endanglelim<<std::endl;
  
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir=FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  
  //std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  double * GetResultsAir=GetAirPropagationPar(LaunchAngleAir,AirTxHeight,IceLayerHeight);
  int FilledLayers=GetResultsAir[4*MaxLayers];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  double LvalueAir;
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[i*4];
    PropagationTimeAir+=GetResultsAir[3+i*4];
  }
  LvalueAir=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*4];    
  delete [] GetResultsAir;
  
  // std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  // std::cout<<"***********Results for Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<std::endl;
  //   std::cout<<"LvalueAir "<<LvalueAir<<std::endl;
  // std::cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<std::endl;


  double * GetResultsIce=GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth,LvalueAir);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  // double IncidentAngleonAntenna=GetResultsIce[1];
  // double LvalueIce=GetResultsIce[2];
  double PropagationTimeIce=GetResultsIce[3];
  delete [] GetResultsIce;

  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<std::endl;
  // std::cout<<"LvalueIce "<<LvalueIce<<std::endl;
  // std::cout<<"PropagationTimeIce "<<PropagationTimeIce*pow(10,9)<<" ns"<<std::endl;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  //double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;

  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice + Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistance "<<TotalHorizontalDistance<<" m"<<std::endl;
  // std::cout<<"TotalPropagationTime "<<TotalPropagationTime*pow(10,9)<<" ns"<<std::endl;

  opticalPathLengthInIce=(PropagationTimeIce*MultiRayAirIceRefraction::spedc)*100;
  opticalPathLengthInAir=(PropagationTimeAir*MultiRayAirIceRefraction::spedc)*100;
  launchAngle=(LaunchAngleAir)*(pi/180);
  horizontalDistanceToIntersectionPoint=TotalHorizontalDistanceinAir*100;
  
  reflectionCoefficientS=Refl_S(IncidentAngleonIce*(pi/180),IceLayerHeight);
  reflectionCoefficientP=Refl_P(IncidentAngleonIce*(pi/180),IceLayerHeight);
  
  //cout<<" in minimser "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<endl;
  
  bool CheckSolution=false;
  double checkminimisation=TotalHorizontalDistance-HorizontalDistance;
  if(fabs(checkminimisation)<1){
    CheckSolution=true;
  }
  
  return CheckSolution;

}

////Just a simple function for interpolating in 1 dimension between two points (xa,ya) and (xb,yb)
double MultiRayAirIceRefraction::oneDLinearInterpolation(double x, double xa, double ya, double xb, double yb){
  double y=ya+(yb-ya)*((x-xa)/(xb-xa));
  return y;
}

////Just a simple function for extrapolating in 1 dimension outside two points
double MultiRayAirIceRefraction::Extrapolate(int Par, int index, double TotalHorizontalDistance, int AntennaNumber){
  double x1,x2,y1,y2,x3,m,c,y3;

  x1=AllTableAllAntData[AntennaNumber][1][index];
  x2=AllTableAllAntData[AntennaNumber][1][index+1];
  y1=AllTableAllAntData[AntennaNumber][Par][index];
  y2=AllTableAllAntData[AntennaNumber][Par][index+1];
  x3=TotalHorizontalDistance;
  
  m=(y2-y1)/(x2-x1);
  c=y1-m*x1;
  y3=m*x3+c;
  
  return y3;
}

////Find the limit to extrapolate for all the parameters using the Air Launch Angle parameter
double MultiRayAirIceRefraction::FindExtrapolationLimit( int index, double TotalHorizontalDistance, int AntennaNumber){
  double x1,x2,y1,y2,m,c,xlimit;
  //double y3,x3;
  
  x1=(AllTableAllAntData[AntennaNumber][1][index]);
  x2=(AllTableAllAntData[AntennaNumber][1][index+1]);
  y1=AllTableAllAntData[AntennaNumber][4][index];
  y2=AllTableAllAntData[AntennaNumber][4][index+1];
  //x3=(TotalHorizontalDistance);
  
  m=(y1-y2)/(x1-x2);
  c=y1-m*x1;
  //y3=m*x3+c;
  xlimit=(90-c)/m;
  
  return xlimit;
}

int MultiRayAirIceRefraction::FindClosestAirTxHeight(double IceLayerHeight, double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber){

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;  
  int CurrentHeightStep=floor((ParValue-LoopStopHeight)/HeightStepSize);

  int mainindex=(TotalHeightSteps-CurrentHeightStep-1);
  int index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
  index1=index1+mainindex;
  if(AllTableAllAntData[AntennaNumber][0][index1]<ParValue){
    mainindex=mainindex-1;
    index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
    index1=index1+mainindex;
  }
  int index2=((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))+index1;

  if(index2==TotalTableEntries-1){
   index1=index1-((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))-1;
   index2=index2-((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))-1;
  }
  
  if(ParValue==(AllTableAllAntData[AntennaNumber][0][index1]+AllTableAllAntData[AntennaNumber][0][index1-1])/2 && index1>0){
    mainindex=(TotalHeightSteps-CurrentHeightStep-1)-1;
    index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
    index1=index1+mainindex;
    index2=((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))+index1;
  }
  
  RStartIndex=index1;
  REndIndex=index2;
  ClosestVal=fabs(AllTableAllAntData[AntennaNumber][0][index1]-ParValue);

  return 0;
}

int MultiRayAirIceRefraction::FindClosestTHD(double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber){

  double minimum=100000000000;
  int index2=0;
  double minval=0;
  //double lastminval=0;
  for(int ipnt=StartIndex;ipnt<EndIndex+1;ipnt++){
    minval=fabs(AllTableAllAntData[AntennaNumber][1][ipnt]-ParValue);
    if(minval<minimum && AllTableAllAntData[AntennaNumber][1][ipnt]>ParValue){
      minimum=minval;
    }else{
      index2=ipnt; 
      break;
    }
    //lastminval=minval;
  }
  double index1=index2-1;

  minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index2]);
  if(minimum>fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1])){
    minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1]);
  }
  
  RStartIndex=index1;
  REndIndex=index2;
  ClosestVal=minimum;
  
  //cout<<"start index is "<<index1<<" last index is "<<index2<<" "<<minimum<<endl;
  return 0;
}

////Interpolate the value of the given parameter for a given TxHeight and THD
int MultiRayAirIceRefraction::GetParValues(double AntennaNumber, double AirTxHeight, double TotalHorizontalDistance, double IceLayerHeight, double &AirTxHeight1, double Par1[6],double &AirTxHeight2, double Par2[6]){

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
  MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
  MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];
  
  int startindex[6];
  int endindex[6];
  double closestvalue[6];
  double x1,x2,y1,y2;
  double MaxTotalHorizontalDistance=0;
  double MinTotalHorizontalDistance=0;
  
  ////Find out the AirTxheight1 bin range for the given Tx height
  FindClosestAirTxHeight(IceLayerHeight,AirTxHeight,0, TotalTableEntries, startindex[0], endindex[0], closestvalue[0], AntennaNumber);
  ////Set the first Tx height
  AirTxHeight1=AllTableAllAntData[AntennaNumber][0][startindex[0]];

  ////Find the maximum and minimum possible values of THD for AirTxHeight1. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
  MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[0]];
  MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[0]];
  //cout<<" minimum distance is 1 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[0]+1]<<endl;
  if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){
    ////Find out the THD bin range for the given THD and AirTxHeight1
    FindClosestTHD(TotalHorizontalDistance, startindex[0], endindex[0], startindex[1], endindex[1], closestvalue[1], AntennaNumber);
    
    if(closestvalue[1]!=0){////If the given THD does not match a bin value
      x1=AllTableAllAntData[AntennaNumber][1][startindex[1]];
      x2=AllTableAllAntData[AntennaNumber][1][endindex[1]];
      for(int ipar=0;ipar<6;ipar++){
	y1=AllTableAllAntData[AntennaNumber][2+ipar][startindex[1]];
	y2=AllTableAllAntData[AntennaNumber][2+ipar][endindex[1]];
	Par1[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
      }
    }
    if(closestvalue[1]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the same and there is no interpolation needed
      startindex[1]=startindex[1]+1;
      endindex[1]=startindex[1];
      for(int ipar=0;ipar<6;ipar++){
	Par1[ipar]=AllTableAllAntData[AntennaNumber][2+ipar][startindex[1]];
      }
    }
    
  }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
    ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
    double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);

    if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
      double THDLimit=FindExtrapolationLimit(startindex[0],TotalHorizontalDistance,AntennaNumber);
      double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
      double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
      if(DistancePercentage<DistancePercentageLimit){
	for(int ipar=0;ipar<6;ipar++){
	  Par1[ipar]=Extrapolate(2+ipar, startindex[0], TotalHorizontalDistance, AntennaNumber);
	}
      }else{
	for(int ipar=0;ipar<6;ipar++){
	  Par1[ipar]=-pow(10,9);
	}
      }
    }
    
  }
  
  ////Find out the Txheight2 bin range for the given Tx height
  startindex[2]=startindex[0]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)+1;
  endindex[2]=startindex[2]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize);
    
  if(closestvalue[0]!=0 && AirTxHeight>MinAirTxHeight && startindex[2]<TotalTableEntries){////If the AirTxHeight does not exactly match a table value and AirTxHeight1 is not equal to AirTxHeight2 AND AirTxHeight1 is not the minimum possible height in the table

    ////Set the second Tx height
    AirTxHeight2=AllTableAllAntData[AntennaNumber][0][startindex[2]];
    ////Find the maximum and minimum possible values of THD for AirTxHeight2. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
    MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[2]];
    MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[2]];
    //cout<<" minimum distance is 2 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[2]+1]<<endl;
    if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){

      FindClosestTHD(TotalHorizontalDistance, startindex[2], endindex[2], startindex[3], endindex[3], closestvalue[2], AntennaNumber);
       
      if(closestvalue[2]!=0){////If the given THD does not match a bin value
	x1=AllTableAllAntData[AntennaNumber][1][startindex[3]];
	x2=AllTableAllAntData[AntennaNumber][1][endindex[3]];
	for(int ipar=0;ipar<6;ipar++){
	  y1=AllTableAllAntData[AntennaNumber][2+ipar][startindex[3]];
	  y2=AllTableAllAntData[AntennaNumber][2+ipar][endindex[3]];
	  Par2[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
	}
      }
      if(closestvalue[2]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the sohame and there is no interpolation needed
	startindex[3]=startindex[3]+1;
	endindex[3]=startindex[3];
	for(int ipar=0;ipar<6;ipar++){
	  Par2[ipar]=AllTableAllAntData[AntennaNumber][2+ipar][startindex[3]];
	}
      }

    }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
      ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
      double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);
      if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
	double THDLimit=FindExtrapolationLimit(startindex[2],TotalHorizontalDistance,AntennaNumber);
	double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
	double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
	if(DistancePercentage<DistancePercentageLimit){
	  for(int ipar=0;ipar<6;ipar++){
	    Par2[ipar]=Extrapolate(2+ipar, startindex[2], TotalHorizontalDistance, AntennaNumber);
	  }
	}else{
	  for(int ipar=0;ipar<6;ipar++){
	    Par2[ipar]=-pow(10,9);
	  }
	}
      }
      
    }
    
  }else{////If AirTxHeight is exactly equal to a bin value then AirTxHeight1 and AirTxHeight2 are equal and we dont need to do interpolation for AirTxHeight2 as we already have done it for AirTxHeight1
    ////Also if AirTxHeight1 is equal to the minimum possible height then stop right there and dont the interpolation for AirTxHeight2
    AirTxHeight2=AirTxHeight1;
    for(int ipar=0;ipar<6;ipar++){
      Par2[ipar]=Par1[ipar];
    }
  }
  
  return 0; 
}

////This functions reads in the antenna tables and interpolates (or extrapolates) from the table to provide output value for raytracing
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint_Table(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, int AntennaNumber, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){

  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice

  //cout<<"Table values "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<endl;
  // double PropagationTimeIce=0;
  // double PropagationTimeAir=0;
  // double LaunchAngleAir=0;
  // double TotalHorizontalDistanceinAir=0;
  //double TotalHorizontalDistance=0;

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
  MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
  MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];
  
  double x1=0,x2=0,y1=0,y2=0;
  double AirTxHeight1,AirTxHeight2;
  double Par1[6];
  double Par2[6];
  double ParInterpolatedValues[6];
  ///0 is OpticalPathIce
  ///1 is OpticalPathAir
  ///2 is LaunchAngleAir
  ///3 is THDAir

  // cout<<"in here "<<AirTxHeight<<" "<<HorizontalDistance<<endl;
  if(AirTxHeight<=MaxAirTxHeight && AirTxHeight>=MinAirTxHeight && AirTxHeight>0){
    //cout<<"in here "<<endl;
   
    GetParValues(AntennaNumber,AirTxHeight,HorizontalDistance,IceLayerHeight,AirTxHeight1,Par1,AirTxHeight2,Par2);
    x1=AirTxHeight1;
    x2=AirTxHeight2;
    for(int ipar=0;ipar<6;ipar++){
      y1=Par1[ipar];
      y2=Par2[ipar];
      double GetInterpolatedParValue=0;
      int checkval=0;
      if(y1==-pow(10,9) || y2==-pow(10,9)){
	checkval=1;
      }
      if(x1!=x2 && checkval==0){
	GetInterpolatedParValue=oneDLinearInterpolation(AirTxHeight,x1,y1,x2,y2);
      }else{////If AirTxHeight is exactly equal to a bin value then AirTxHeight1 and AirTxHeight2 are equal and we dont need to do interpolation for AirTxHeight2 as we already have done it for AirTxHeight1
	if(x1==x2 && y1==y2){
	  GetInterpolatedParValue=Par1[ipar];
	}
	if(y2==-pow(10,9) && y1==-pow(10,9)){
	  ipar=5;
	}

      }
      ParInterpolatedValues[ipar]=GetInterpolatedParValue;
    }
    
  }
  
  opticalPathLengthInIce=ParInterpolatedValues[0]*100;
  opticalPathLengthInAir=ParInterpolatedValues[1]*100;
  launchAngle=ParInterpolatedValues[2]*(pi/180);
  horizontalDistanceToIntersectionPoint=ParInterpolatedValues[3]*100;
  reflectionCoefficientS=ParInterpolatedValues[4];
  reflectionCoefficientP=ParInterpolatedValues[5];

  bool CheckSolBool=false;
  if( (y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9)) ){
    CheckSolBool=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(SrcHeightASL, HorizontalDistanceToRx, RxDepthBelowIceBoundary, IceLayerHeight*100, opticalPathLengthInIce, opticalPathLengthInAir, launchAngle, horizontalDistanceToIntersectionPoint,reflectionCoefficientS,reflectionCoefficientP);
    //cout<<"in here 2 "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<y1<<" "<<y2<<endl;
  }
  
  bool CheckSolution=true;

  if(y2==-pow(10,9) && y1==-pow(10,9)){
    // std::cout<<"The given Total Horizontal Distance bigger than the maximum percentage of the Total Horizontal Distance range for the given AirTxHeight! Cannot extrapolate! "<<HorizontalDistance<<" "<<AirTxHeight<<" "<<IceLayerHeight<<std::endl;
    CheckSolution=false;
  }
  if( ((y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9))) && CheckSolBool==false ){
    // std::cout<<"The minimiser failed for the weird extrapolation case"<<std::endl;
    CheckSolution=false;    
  }
  if(AirTxHeight>MaxAirTxHeight){
//    std::cout<<"AirTx Height is greater than the maximum possible height in GDAS!"<<std::endl;
    CheckSolution=false;
  }
  if(AirTxHeight<MinAirTxHeight){
//    std::cout<<"AirTx Height is less than the minimum possible height in the table!"<<std::endl;
    CheckSolution=false;
  }
  if(AirTxHeight<0){
//    std::cout<<"AirTx Height is less than the zero!"<<std::endl;
    CheckSolution=false;
  }

  if(CheckSolution==false){
    opticalPathLengthInIce=0;
    opticalPathLengthInAir=0;
    launchAngle=0;
    horizontalDistanceToIntersectionPoint=0;
  }
  
  return CheckSolution;
}


int MultiRayAirIceRefraction::MakeRayTracingTable(double AntennaDepth, double IceLayerHeight, int AntennaNumber){

  std::vector<std::vector <double>> AllTableData;
  
  ////convert cm to m
  AntennaDepth=AntennaDepth/100;
  IceLayerHeight=IceLayerHeight/100;
  
  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();
  
  ////Print out the entry number, the Tx height, ice layer height, Tx height above the icelayer height, total horizontal distance on surface, total horizontal distance in ice, RayLaunchAngle at Tx, incident angle on ice and recievd angle in ice at the antenna inside this file
  MakeAtmosphere();
 
  ////Define variables for the loop over Tx height and ray launch angle
  double RayLaunchAngleInAir=0;////Set zero for now and 0 deg straight up. This variable defines the initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 deg is straight up
  //double AirTxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];////Maximum height available with the refractive index data
  double AirTxHeight=100000;///Maximum height available with the refractive index data
  
  // ////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
  // double AngleStepSize=0.5;
  // double LoopStartAngle=92;
  // double LoopStopAngle=178;
  // int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize);

  ////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
  LoopStartHeight=AirTxHeight;
  LoopStopHeight=IceLayerHeight;
  TotalHeightSteps=floor((LoopStartHeight-LoopStopHeight)/HeightStepSize)+1;

  ////Temporary vectors
  std::vector <double> temp1;
  std::vector <double> temp2;
  std::vector <double> temp3;
  std::vector <double> temp4;
  std::vector <double> temp5;
  std::vector <double> temp6;
  std::vector <double> temp7;
  std::vector <double> temp8;
  
  int ifileentry=0;
  //cout<<" our final angle is "<<LoopStartAngle+AngleStepSize*TotalAngleSteps<<endl;
  ////Start looping over the Tx Height and Launch angle values
  for(int ihei=0;ihei<TotalHeightSteps;ihei++){
    AirTxHeight=LoopStartHeight-HeightStepSize*ihei;
    if(AirTxHeight!=IceLayerHeight && ihei==TotalHeightSteps-1){
      AirTxHeight=IceLayerHeight;
    }
    if(AirTxHeight>0){
      for(int iang=0;iang<TotalAngleSteps;iang++){
	RayLaunchAngleInAir=LoopStartAngle+AngleStepSize*iang;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Section for propogating the ray through the atmosphere
	
	////Find out how many atmosphere layers are above the source or Tx which we do not need
	int skiplayer=0;
	for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
	  //std::cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<std::endl;
	  if(AirTxHeight<ATMLAY[ilayer]/100 && AirTxHeight>=ATMLAY[ilayer-1]/100){
	    //std::cout<<"Tx Height is in this layer with a height range of "<<MultiRayAirIceRefraction::ATMLAY[ilayer]/100<<" m to "<<MultiRayAirIceRefraction::ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<std::endl;
	    ilayer=-100;
	  }
	  if(ilayer>-1){
	    skiplayer++;
	  }
	}
	int SkipLayersAbove=skiplayer;
	//std::cout<<"The tota number of layers that need to be skipped from above is "<<skiplayer<<std::endl;
      
	////Find out how many atmosphere layers are below the ice height which we do not need
	skiplayer=0;
	for(int ilayer=0;ilayer<MaxLayers;ilayer++){
	  //std::cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer+1]/100<<std::endl;
	  if(IceLayerHeight>=ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
	    //std::cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<std::endl;
	    ilayer=100;
	  }
	  if(ilayer<MaxLayers){
	    skiplayer++;
	  }
	}
	int SkipLayersBelow=skiplayer;
	//std::cout<<"The total number of layers that need to be skipped from below is "<<skiplayer<<std::endl;
      
	////Define variables for ray propogation through mutliple layers in the atmosphere
	double Start_nh=0;
	double StartHeight=0;
	double StopHeight=0;
	double StartAngle=0;
	double TotalHorizontalDistanceInAir=0;
	double TimeInAir=0;
      
	////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
	//std::cout<<"Fitting the atmosphere refrative index profile with multiple layers and propogate the ray"<<std::endl;
	for(int ilayer=MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
	  //std::cout<<B_air<<std::endl;
	  ////Set the starting height of the ray for propogation for that layer
	  if(ilayer==MaxLayers-SkipLayersAbove-1){
	    ////If this is the first layer then set the start height to be the height of the source
	    StartHeight=AirTxHeight;
	  }else{
	    ////If this is any layer after the first layer then set the start height to be the starting height of the layer
	    StartHeight=ATMLAY[ilayer+1]/100-0.00001;
	  }
	
	  ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
	  Start_nh=Getnz_air(StartHeight);
	
	  ////Set the staopping height of the ray for propogation for that layer
	  if(ilayer==(SkipLayersBelow-1)+1){
	    ////If this is the last layer then set the stopping height to be the height of the ice layer
	    StopHeight=IceLayerHeight;
	  }else{
	    ////If this is NOT the last layer then set the stopping height to be the end height of the layer
	    StopHeight=ATMLAY[ilayer]/100;
	  }
	
	  ////If this is the first layer then set the initial launch angle of the ray through the layers. I calculate the final launch angle by doing 180-RayLaunchAngleInAir since my raytracer only works with 0 to 90 deg. Setting an angle of 95 deg w.r.t to the vertical where 0 is up means that my raytraces takes in an launch angle of 85.
	  if(ilayer==MaxLayers-SkipLayersAbove-1){
	    StartAngle=180-RayLaunchAngleInAir;
	  }
	  //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<std::endl;
	
	  ////Get the hit parameters from the function. The output is:
	  //// How much horizontal distance did the ray travel in the layer
	  //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
	  //// The value of the L parameter for that layer
	  double* GetHitPar=GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
	  TotalHorizontalDistanceInAir+=GetHitPar[0];
	  StartAngle=GetHitPar[1];
	  TimeInAir+=GetHitPar[3];
	
	  ////dont forget to delete the pointer!
	  delete []GetHitPar;
	}
      
	double IncidentAngleonIce=StartAngle;
	//std::cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<std::endl;
      
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Section for propogating the ray through the ice
            
	////Set the starting depth of the ray for propogation to at the ice surface
	double StartDepth=0.0;
	////Since we have the starting height of the ice layer we can find out the refactive index of air at that height from data using spline interpolation
	//double Start_nh=gsl_spline_eval(spline, IceLayerHeight, accelerator);
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
	double TotalHorizontalDistanceInIce=GetHitPar[0];
	double RecievdAngleInIce=GetHitPar[1];
	//double LvalueIce=GetHitPar[2];
	double TimeInIce=GetHitPar[3];
      
	//std::cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceInIce<<std::endl;
      
	//if(isnan(TotalHorizontalDistanceInAir)==false){
	
	  ////define dummy/temporary variables for storing data
	double dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7,dummy8,dummy9,dummy10,dummy11,dummy12,dummy13,dummy14,dummy15,dummy16;  
	  
	  dummy1=ifileentry;
	  dummy2=AirTxHeight;
	  dummy3=TotalHorizontalDistanceInAir + TotalHorizontalDistanceInIce;
	  dummy4=TotalHorizontalDistanceInAir;
	  dummy5=TotalHorizontalDistanceInIce;
	  dummy6=(TimeInIce+TimeInAir)*MultiRayAirIceRefraction::spedc;
	  dummy7=TimeInAir*MultiRayAirIceRefraction::spedc;
	  dummy8=TimeInIce*MultiRayAirIceRefraction::spedc;
	  dummy9=(TimeInIce+TimeInAir)*pow(10,9);
	  dummy10=TimeInAir*pow(10,9);
	  dummy11=TimeInIce*pow(10,9);
	  dummy12=RayLaunchAngleInAir;
	  dummy13=IncidentAngleonIce;
	  dummy14=RecievdAngleInIce;
	  dummy15=MultiRayAirIceRefraction::Refl_S(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);
	  dummy16=MultiRayAirIceRefraction::Refl_P(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);

	  //	  cout<<dummy1<<" "<<dummy2<<" "<<dummy3<<" "<<dummy4<<" "<<dummy5<<" "<<dummy6<<" "<<dummy7<<" "<<dummy8<<" "<<dummy9<<" "<<dummy10<<" "<<dummy11<<" "<<dummy12<<" "<<dummy13<<" "<<dummy14<<endl;;  
	  
	  temp1.push_back(dummy2);///AirTx Height
	  temp2.push_back(dummy3);///THDTotal
	  temp3.push_back(dummy8);///OpticalPathIce
	  temp4.push_back(dummy7);///OpticalPathAir
	  temp5.push_back(dummy12);///LaunchAngleAir
	  temp6.push_back(dummy4);///THDAir
	  temp7.push_back(dummy15);///THDAir
	  temp8.push_back(dummy16);///THDAir

	  ifileentry++;
	  //}
	delete[] GetHitPar;
     
      }//// end of iang loop
    }////if condition to make sure AirTxHeight does not go below zero
  }//// end of ihei loop

  
  AllTableData.push_back(temp1);
  AllTableData.push_back(temp2);
  AllTableData.push_back(temp3);
  AllTableData.push_back(temp4);
  AllTableData.push_back(temp5);
  AllTableData.push_back(temp6);
  AllTableData.push_back(temp7);
  AllTableData.push_back(temp8);

  AllTableAllAntData.push_back(AllTableData);
  
  AllTableData.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();
  temp4.clear();
  temp5.clear();
  temp6.clear();
  temp7.clear();
  temp8.clear();
  
  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  durationb=durationb/1000000;
  std::cout<<"total time taken by the script to generate the table: "<<durationb<<" s"<<std::endl;
  return 0;
  
}
