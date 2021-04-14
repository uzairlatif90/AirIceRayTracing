#include "RayTracingFunctions.h"

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
int RayTracingFunctions::readATMpar(){
  
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

  abc[4][0]=abc[3][0];
  abc[4][1]=abc[3][1];
  abc[4][2]=abc[3][2];
  
  return 0;
}

int RayTracingFunctions::readnhFromFile(){

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
    ////define dummy/temporary vectors for storing data.
    std::vector <double> temp1,temp2,temp3;
    
    while (getline(ain,line)){
      ain>>dummy1>>dummy2;
      
      if(dummy1>-1){////start storing height at above and equal to 0 m
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

  MaxLayers=h_data.size()+1;////store the total number of layers present in the data
  
  return 0;
}

////Get the value of the B parameter for the ice refractive index model
double RayTracingFunctions::GetB_ice(double z){
  double zabs=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

////Get the value of the C parameter for the ice refractive index model
double RayTracingFunctions::GetC_ice(double z){
  double zabs=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth in ice
double RayTracingFunctions::Getnz_ice(double z){
  z=fabs(z);
  return RayTracingFunctions::A_ice+RayTracingFunctions::GetB_ice(z)*exp(-RayTracingFunctions::GetC_ice(z)*z);
}

int RayTracingFunctions::FillInAirRefractiveIndex(){
  
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
  }

  return 0;   
}

////Get the value of the B parameter for the air refractive index model
double RayTracingFunctions::GetB_air(double z){
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
  //B=1e-9;
  return B;
}

////Get the value of the C parameter for the air refractive index model
double RayTracingFunctions::GetC_air(double z){
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
  //C=1e-9;
  return C;
}

////Get the value of refractive index model for a given depth in air
double RayTracingFunctions::Getnz_air(double z){
  double zabs=fabs(z);

  return RayTracingFunctions::A_air+RayTracingFunctions::GetB_air(zabs)*exp(-RayTracingFunctions::GetC_air(zabs)*zabs);
}

////E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R
double RayTracingFunctions::Refl_S(double thetai, double IceLayerHeight){
  double Nair=RayTracingFunctions::Getnz_air(IceLayerHeight);
  double Nice=RayTracingFunctions::Getnz_ice(0); 
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
double RayTracingFunctions::Refl_P(double thetai, double IceLayerHeight){
  double Nair=RayTracingFunctions::Getnz_air(IceLayerHeight);
  double Nice=RayTracingFunctions::Getnz_ice(0); 
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
double RayTracingFunctions::FindFunctionRoot(gsl_function F,double x_lo, double x_hi,const gsl_root_fsolver_type *T,double tolerance)
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

////Analytical solution describing the ray path in ice
double RayTracingFunctions::fDnfR(double x,void *params){
  
  struct RayTracingFunctions::fDnfR_params *p= (struct RayTracingFunctions::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
double RayTracingFunctions::fdxdz(double x,void *params){
  
  struct RayTracingFunctions::fdxdz_params *p= (struct RayTracingFunctions::fdxdz_params *) params;
  double Lang = p->lang;
  double Z0 = p->z0;
  double Z1 = p->z1;
  int AirOrIce = p->airorice;

  double output=0;
  if(AirOrIce==0){
    output=tan(asin( (RayTracingFunctions::Getnz_ice(Z0)*sin(x))/RayTracingFunctions::Getnz_ice(Z1) ) ) - tan(Lang);
  }
  
  if(AirOrIce==1){
    output=tan(asin( (RayTracingFunctions::Getnz_air(Z0)*sin(x))/RayTracingFunctions::Getnz_air(Z1) ) ) - tan(Lang);
  }
  return output;
}

////The function used to calculate ray propogation time in ice
double RayTracingFunctions::ftimeD(double x,void *params){

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

double RayTracingFunctions::GetRayOpticalPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct RayTracingFunctions::fDnfR_params params2a;
  struct RayTracingFunctions::fDnfR_params params2b;
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

double RayTracingFunctions::GetRayPropagationTime(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct RayTracingFunctions::ftimeD_params params3a;
  struct RayTracingFunctions::ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), RayTracingFunctions::spedc, Lvalue,0};
    params3b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), RayTracingFunctions::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, GetB_air(RxDepth), -GetC_air(RxDepth), RayTracingFunctions::spedc, Lvalue,1};
    params3b = {A, GetB_air(TxDepth), -GetC_air(TxDepth), RayTracingFunctions::spedc, Lvalue,1};
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
double *RayTracingFunctions::GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){

  double *output=new double[4];
  
  //double x0=0;////Starting horizontal point of the ray. Always set at zero
  double x1=0;////Variable to store the horizontal distance that will be traveled by the ray
  
  double ReceiveAngle=0;////Angle from the vertical at which the target will recieve the ray
  double Lvalue=0;//// L parameter of the ray for that layer
  double RayTimeIn2ndLayer=0;////Time of propagation in 2ndLayer 
  //double AngleOfEntryIn2ndLayer=0;////Angle at which the ray enters the layer

  double SurfaceRayIncidentAngle=IncidentAng*(RayTracingFunctions::pi/180.0);////Angle at which the ray is incident on the second layer
  double RayAngleInside2ndLayer=0;////Use Snell's Law to find the angle of transmission in the 2ndlayer

  double A=0;
  double nzRx=0;
  double nzTx=0;
  double GSLFnLimit=0;

  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    A=RayTracingFunctions::A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    A=RayTracingFunctions::A_air;
    nzRx=Getnz_air(RxDepth);
    nzTx=Getnz_air(TxDepth);
  }

  ////LimitAngle sets a limit on the range to which the GSL minimisation will work. This limit comes from the fact that in fdxdx() you have tan(asin(x)) which goes to infinity at x=1. In our case x=(nz(Z0)*sin(Angle))/nz(Z1) . Solving for Angle gives us our limit.
  double LimitAngle=asin(nzTx/nzRx);
  
  GSLFnLimit=LimitAngle;
  RayAngleInside2ndLayer=asin((n_layer1/nzTx)*sin(SurfaceRayIncidentAngle));////Use Snell's Law to find the angle of transmission in the 2ndlayer
  
  ////calculate the angle at which the target receives the ray
  gsl_function F1;
  struct RayTracingFunctions::fdxdz_params params1 = {RayAngleInside2ndLayer, RxDepth, TxDepth, AirOrIce};
  F1.function = &fdxdz;
  F1.params = &params1;
  ReceiveAngle=RayTracingFunctions::FindFunctionRoot(F1,0*(RayTracingFunctions::pi/180),GSLFnLimit, gsl_root_fsolver_falsepos,0.000000001);
  //std::cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/RayTracingFunctions::pi)<<" deg"<<std::endl;
  
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
  // F4.function = &RayTracingFunctions::fDnfR;
  // F4.params = &params2b;
  // gsl_deriv_central (&F4, TxDepth, 1e-8, &result, &abserr);
  // AngleOfEntryIn2ndLayer=atan(result)*(180.0/RayTracingFunctions::pi);
  // if(TxDepth==RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=180-ReceiveAngle;
  // }
  // if(TxDepth!=RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=90;
  // }
  //std::cout<<"AngleOfEntryIn2ndLayer= "<<AngleOfEntryIn2ndLayer<<" ,RayAngleInside2ndLayer="<<RayAngleInside2ndLayer*(180/RayTracingFunctions::pi)<<std::endl;

  output[0]=x1;
  output[1]=ReceiveAngle*(180/RayTracingFunctions::pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  
  return output;
}

////This function flattens out 2d std::vectors into 1d std::vectors
std::vector<double> RayTracingFunctions::flatten(const std::vector<std::vector<double>>& v) {
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
double * RayTracingFunctions::GetAirPropagationPar(double LaunchAngleAir, double AirTxHeight, double IceLayerHeight){
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
      StartAngle=180-LaunchAngleAir;
    }
    //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;
    
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      double* GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
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
      RecAng=RecAng*(180/RayTracingFunctions::pi);
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
double * RayTracingFunctions::GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth, double Lvalue){
  double *output=new double[4];

  double StartAngle=IncidentAngleonIce;
  double StartDepth=0.0;
  double StopDepth=AntennaDepth;
  double nzStopDepth=Getnz_ice(StopDepth);
  
  double TotalHorizontalDistance=GetRayOpticalPath(A_ice, StopDepth, StartDepth, Lvalue, 0);
  double ReceiveAngle=asin(Lvalue/nzStopDepth)*(180/RayTracingFunctions::pi);
  double PropagationTime=GetRayPropagationTime(A_ice, StopDepth, StartDepth, Lvalue, 0);

  output[0]=TotalHorizontalDistance;
  output[1]=ReceiveAngle;
  output[2]=Lvalue;
  output[3]=PropagationTime;

  return output;
}


////This function will be minimised to get the launch angle in air. It looks at the total horizontal distance and tries to find a ray in air and ice which can cover that horizontal distance perfectly. It takes into account the discontinuity in refractive index at the air ice boundary
double RayTracingFunctions::MinimizeforLaunchAngle(double x, void *params){
  
  struct RayTracingFunctions::MinforLAng_params *p= (struct RayTracingFunctions::MinforLAng_params *) params;
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
int RayTracingFunctions::MakeAtmosphere(){
   
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
