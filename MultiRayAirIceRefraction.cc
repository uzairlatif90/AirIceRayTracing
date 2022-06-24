#include "MultiRayAirIceRefraction.h"

////Store the maximum possible height allowed by GDAS tables
double MaxAirTxHeight=0;
////Store the minimum possible height allowed antenna tables
double MinAirTxHeight=0;

////Multidimensional vector to store data read in from the antenna tables for interpolation puprposes
std::vector<std::vector<std::vector <double>>> AllTableAllAntData;

////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
double AngleStepSize=0.1;
double LoopStartAngle=90.1;
double LoopStopAngle=180.0;
int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize)+1;

////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
double HeightStepSize=10;
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
	//cout<<n1<<" "<<dummya[0]<<" , "<<dummya[1]<<" , "<<dummya[2]<<" , "<<dummya[3]<<" , "<<dummya[4]<<endl;
      }

      ////Store the values in their respective arrays
      if(n1==0){
	for (int i=0; i<5; i++){ MultiRayAirIceRefraction::ATMLAY[i]=dummya[i]; }
      }    
      if(n1==1){
	for (int i=0; i<5; i++){ MultiRayAirIceRefraction::abc[i][0]=dummya[i]; }
      }
      if(n1==2){
	for (int i=0; i<5; i++){ MultiRayAirIceRefraction::abc[i][1]=dummya[i]; }
      }
      if(n1==3){
	for (int i=0; i<5; i++){ MultiRayAirIceRefraction::abc[i][2]=dummya[i]; }
      }
      n1++;
    }////end the while loop
    
    ain.close();
  }////if condition to check if file is open

  MultiRayAirIceRefraction::abc[4][0]=MultiRayAirIceRefraction::abc[3][0];
  MultiRayAirIceRefraction::abc[4][1]=MultiRayAirIceRefraction::abc[3][1];
  MultiRayAirIceRefraction::abc[4][2]=MultiRayAirIceRefraction::abc[3][2];

  MultiRayAirIceRefraction::ATMLAY[4]=150000*100;////set max possible height that can be used in cm

  //std::cout<<"abc values are "<<MultiRayAirIceRefraction::abc[4][0]<<" "<<MultiRayAirIceRefraction::abc[4][1]<<" "<<MultiRayAirIceRefraction::abc[4][2]<<std::endl;
  
  return 0;
}

int MultiRayAirIceRefraction::readnhFromFile(){

  MultiRayAirIceRefraction::nh_data.clear();
  MultiRayAirIceRefraction::lognh_data.clear();
  MultiRayAirIceRefraction::h_data.clear();
  
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
	
	if(dummy1*100>=MultiRayAirIceRefraction::ATMLAY[layer]){////change the layer once the data of all the heights of that layer has been read in
	  if(layer>0){////now since the layer has finished and the temporary vectors have been filled in. Now we push the vectors in the main 2d height and refractice index vectors
	    MultiRayAirIceRefraction::h_data.push_back(temp1);
	    MultiRayAirIceRefraction::nh_data.push_back(temp2);
	    MultiRayAirIceRefraction::lognh_data.push_back(temp3);

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
      MultiRayAirIceRefraction::h_data.push_back(temp1);
      MultiRayAirIceRefraction::nh_data.push_back(temp2);
      MultiRayAirIceRefraction::lognh_data.push_back(temp3);
      ////clear the vectors now for storing the next layer
      temp1.clear();
      temp2.clear();
      temp3.clear();
    }
    layer++;
    
    ain.close();
  }////if condition to check if file is open

  ////The file reading condition "while (getline(ain,line))" reads the last the datapoint of the file twice. This is to to remove the last repeat data point in all the data arrays
  MultiRayAirIceRefraction::h_data[MultiRayAirIceRefraction::h_data.size()-1].erase(MultiRayAirIceRefraction::h_data[MultiRayAirIceRefraction::h_data.size()-1].end() - 1);
  MultiRayAirIceRefraction::nh_data[MultiRayAirIceRefraction::nh_data.size()-1].erase(MultiRayAirIceRefraction::nh_data[MultiRayAirIceRefraction::nh_data.size()-1].end() - 1);
  MultiRayAirIceRefraction::lognh_data[MultiRayAirIceRefraction::lognh_data.size()-1].erase(MultiRayAirIceRefraction::lognh_data[MultiRayAirIceRefraction::lognh_data.size()-1].end() - 1);

  MultiRayAirIceRefraction::MaxLayers=MultiRayAirIceRefraction::h_data.size()+1;////store the total number of layers present in the data

  //cout<<"max layers are "<<MaxLayers<<endl;
  
  return 0;
}

////Get the value of the B parameter for the ice refractive index model
double MultiRayAirIceRefraction::GetB_ice(double z){
  double zabs=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

////Get the value of the C parameter for the ice refractive index model
double MultiRayAirIceRefraction::GetC_ice(double z){
  double zabs=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth in ice
double MultiRayAirIceRefraction::Getnz_ice(double z){
  z=fabs(z);
  return MultiRayAirIceRefraction::A_ice+MultiRayAirIceRefraction::GetB_ice(z)*exp(-MultiRayAirIceRefraction::GetC_ice(z)*z);
}

int MultiRayAirIceRefraction::FillInAirRefractiveIndex(){
  
  double N0=0;
  for(int ilayer=0;ilayer<5;ilayer++){
    double hlow=MultiRayAirIceRefraction::ATMLAY[ilayer]/100;
    MultiRayAirIceRefraction::C_air[ilayer]=1.0/(MultiRayAirIceRefraction::abc[ilayer][2]/100);
    if(ilayer>0){
      N0=MultiRayAirIceRefraction::A_air+MultiRayAirIceRefraction::B_air[ilayer-1]*exp(-hlow*MultiRayAirIceRefraction::C_air[ilayer-1]);
    }
    if(ilayer==0){
      N0=gsl_spline_eval(spline, 0, accelerator);
    }
    MultiRayAirIceRefraction::B_air[ilayer]=((N0-1)/exp(-hlow*MultiRayAirIceRefraction::C_air[ilayer]));
  }

  // for(int ilayer=0;ilayer<5;ilayer++){
  //   cout<<A_air<<" "<<B_air[ilayer]<<" "<<C_air[ilayer]<<endl;
  // }

  return 0;   
}

////Get the value of the B parameter for the air refractive index model
double MultiRayAirIceRefraction::GetB_air(double z){
  double zabs=fabs(z);
  double B=0;
  int whichlayer=0;
 
  for(int ilayer=0;ilayer<MultiRayAirIceRefraction::MaxLayers-1;ilayer++){

    if(zabs<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100 && zabs>=MultiRayAirIceRefraction::ATMLAY[ilayer]/100){
      whichlayer=ilayer;
      ilayer=100;
    }  
  }
  if(zabs>=MultiRayAirIceRefraction::ATMLAY[MultiRayAirIceRefraction::MaxLayers-1]/100){
    whichlayer=MultiRayAirIceRefraction::MaxLayers-1;
  }
 
  B=MultiRayAirIceRefraction::B_air[whichlayer];
  //B=1e-9;
  return B;
}

////Get the value of the C parameter for the air refractive index model
double MultiRayAirIceRefraction::GetC_air(double z){
  double zabs=fabs(z);
  double C=0;
  int whichlayer=0;
  
  for(int ilayer=0;ilayer<MultiRayAirIceRefraction::MaxLayers-1;ilayer++){
    if(zabs<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100 && zabs>=MultiRayAirIceRefraction::ATMLAY[ilayer]/100){
      whichlayer=ilayer;
      ilayer=100;
    }
  }
  
  if(zabs>=MultiRayAirIceRefraction::ATMLAY[MultiRayAirIceRefraction::MaxLayers-1]/100){
    whichlayer=MultiRayAirIceRefraction::MaxLayers-1;
  }
  C=MultiRayAirIceRefraction::C_air[whichlayer];
  //C=1e-9;
  return C;
}

////Get the value of refractive index model for a given depth in air
double MultiRayAirIceRefraction::Getnz_air(double z){
  double zabs=fabs(z);

  return MultiRayAirIceRefraction::A_air+MultiRayAirIceRefraction::GetB_air(zabs)*exp(-MultiRayAirIceRefraction::GetC_air(zabs)*zabs);
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
  int iter = 0, max_iter = 20;
  //const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  //double tolerance=0.000000001;
  
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  //printf ("using %s method\n", gsl_root_fsolver_name (s));
  //printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  //cout<<" we are here "<<endl;
  //cout<<x_lo<<" "<<x_hi<<" "<<r<<" "<<iter<<endl; 
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

  //cout<<"we are here now "<<endl;
  
  gsl_root_fsolver_free (s);

  return r;
}

////Analytical solution describing the ray path in ice
double MultiRayAirIceRefraction::fDnfR(double x,void *params){
  
  struct MultiRayAirIceRefraction::fDnfR_params *p= (struct MultiRayAirIceRefraction::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

// ////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
// double MultiRayAirIceRefraction::fdxdz(double x,void *params){
  
//   struct MultiRayAirIceRefraction::fdxdz_params *p= (struct MultiRayAirIceRefraction::fdxdz_params *) params;
//   double Lang = p->lang;
//   double Z0 = p->z0;
//   double Z1 = p->z1;
//   int AirOrIce = p->airorice;

//   double output=0,dumx=0;
//   if(AirOrIce==0){
//     dumx=(MultiRayAirIceRefraction::Getnz_ice(Z0)*sin(x))/MultiRayAirIceRefraction::Getnz_ice(Z1);
//   }
//   if(AirOrIce==1){
//     dumx=(MultiRayAirIceRefraction::Getnz_air(Z0)*sin(x))/MultiRayAirIceRefraction::Getnz_air(Z1);
//   }
//   //output=((dumx/sqrt(1-dumx*dumx)) - tan(Lang));
//   //cout<<"output is "<<output<<" "<<x<<endl;
//   output=dumx - sin(Lang);
  
//   return output;
// }

////The function used to calculate ray propogation time in ice
double MultiRayAirIceRefraction::ftimeD(double x,void *params){

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

/* The function is used to calculate ray geometric path in ice */
double MultiRayAirIceRefraction::fpathD(double x,void *params){

  struct MultiRayAirIceRefraction::ftimeD_params *p= (struct MultiRayAirIceRefraction::ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;

  //integral sec(sin^(-1)(L/(A + B e^(C x)))) dx = (log((A + B e^(C x)) (sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + 1)) - (A log(A sqrt(A^2 - L^2) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + B sqrt(A^2 - L^2) e^(C x) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + A^2 + A B e^(C x) - L^2))/sqrt(A^2 - L^2) + (A C x)/sqrt(A^2 - L^2))/C;
  
  return (log((A + B*exp(C*x))*(sqrt((A*A + 2*A*B*exp(C*x) + B*B*exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x))) ) + 1)) - (A*log(A*sqrt(A*A - L*L)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/(( A + B*exp(C*x))*(A + B*exp(C*x)))) + B*sqrt(A*A - L*L)*exp(C*x)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x)))) + A*A + A*B*exp(C*x) - L*L))/sqrt(A*A - L*L) + (A*C*x)/sqrt(A*A - L*L))/C ;

}

double MultiRayAirIceRefraction::GetRayHorizontalPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct MultiRayAirIceRefraction::fDnfR_params params2a;
  struct MultiRayAirIceRefraction::fDnfR_params params2b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params2a = {A, MultiRayAirIceRefraction::GetB_ice(RxDepth), -MultiRayAirIceRefraction::GetC_ice(RxDepth), Lvalue};
    params2b = {A, MultiRayAirIceRefraction::GetB_ice(TxDepth), -MultiRayAirIceRefraction::GetC_ice(TxDepth), Lvalue};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params2a = {A, MultiRayAirIceRefraction::GetB_air(RxDepth), -MultiRayAirIceRefraction::GetC_air(RxDepth), Lvalue};
    params2b = {A, MultiRayAirIceRefraction::GetB_air(TxDepth), -MultiRayAirIceRefraction::GetC_air(TxDepth), Lvalue};
  }
  double x1=+MultiRayAirIceRefraction::fDnfR(RxDepth,&params2a)-MultiRayAirIceRefraction::fDnfR(TxDepth,&params2b);
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
    params3a = {A, MultiRayAirIceRefraction::GetB_ice(RxDepth), -MultiRayAirIceRefraction::GetC_ice(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
    params3b = {A, MultiRayAirIceRefraction::GetB_ice(TxDepth), -MultiRayAirIceRefraction::GetC_ice(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, MultiRayAirIceRefraction::GetB_air(RxDepth), -MultiRayAirIceRefraction::GetC_air(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
    params3b = {A, MultiRayAirIceRefraction::GetB_air(TxDepth), -MultiRayAirIceRefraction::GetC_air(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
  }
  double RayTimeIn2ndLayer=+MultiRayAirIceRefraction::ftimeD(RxDepth,&params3a)-MultiRayAirIceRefraction::ftimeD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  
  return RayTimeIn2ndLayer;
}

double MultiRayAirIceRefraction::GetRayGeometricPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct MultiRayAirIceRefraction::ftimeD_params params3a;
  struct MultiRayAirIceRefraction::ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, MultiRayAirIceRefraction::GetB_ice(RxDepth), -MultiRayAirIceRefraction::GetC_ice(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
    params3b = {A, MultiRayAirIceRefraction::GetB_ice(TxDepth), -MultiRayAirIceRefraction::GetC_ice(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, MultiRayAirIceRefraction::GetB_air(RxDepth), -MultiRayAirIceRefraction::GetC_air(RxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
    params3b = {A, MultiRayAirIceRefraction::GetB_air(TxDepth), -MultiRayAirIceRefraction::GetC_air(TxDepth), MultiRayAirIceRefraction::spedc, Lvalue,1};
  }
  double RayGeometricPath=+MultiRayAirIceRefraction::fpathD(RxDepth,&params3a)-MultiRayAirIceRefraction::fpathD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayGeometricPath*=-1;
  }
  
  return RayGeometricPath;
}

////Get the distance on the 2ndLayer at which the ray should hit given an incident angle such that it hits an target at depth of z0 m in the second layer.
//// n_layer1 is the refractive index value of the previous layer at the boundary of the two mediums
//// A,B and C values are the values required for n(z)=A+B*exp(Cz) for the second layer
//// TxDepth is the starting height or depth
//// RxDepth is the final height or depth
//// AirOrIce variable is used to determine whether we are working in air or ice as that sets the range for the GSL root finder.
double *MultiRayAirIceRefraction::GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){

  //std::cout<<"in new function "<<n_layer1<<" "<<RxDepth<<" "<<TxDepth<<" "<<IncidentAng<<" "<<AirOrIce<<std::endl;
  
  double *output=new double[5];

  //auto t1a = std::chrono::high_resolution_clock::now();
  //double x0=0;////Starting horizontal point of the ray. Always set at zero
  double x1=0;////Variable to store the horizontal distance that will be traveled by the ray
  double x1_Geo=0;
  
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

  // auto t2a = std::chrono::high_resolution_clock::now();
  // auto t1b = std::chrono::high_resolution_clock::now();
  
  ////LimitAngle sets a limit on the range to which the GSL minimisation will work. This limit comes from the fact that in fdxdx() you have tan(asin(x)) which goes to infinity at x=1. In our case x=(nz(Z0)*sin(Angle))/nz(Z1) . Solving for Angle gives us our limit.
  double LimitAngle=asin(nzTx/nzRx);
  
  GSLFnLimit=LimitAngle;
  RayAngleInside2ndLayer=asin((n_layer1/nzTx)*sin(SurfaceRayIncidentAngle));////Use Snell's Law to find the angle of transmission in the 2ndlayer
  
  ////calculate the angle at which the target receives the ray
  // gsl_function F1;
  // struct MultiRayAirIceRefraction::fdxdz_params params1 = {RayAngleInside2ndLayer, RxDepth, TxDepth, AirOrIce};
  // F1.function = &fdxdz;
  // F1.params = &params1;
  // //cout<<"limits are "<<RayAngleInside2ndLayer*(MultiRayAirIceRefraction::pi/180)<<" "<<GSLFnLimit*(180.0/MultiRayAirIceRefraction::pi)<<endl;
  // ReceiveAngle=MultiRayAirIceRefraction::FindFunctionRoot(F1,0.0*(MultiRayAirIceRefraction::pi/180),GSLFnLimit, gsl_root_fsolver_brent,0.00000001);

  double Lang = RayAngleInside2ndLayer;
  double Z0 = RxDepth;
  double Z1 = TxDepth;

  if(AirOrIce==0){
    ReceiveAngle= asin((MultiRayAirIceRefraction::Getnz_ice(Z1)*sin(Lang))/MultiRayAirIceRefraction::Getnz_ice(Z0));
  }
  if(AirOrIce==1){
    ReceiveAngle= asin((MultiRayAirIceRefraction::Getnz_air(Z1)*sin(Lang))/MultiRayAirIceRefraction::Getnz_air(Z0));   
  }
  
  //std::cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/MultiRayAirIceRefraction::pi)<<" deg"<<std::endl;
  
  ////calculate the distance of the point of incidence on the 2ndLayer surface and also the value of the L parameter of the solution
  Lvalue=nzRx*sin(ReceiveAngle);

  // auto t2b = std::chrono::high_resolution_clock::now();
  // auto t1c = std::chrono::high_resolution_clock::now();
  
  x1=GetRayHorizontalPath(A, RxDepth, TxDepth, Lvalue, AirOrIce);
  //std::cout<<"The hit point horizontal distance is from the Rx target "<<x1<<" m  on the surface"<<std::endl;

  // auto t2c = std::chrono::high_resolution_clock::now();
  // auto t1d = std::chrono::high_resolution_clock::now();
  
  ////calculate the propagation time in 2ndLayer 
  RayTimeIn2ndLayer=GetRayPropagationTime(A, RxDepth, TxDepth, Lvalue, AirOrIce);
  //std::cout<<"The propagation time in 2ndLayer is: "<<RayTimeIn2ndLayer<<" s"<<std::endl;

  if(AirOrIce==0){
    x1_Geo=GetRayGeometricPath(A, RxDepth, TxDepth, Lvalue, AirOrIce);
    //x1_Geo=0;
  }
  if(AirOrIce==1){
    x1_Geo=GetRayGeometricPath(A, RxDepth, TxDepth, Lvalue, AirOrIce);
    //x1_Geo=0;
  }
  //auto t2d = std::chrono::high_resolution_clock::now();
  
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

  // auto durationa = std::chrono::duration_cast<std::chrono::nanoseconds>( t2a - t1a ).count();
  // auto durationb = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  // auto durationc = std::chrono::duration_cast<std::chrono::nanoseconds>( t2c - t1c ).count();
  // auto durationd = std::chrono::duration_cast<std::chrono::nanoseconds>( t2d - t1d ).count();

  // std::cout<<"total time taken by the script to do a: "<<durationa<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do b: "<<durationb<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do c: "<<durationc<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do d: "<<durationd<<" ns"<<std::endl;  
  
  output[0]=x1;
  output[1]=ReceiveAngle*(180/MultiRayAirIceRefraction::pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  output[4]=x1_Geo;
  
  return output;
}

////This function flattens out 2d std::vectors into 1d std::vectors
std::vector<double> MultiRayAirIceRefraction::flatten(const std::vector<std::vector<double> >& v) {
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
  double *output=new double[5*MultiRayAirIceRefraction::MaxLayers+2];

  //auto t1a = std::chrono::high_resolution_clock::now();  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MultiRayAirIceRefraction::MaxLayers;ilayer>-1;ilayer--){
    if(AirTxHeight<MultiRayAirIceRefraction::ATMLAY[ilayer]/100 && AirTxHeight>=MultiRayAirIceRefraction::ATMLAY[ilayer-1]/100){
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
  for(int ilayer=0;ilayer<MultiRayAirIceRefraction::MaxLayers;ilayer++){
    if(IceLayerHeight>=MultiRayAirIceRefraction::ATMLAY[ilayer]/100 && IceLayerHeight<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100){
      //cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
      ilayer=100;
    }
    if(ilayer<MultiRayAirIceRefraction::MaxLayers){
      skiplayer++;
    }
  }
  int SkipLayersBelow=skiplayer;

  //auto t2a = std::chrono::high_resolution_clock::now();
  //auto t1b = std::chrono::high_resolution_clock::now();  
  
  double StartAngle=0;
  double StartHeight=0;
  double Start_nh=0;
  double StopHeight=0;

  std::vector <double> TotalHorizontalDistance;
  std::vector <double> TotalGeometricPath;
  std::vector <double> ReceiveAngle;
  std::vector <double> Lvalue;
  std::vector <double> PropagationTime;
 
  //int ipoints=0;
  for(int ilayer=MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
    
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100-0.00001;
    }
    
    ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
    Start_nh=MultiRayAirIceRefraction::Getnz_air(StartHeight);//gsl_spline_eval(spline, StartHeight, accelerator);
    
    ////Set the stopping height of the ray for propogation for that layer
    if(ilayer==(SkipLayersBelow-1)+1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      StopHeight=IceLayerHeight;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      StopHeight=MultiRayAirIceRefraction::ATMLAY[ilayer]/100;
    }
    
    ////If this is the first layer then set the initial launch angle of the ray through the layers
    if(ilayer==MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){
      StartAngle=180-LaunchAngleAir;
    }
    //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A_air<<" ,B="<<B_air[ilayer]<<" ,C="<<C_air[ilayer]<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<" , UserLaunchAngle "<<LaunchAngleAir<<std::endl;
    
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    if(ilayer==MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){ 
      //auto t1c = std::chrono::high_resolution_clock::now();  
      //cout<<"in layer "<<ilayer<<endl;
      double* GetHitPar=MultiRayAirIceRefraction::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
      //auto t2c = std::chrono::high_resolution_clock::now();
      
      TotalHorizontalDistance.push_back(GetHitPar[0]);
      ReceiveAngle.push_back(GetHitPar[1]);
      Lvalue.push_back(GetHitPar[2]);
      PropagationTime.push_back(GetHitPar[3]);
      TotalGeometricPath.push_back(GetHitPar[4]);
      StartAngle=GetHitPar[1];
      delete []GetHitPar;
      //cout<<ilayer<<" "<<path<<endl;
      //auto durationc = std::chrono::duration_cast<std::chrono::nanoseconds>( t2c - t1c ).count();
      //std::cout<<"total time taken by the script to do c: "<<durationc<<" ns"<<std::endl;
    }
    if(ilayer<MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){
      Lvalue.push_back(Lvalue[0]);
      double nzStopHeight=Getnz_air(StopHeight);
      double RecAng=asin(Lvalue[0]/nzStopHeight);
      RecAng=RecAng*(180/MultiRayAirIceRefraction::pi);
      ReceiveAngle.push_back(RecAng);
      double THD=GetRayHorizontalPath(A_air, StopHeight, StartHeight, Lvalue[0], 1);
      TotalHorizontalDistance.push_back(THD);
      double PropTime=GetRayPropagationTime(A_air, StopHeight, StartHeight, Lvalue[0], 1);
      PropagationTime.push_back(PropTime);
      double GeoPath=GetRayGeometricPath(A_air, StopHeight, StartHeight, Lvalue[0], 1);
      TotalGeometricPath.push_back(GeoPath);
      StartAngle=RecAng;
      //cout<<ilayer<<" "<<path<<endl;
    }
    
    //cout<<ilayer<<" "<<TotalHorizontalDistance[ipoints]<<" "<<ReceiveAngle[ipoints]<<" "<<Lvalue[ipoints]<<" "<<PropagationTime[ipoints]<<endl;
    
    //ipoints++;
    ////dont forget to delete the pointer!
    
  }
  //cout<<"GeoPath is "<<path<<endl;
  
  for(int i=0;i<Lvalue.size();i++){
    output[0+i*5]=TotalHorizontalDistance[i];
    output[1+i*5]=ReceiveAngle[i];
    output[2+i*5]=Lvalue[i];
    output[3+i*5]=PropagationTime[i];
    output[4+i*5]=TotalGeometricPath[i];
    //output[4+i*5]=0;

    //cout<<"lval size "<<i<<" "<<output[MultiRayAirIceRefraction::MaxLayers*i+4]<<" "<<MultiRayAirIceRefraction::MaxLayers*i+4<<" "<<output[MultiRayAirIceRefraction::MaxLayers*i+3]<<endl;
  }
  output[5*MultiRayAirIceRefraction::MaxLayers+1]=Lvalue.size();
  //auto t2b = std::chrono::high_resolution_clock::now();

  // for(int i=0;i<5*MultiRayAirIceRefraction::MaxLayers+2;i++){
  //   cout<<"check array A "<<i<<" "<<output[i]<<endl;
  // }
  
  // auto durationa = std::chrono::duration_cast<std::chrono::nanoseconds>( t2a - t1a ).count();
  // auto durationb = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  
  // std::cout<<"total time taken by the script to do a: "<<durationa<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do b: "<<durationb<<" ns"<<std::endl;

  return output;
}

////Get Propogation parameters for ray propagating in ice
double * MultiRayAirIceRefraction::GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth, double Lvalue){
  double *output=new double[5];

  double StartAngle=IncidentAngleonIce;
  double StartDepth=0.0;
  double StopDepth=AntennaDepth;
  double nzStopDepth=Getnz_ice(StopDepth);
  
  double TotalHorizontalDistance=GetRayHorizontalPath(A_ice, StopDepth, StartDepth, Lvalue, 0);
  double ReceiveAngle=asin(Lvalue/nzStopDepth)*(180/MultiRayAirIceRefraction::pi);
  double PropagationTime=GetRayPropagationTime(A_ice, StopDepth, StartDepth, Lvalue, 0);
  double TotalGeometricPath=GetRayGeometricPath(A_ice, StopDepth, StartDepth, Lvalue, 0);
  
  output[0]=TotalHorizontalDistance;
  output[1]=ReceiveAngle;
  output[2]=Lvalue;
  output[3]=PropagationTime;
  output[4]=TotalGeometricPath;

  return output;
}


////This function will be minimised to get the launch angle in air. It looks at the total horizontal distance and tries to find a ray in air and ice which can cover that horizontal distance perfectly. It takes into account the discontinuity in refractive index at the air ice boundary
double MultiRayAirIceRefraction::MinimizeforLaunchAngle(double x, void *params){

  struct MultiRayAirIceRefraction::MinforLAng_params *p= (struct MultiRayAirIceRefraction::MinforLAng_params *) params;
  double AirTxHeight = p->airtxheight;
  double IceLayerHeight = p->icelayerheight;
  double AntennaDepth = p->antennadepth;
  double HorizontalDistance = p->horizontaldistance;
  //std::cout<<"values are "<<AirTxHeight<<" "<<IceLayerHeight<<" "<<AntennaDepth<<" "<<HorizontalDistance<<std::endl;

  //auto t1a = std::chrono::high_resolution_clock::now();
  
  double TotalHorizontalDistanceinAir=0;
  double IncidentAngleonIce=0;
  double Lvalue=0;  
  double * GetResultsAir=GetAirPropagationPar(x,AirTxHeight,IceLayerHeight);
  TotalHorizontalDistanceinAir=0;
  int FilledLayers=GetResultsAir[5*MultiRayAirIceRefraction::MaxLayers+1];
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[0+i*5];
  }
  IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*5];
  Lvalue=GetResultsAir[2];
  delete [] GetResultsAir;

  // auto t2a = std::chrono::high_resolution_clock::now();
  // auto t1b = std::chrono::high_resolution_clock::now();
  
  double TotalHorizontalDistanceinIce=0;
  if(AntennaDepth!=0){
    double * GetResultsIce=GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth, Lvalue);
    TotalHorizontalDistanceinIce=GetResultsIce[0];
    delete [] GetResultsIce;
  }else{
    TotalHorizontalDistanceinIce=0;
  }

  //std::cout<<TotalHorizontalDistanceinIce<<" "<<TotalHorizontalDistanceinAir<<" "<< HorizontalDistance<<std::endl;
  double checkmin=(HorizontalDistance-(TotalHorizontalDistanceinIce + TotalHorizontalDistanceinAir) );

  // auto t2b = std::chrono::high_resolution_clock::now();
  // auto durationa = std::chrono::duration_cast<std::chrono::nanoseconds>( t2a - t1a ).count();
  // auto durationb = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  
  // std::cout<<"total time taken by the script to do a: "<<durationa<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do b: "<<durationb<<" ns"<<std::endl;
  
  return checkmin;
}

///This function loads in the GDAS atmosphere file. It calls the other functions to load in the tabulated refractive index values and the sea level refractive index value from the file. It also reads the mass overburden A,B and C values from the file
int MultiRayAirIceRefraction::MakeAtmosphere(){
   
  ////Fill in the n(h) and h arrays and ATMLAY and a,b and c (these 3 are the mass overburden parameters) from the data file
  MultiRayAirIceRefraction::readATMpar();
  MultiRayAirIceRefraction::readnhFromFile();
  
  ////Flatten out the height and the refractive index std::vectors to be used for setting the up the spline interpolation.
  std::vector <double> flattened_h_data=flatten(MultiRayAirIceRefraction::h_data);
  std::vector <double> flattened_nh_data=flatten(MultiRayAirIceRefraction::nh_data);

  ////Set up the GSL cubic spline interpolation. This used for interpolating values of refractive index at different heights.
  MultiRayAirIceRefraction::accelerator =  gsl_interp_accel_alloc();
  MultiRayAirIceRefraction::spline = gsl_spline_alloc (gsl_interp_cspline,flattened_h_data.size());
  gsl_spline_init(MultiRayAirIceRefraction::spline, flattened_h_data.data(), flattened_nh_data.data(), flattened_h_data.size());
 
  MultiRayAirIceRefraction::FillInAirRefractiveIndex();

  // flattened_h_data.clear();
  // flattened_nh_data.clear();
  
  return 0;
}

////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& geometricalPathLengthInIce, double& geometricalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){
  
  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice

  double thR=0;
  if(AntennaDepth<0){
    //324 3124 -3000 +200
    //cout<<"in here 1 "<<AirTxHeight-IceLayerHeight-AntennaDepth<<" "<<AirTxHeight<<" "<<IceLayerHeight<<" "<<AntennaDepth<<" atan "<<HorizontalDistance<<" "<<(AirTxHeight-IceLayerHeight-AntennaDepth)<<" "<<HorizontalDistance/(AirTxHeight-IceLayerHeight-AntennaDepth)<<endl;
    thR=180-(atan( HorizontalDistance/(AirTxHeight-IceLayerHeight-AntennaDepth) )*(180.0/MultiRayAirIceRefraction::pi) );
  }
  if(AntennaDepth>=0){
    //cout<<"in here 2 "<<AirTxHeight<<" "<<(IceLayerHeight+AntennaDepth)<<endl;
    thR=180-(atan( HorizontalDistance/(AirTxHeight-(IceLayerHeight+AntennaDepth)) )*(180.0/MultiRayAirIceRefraction::pi) );
  }

  //cout<<"in here too 1"<<endl;
  double dummy[20];
  MultiRayAirIceRefraction::Air2IceRayTracing(AirTxHeight, HorizontalDistance, IceLayerHeight, AntennaDepth,thR, dummy);
  //cout<<"in here too 2"<<endl;

  opticalPathLengthInIce=dummy[5]*100;
  opticalPathLengthInAir=dummy[6]*100;
  geometricalPathLengthInIce=dummy[15]*100;
  geometricalPathLengthInAir=dummy[14]*100;
  
  launchAngle=dummy[10]*(MultiRayAirIceRefraction::pi/180);
  horizontalDistanceToIntersectionPoint=dummy[2]*100;
  reflectionCoefficientS=dummy[12];
  reflectionCoefficientP=dummy[13];
  
  bool CheckSolution=false;
  double checkminimisation=dummy[1]-HorizontalDistance;

  //cout<<"raytrace arguments are "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<AntennaDepth<<" "<<IceLayerHeight<<" "<<thR<<endl;
  //cout<<"in here too 3"<<endl;
  if((fabs(dummy[1]-HorizontalDistance)/HorizontalDistance<0.01 && HorizontalDistance<=100) || (fabs(dummy[1]-HorizontalDistance)<1 && HorizontalDistance>100)){
    CheckSolution=true;
  }
  if(dummy[1]<0){
    CheckSolution=false;
  }
  //cout<<"in here too 4"<<endl;

  //cout<<"raytrace results are "<<dummy[1]<<" "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<reflectionCoefficientS<<" "<<reflectionCoefficientP<<" "<<CheckSolution<<endl;
  
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

void MultiRayAirIceRefraction::FindClosestAirTxHeight(double ParValue, int &RStartIndex1, int &REndIndex1, double &ClosestVal1,  int &RStartIndex2, int &REndIndex2, double &ClosestVal2, int AntennaNumber){

  int CurrentHeightStep=floor((ParValue-LoopStopHeight)/HeightStepSize);

  int Index=TotalHeightSteps-CurrentHeightStep-1;
  int MaxAngleBin=Index*TotalAngleSteps+TotalAngleSteps-1;
  int MinAngleBin=Index*TotalAngleSteps+0;

  // for(int i=MinAngleBin-2;i<MaxAngleBin+2;i++){
  //   cout<<"check values "<<i<<" "<<AllTableAllAntData[AntennaNumber][0][i]<<" "<<AllTableAllAntData[AntennaNumber][1][i]<<endl;
  // }

  // for(int i=MinAngleBin-2-TotalAngleSteps;i<MaxAngleBin+2-TotalAngleSteps;i++){
  //   cout<<"check values too "<<i<<" "<<AllTableAllAntData[AntennaNumber][0][i]<<" "<<AllTableAllAntData[AntennaNumber][1][i]<<endl;
  // }

  ////FindFirstBin
  double val=-0.001;
  int StartBin=MaxAngleBin;
  bool wentinside=false;
  while((val!=0 && val<0.01) || std::isnan(val)==true){
    val=AllTableAllAntData[AntennaNumber][1][StartBin];
    StartBin--;
    wentinside=true;
  }
  if(wentinside==true){
    StartBin=StartBin+1;
  }

  val=-0.001;
  int EndBin=MinAngleBin;
  wentinside=false;
  while((val!=0 && val<0.01) || std::isnan(val)==true){
    val=AllTableAllAntData[AntennaNumber][1][EndBin];
    EndBin++;
    wentinside=true;
  }
  if(wentinside==true){
    EndBin=EndBin-1;
  }
  
  RStartIndex1=EndBin;
  REndIndex1=StartBin;
  ClosestVal1=fabs(AllTableAllAntData[AntennaNumber][0][Index]-ParValue);
  //cout<<"bin values are "<<EndBin<<" "<<StartBin<<" "<<AllTableAllAntData[AntennaNumber][0][Index]<<" "<<ParValue<<endl;

  // ////FindSecondBin
  // val=-0.001;
  // StartBin=MaxAngleBin-TotalAngleSteps;
  // if(StartBin<0){
  //   StartBin=MaxAngleBin+TotalAngleSteps;
  // }
  // wentinside=false;
  // while((val!=0 && val<0.01) || std::isnan(val)==true){
  //   val=AllTableAllAntData[AntennaNumber][1][StartBin];
  //   StartBin--;
  //   wentinside=true;
  // }
  // if(wentinside==true){
  //   StartBin=StartBin+1;
  // }

  // val=-0.001;
  // EndBin=MinAngleBin-TotalAngleSteps;
  // if(EndBin<0){
  //   EndBin=MinAngleBin+TotalAngleSteps;
  // }
  // wentinside=false;
  // while((val!=0 && val<0.01) || std::isnan(val)==true){
  //   val=AllTableAllAntData[AntennaNumber][1][EndBin];
  //   EndBin++;
  //   wentinside=true;
  // }
  // if(wentinside==true){
  //   EndBin=EndBin-1;
  // }
  
  // RStartIndex2=EndBin;
  // REndIndex2=StartBin;

  RStartIndex2=RStartIndex1-TotalAngleSteps;
  REndIndex2=REndIndex1-TotalAngleSteps;
 
  if(RStartIndex2<0){
    RStartIndex2=RStartIndex1+TotalAngleSteps;
  }
  if(REndIndex2<0){
    REndIndex2=REndIndex1+TotalAngleSteps;
  }
  
  ClosestVal2=fabs(AllTableAllAntData[AntennaNumber][0][Index]-ParValue);
  //cout<<"bin values are "<<EndBin<<" "<<StartBin<<" "<<AllTableAllAntData[AntennaNumber][0][Index]<<" "<<ParValue<<endl;
  
}

int MultiRayAirIceRefraction::FindClosestTHD(double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber){
  int MidIndex=floor((StartIndex+EndIndex)/2);

  for(int i=0;i<8;i++){
    if(EndIndex-StartIndex>=3){
      MidIndex=floor((StartIndex+EndIndex)/2);
      //cout<<"the initial values are  "<<AllTableAllAntData[AntennaNumber][1][StartIndex]<<" "<<AllTableAllAntData[AntennaNumber][1][MidIndex]<<" "<<AllTableAllAntData[AntennaNumber][1][EndIndex]<<" "<<ParValue<<endl;
      //cout<<"the first start and stop index are "<<StartIndex<<" "<<MidIndex<<" "<<EndIndex<<endl;
      if(AllTableAllAntData[AntennaNumber][1][MidIndex]-ParValue>0){
	//cout<<"we are here 1 "<<endl;
	StartIndex=MidIndex;
      }
      if(AllTableAllAntData[AntennaNumber][1][MidIndex]-ParValue<0){
	//cout<<"we are here 2 "<<endl;
	EndIndex=MidIndex;
      }
    }
  }
  
  double minimum=100000000000;
  int index2=0;
  double minval=0;
  for(int ipnt=StartIndex;ipnt<EndIndex+1;ipnt++){
    minval=fabs(AllTableAllAntData[AntennaNumber][1][ipnt]-ParValue);
    if(minval<minimum && AllTableAllAntData[AntennaNumber][1][ipnt]>ParValue){
      minimum=minval;
    }else{
      index2=ipnt; 
      break;
    }
  }
  double index1=index2-1;

  minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index2]);
  if(minimum>fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1])){
    minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1]);
  }
  
  RStartIndex=index1;
  REndIndex=index2;
  ClosestVal=minimum;
  
  return 0;
}

////Interpolate the value of the given parameter for a given TxHeight and THD
int MultiRayAirIceRefraction::GetParValues(double AntennaNumber, double AirTxHeight, double TotalHorizontalDistance, double IceLayerHeight, double &AirTxHeight1, double Par1[20],double &AirTxHeight2, double Par2[20]){

  //cout<<"inside interpolator"<<endl;
  //std::cout<<" "<<std::endl;

  //cout<<"arguments to interpolator are "<<AntennaNumber<<" "<<AirTxHeight<<" "<<TotalHorizontalDistance<<" "<<IceLayerHeight<<endl;
  
  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
  MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
  MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];

  //cout<<" TotalTableEntries "<<TotalTableEntries<<" MaxAirTxHeight "<<MaxAirTxHeight<<" MinAirTxHeight "<<MinAirTxHeight<<endl;
  
  int startindex[20];
  int endindex[20];
  double closestvalue[20];
  double x1,x2,y1,y2;
  double MaxTotalHorizontalDistance=0;
  double MinTotalHorizontalDistance=0;

  ////Find out the AirTxheight1 bin range for the given Tx height
  FindClosestAirTxHeight(AirTxHeight, startindex[0], endindex[0], closestvalue[0],startindex[2], endindex[2], closestvalue[2], AntennaNumber);

  //cout<<" FindClosestAirTxHeight: IceLayerHeight "<<IceLayerHeight<<" AirTxHeight "<<AirTxHeight<<" "<<0<<" TotalTableEntries "<<TotalTableEntries<<" startindex[0] "<<startindex[0]<<" endindex[0] "<<endindex[0]<<" closestvalue[0] "<<closestvalue[0]<<" AntennaNumber "<<AntennaNumber<<endl;
  
  ////Set the first Tx height
  AirTxHeight1=AllTableAllAntData[AntennaNumber][0][startindex[0]];

  //cout<<"AirTxHeight1 "<<AirTxHeight1<<endl;
  
  ////Find the maximum and minimum possible values of THD for AirTxHeight1. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
  MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[0]];
  MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[0]];

  //cout<<" MaxTotalHorizontalDistance "<<MaxTotalHorizontalDistance<<" MinTotalHorizontalDistance "<<MinTotalHorizontalDistance<<endl;
  
  //cout<<" minimum distance is 1 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[0]+1]<<endl;
  if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){
    ////Find out the THD bin range for the given THD and AirTxHeight1
    FindClosestTHD(TotalHorizontalDistance, startindex[0], endindex[0], startindex[1], endindex[1], closestvalue[1], AntennaNumber);
    
    if(closestvalue[1]!=0){////If the given THD does not match a bin value
      x1=AllTableAllAntData[AntennaNumber][1][startindex[1]];
      x2=AllTableAllAntData[AntennaNumber][1][endindex[1]];
      for(int ipar=0;ipar<9;ipar++){
	y1=AllTableAllAntData[AntennaNumber][1+ipar][startindex[1]];
	y2=AllTableAllAntData[AntennaNumber][1+ipar][endindex[1]];
	Par1[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
      }
    }
    if(closestvalue[1]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the same and there is no interpolation needed
      startindex[1]=startindex[1]+1;
      endindex[1]=startindex[1];
      for(int ipar=0;ipar<9;ipar++){
	Par1[ipar]=AllTableAllAntData[AntennaNumber][1+ipar][startindex[1]];
      }
    }
    
  }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
  ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
    //cout<<" we are here 1 "<<endl;
    //   double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);

  //   if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
  //     double THDLimit=FindExtrapolationLimit(startindex[0],TotalHorizontalDistance,AntennaNumber);
  //     double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
  //     double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
  //     if(DistancePercentage<DistancePercentageLimit){
  // 	for(int ipar=0;ipar<7;ipar++){
  // 	  Par1[ipar]=Extrapolate(1+ipar, startindex[0], TotalHorizontalDistance, AntennaNumber);
  // 	}
  //     }else{
	for(int ipar=0;ipar<9;ipar++){
	  Par1[ipar]=-pow(10,9);
	}
  //     }
  //   }
    
  }
  
  // ////Find out the Txheight2 bin range for the given Tx height
  // startindex[2]=startindex[0]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)+1;
  // endindex[2]=startindex[2]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize);
  
  // for(int i=startindex[2]-2;i<endindex[2]+2;i++){
  //   cout<<"check values too  "<<i<<" "<<AllTableAllAntData[AntennaNumber][0][i]<<" "<<AllTableAllAntData[AntennaNumber][1][i]<<endl;
  // }
  
  if(closestvalue[0]!=0 && AirTxHeight>MinAirTxHeight && startindex[2]<TotalTableEntries){////If the AirTxHeight does not exactly match a table value and AirTxHeight1 is not equal to AirTxHeight2 AND AirTxHeight1 is not the minimum possible height in the table

    ////Set the second Tx height
    AirTxHeight2=AllTableAllAntData[AntennaNumber][0][startindex[2]];
    ////Find the maximum and minimum possible values of THD for AirTxHeight2. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
    MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[2]];
    MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[2]];
    //cout<<"part 2 "<<TotalHorizontalDistance<<" MaxTotalHorizontalDistance "<<MaxTotalHorizontalDistance<<" MinTotalHorizontalDistance "<<MinTotalHorizontalDistance<<endl;
    //cout<<" minimum distance is 2 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[2]+1]<<endl;
    if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){
      //cout<<"we are here now "<<endl;
      
      FindClosestTHD(TotalHorizontalDistance, startindex[2], endindex[2], startindex[3], endindex[3], closestvalue[2], AntennaNumber);
       
      if(closestvalue[2]!=0){////If the given THD does not match a bin value
	x1=AllTableAllAntData[AntennaNumber][1][startindex[3]];
	x2=AllTableAllAntData[AntennaNumber][1][endindex[3]];
	for(int ipar=0;ipar<9;ipar++){
	  y1=AllTableAllAntData[AntennaNumber][1+ipar][startindex[3]];
	  y2=AllTableAllAntData[AntennaNumber][1+ipar][endindex[3]];
	  Par2[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
	}
      }
      if(closestvalue[2]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the sohame and there is no interpolation needed
	startindex[3]=startindex[3]+1;
	endindex[3]=startindex[3];
	for(int ipar=0;ipar<9;ipar++){
	  Par2[ipar]=AllTableAllAntData[AntennaNumber][1+ipar][startindex[3]];
	}
      }

    }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
      //   ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
      //cout<<" we are here 2 "<<endl;
    
      //   double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);
      //   if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
      // 	double THDLimit=FindExtrapolationLimit(startindex[2],TotalHorizontalDistance,AntennaNumber);
      // 	double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
      // 	double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
      // 	if(DistancePercentage<DistancePercentageLimit){
      // 	  for(int ipar=0;ipar<7;ipar++){
      // 	    Par2[ipar]=Extrapolate(1+ipar, startindex[2], TotalHorizontalDistance, AntennaNumber);
      // 	  }
      // 	}else{
      for(int ipar=0;ipar<9;ipar++){
	Par2[ipar]=-pow(10,9);
      }
    }
    //   }
      
    // }
    //cout<<" we are here 3 "<<endl;
  }else{////If AirTxHeight is exactly equal to a bin value then AirTxHeight1 and AirTxHeight2 are equal and we dont need to do interpolation for AirTxHeight2 as we already have done it for AirTxHeight1
    ////Also if AirTxHeight1 is equal to the minimum possible height then stop right there and dont the interpolation for AirTxHeight2
    //cout<<" we are here 4 "<<endl;
    AirTxHeight2=AirTxHeight1;
    for(int ipar=0;ipar<9;ipar++){
      Par2[ipar]=Par1[ipar];
    }
  }

  //cout<<" we are here 5 "<<endl;
  
  return 0; 
}

////This functions reads in the antenna tables and interpolates (or extrapolates) from the table to provide output value for raytracing
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint_Table(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, int AntennaNumber, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& geometricalPathLengthInIce, double& geometricalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){

  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice
  
  double opticalPathLengthInIceB=0;
  double opticalPathLengthInAirB=0;
  double launchAngleB=0;
  double horizontalDistanceToIntersectionPointB=0;
  double reflectionCoefficientSB=0;
  double reflectionCoefficientPB=0;
  
  // double thR=180-(atan( HorizontalDistance/(AirTxHeight-IceLayerHeight+AntennaDepth) )*(180.0/MultiRayAirIceRefraction::pi) );
  //std::cout<<"Table values "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<std::endl;
  // // double PropagationTimeIce=0;
  // // double PropagationTimeAir=0;
  // // double LaunchAngleAir=0;
  // // double TotalHorizontalDistanceinAir=0;
  // //double TotalHorizontalDistance=0;

  // // ////0 is AirTxHeight, 1 is THD, 2 is Optical Path in Ice, 3 is Optical Path in Air, 4 is Launch Angle in Air, 5 is THD Air, 6 is Refl Coeff S, 7 is Refl Coeff For P
  // bool CheckSolution=true;  
  // if(AirTxHeight<100000+1 && AirTxHeight>IceLayerHeight && thR<180 && thR>90){
  //   double THD=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 1);
  //   opticalPathLengthInIce=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 2)*100;
  //   opticalPathLengthInAir=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 3)*100;
  //   launchAngle=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 4)*(MultiRayAirIceRefraction::pi/180);
  //   horizontalDistanceToIntersectionPoint=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 5)*100;
  //   reflectionCoefficientS=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 6);
  //   reflectionCoefficientP=MultiRayAirIceRefraction::GetInterpolatedValue(AirTxHeight, thR, 7);
   
  //   if(THD<0.001 || std::isnan(THD)==true){
  //     CheckSolution=false;
  //   }
  // }else{
  //   CheckSolution=false;
  // }

  // //std::cout<<"results are "<<CheckSolution<<" "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<reflectionCoefficientS<<" "<<reflectionCoefficientP<<" params are "<<HorizontalDistance<<" "<<AirTxHeight<<" "<<IceLayerHeight<<" "<<AntennaDepth<<" thR "<<thR<<std::endl;

  //AntennaNumber=0;
  bool CheckSolution=true;  
   
    int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
    MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
    MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];
    
    double x1=0,x2=0,y1=0,y2=0;
    double AirTxHeight1,AirTxHeight2;
    double Par1[10];
    double Par2[10];
    double ParInterpolatedValues[10];
    ///0 is OpticalPathIce
    ///1 is OpticalPathAir
    ///2 is LaunchAngleAir
    ///3 is THDAir

    if(AirTxHeight<=MaxAirTxHeight && AirTxHeight>=MinAirTxHeight && AirTxHeight>0){
      //cout<<"in here "<<endl;
   
      GetParValues(AntennaNumber,AirTxHeight,HorizontalDistance,IceLayerHeight,AirTxHeight1,Par1,AirTxHeight2,Par2);

      //cout<<"we are here 6"<<endl;
      x1=AirTxHeight1;
      x2=AirTxHeight2;
      for(int ipar=0;ipar<9;ipar++){
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
	    ipar=8;
	  }

        }
        ParInterpolatedValues[ipar]=GetInterpolatedParValue;
      }
      //cout<<" we are here 7 "<<endl;
    
    }
    
    double THD=ParInterpolatedValues[0];
    opticalPathLengthInIce=ParInterpolatedValues[1]*100;
    opticalPathLengthInAir=ParInterpolatedValues[2]*100;
    geometricalPathLengthInIce=ParInterpolatedValues[8]*100;
    geometricalPathLengthInAir=ParInterpolatedValues[7]*100;
    launchAngle=ParInterpolatedValues[3]*(MultiRayAirIceRefraction::pi/180);
    horizontalDistanceToIntersectionPoint=ParInterpolatedValues[4]*100;
    reflectionCoefficientS=ParInterpolatedValues[5];
    reflectionCoefficientP=ParInterpolatedValues[6];
    //cout<<"inter values are "<<ParInterpolatedValues[0]<<" "<<ParInterpolatedValues[1]<<" "<<ParInterpolatedValues[2]<<" "<<ParInterpolatedValues[3]<<" "<<ParInterpolatedValues[4]<<" "<<ParInterpolatedValues[5]<<" "<<ParInterpolatedValues[6]<<" "<<ParInterpolatedValues[7]<<" "<<ParInterpolatedValues[8]<<endl;
    bool CheckSolBool=false;
    if( (y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9)) ){
      CheckSolBool=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(SrcHeightASL*100, HorizontalDistanceToRx*100, RxDepthBelowIceBoundary*100, IceLayerHeight*100, geometricalPathLengthInIce, geometricalPathLengthInAir, opticalPathLengthInIce, opticalPathLengthInAir, launchAngle, horizontalDistanceToIntersectionPoint,reflectionCoefficientS,reflectionCoefficientP);
      //cout<<"in here 2 "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<y1<<" "<<y2<<endl;
    }

    if(y2==-pow(10,9) && y1==-pow(10,9)){
      //std::cout<<"The given Total Horizontal Distance bigger than the maximum percentage of the Total Horizontal Distance range for the given AirTxHeight! Cannot extrapolate! "<<HorizontalDistance<<" "<<AirTxHeight<<" "<<IceLayerHeight<<std::endl;
      CheckSolution=false;
    }
    if( ((y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9))) && CheckSolBool==false ){
      //std::cout<<"The minimiser failed for the weird extrapolation case"<<std::endl;
      CheckSolution=false;    
    }
    if(AirTxHeight>MaxAirTxHeight){
      //std::cout<<"AirTx Height is greater than the maximum possible height in GDAS!"<<std::endl;
      CheckSolution=false;
    }
    if(AirTxHeight<MinAirTxHeight){
      //std::cout<<"AirTx Height is less than the minimum possible height in the table!"<<std::endl;
      CheckSolution=false;
    }
    if(AirTxHeight<0){
      //std::cout<<"AirTx Height is less than the zero!"<<std::endl;
      CheckSolution=false;
    }
    if(launchAngle<0){
      //std::cout<<"launchAngle is less than the zero!"<<std::endl;
      CheckSolution=false;
    }
    
    if((fabs(THD-HorizontalDistance)/HorizontalDistance>0.01 && HorizontalDistance<=100) || (fabs(THD-HorizontalDistance)>1 && HorizontalDistance>100)){
      CheckSolution=false;
    }
    
    //std::cout<<"results 1 "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<CheckSolution<<" "<<THD<<" "<<HorizontalDistance<<std::endl;
    if(CheckSolution==false){
      opticalPathLengthInIce=0;
      opticalPathLengthInAir=0;
      launchAngle=0;
      horizontalDistanceToIntersectionPoint=0;
    }
    
    //std::cout<<"results 1 "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<std::endl;
    //std::cout<<"results 2 "<<opticalPathLengthInIceB<<" "<<opticalPathLengthInAirB<<" "<<launchAngleB<<" "<<horizontalDistanceToIntersectionPointB<<std::endl;
  
  return CheckSolution;
}

void MultiRayAirIceRefraction::Air2IceRayTracing(double AirTxHeight, double HorizontalDistance, double IceLayerHeight,double AntennaDepth, double StraightAngle, double dummy[20]){
  //cout<<"max layers are "<<MultiRayAirIceRefraction::MaxLayers<<endl;
    
  //std::cout<<"parameters are "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<std::endl;
  ////For recording how much time the process took
  //auto t1b = std::chrono::high_resolution_clock::now();  

  gsl_function F1;
  struct MultiRayAirIceRefraction::MinforLAng_params params1;
  if(AntennaDepth>=0){
    //cout<<"here 1 "<<endl;
    IceLayerHeight=AntennaDepth+IceLayerHeight;
    AntennaDepth=0;
    params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  }
  if(AntennaDepth<0){
    //cout<<"here 2 "<<endl;
    params1 = { AirTxHeight, IceLayerHeight, -AntennaDepth, HorizontalDistance};
  }
  F1.function = & MultiRayAirIceRefraction::MinimizeforLaunchAngle;
  F1.params = &params1;
 
  ////Set the initial angle limits for the minimisation
  double startanglelim=90;
  double endanglelim=180;

  startanglelim=StraightAngle-16;
  endanglelim=StraightAngle;
  //cout<<"angles are "<<startanglelim<<" "<<endanglelim<<" "<<StraightAngle<<endl;
  if(startanglelim<90.001){
    startanglelim=90.001;
    ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
    bool checknan=false;
    double TotalHorizontalDistanceinAirt=0;
    int FilledLayerst=0;
    while(checknan==false && startanglelim>89.9){
      double *GetResultsAirTest1= MultiRayAirIceRefraction::GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
      TotalHorizontalDistanceinAirt=0;
      FilledLayerst=GetResultsAirTest1[5*MultiRayAirIceRefraction::MaxLayers+1];
      for(int i=0;i<FilledLayerst;i++){
	TotalHorizontalDistanceinAirt+=GetResultsAirTest1[0+i*5];
      }
      delete []GetResultsAirTest1;
    
      if((isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0) || startanglelim>endanglelim-0.1){    
	checknan=true;
      }else{
	startanglelim=startanglelim+0.05;
      }
    }
  }
  //cout<<"angles are "<<startanglelim<<" "<<endanglelim<<" "<<StraightAngle<<endl;
  if(endanglelim<90.001 && endanglelim>90.00){
    //startanglelim=90.05;
    endanglelim=90.05;
  }
  //cout<<"angles are "<<startanglelim<<" "<<endanglelim<<" "<<StraightAngle<<endl;
  //cout<<"in here too too 1"<<endl;
  
  //auto t1b_air = std::chrono::high_resolution_clock::now();
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir= MultiRayAirIceRefraction::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  //std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;

  //cout<<"in here too too 1a"<<endl;
  
  double * GetResultsAir= MultiRayAirIceRefraction::GetAirPropagationPar(LaunchAngleAir,AirTxHeight,IceLayerHeight);
  // for(int i=0;i<5*MultiRayAirIceRefraction::MaxLayers+2;i++){
  //   cout<<"check array B "<<i<<" "<<GetResultsAir[i]<<endl;
  // }
  int FilledLayers=GetResultsAir[5*MultiRayAirIceRefraction::MaxLayers+1];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  double TotalGeometricPathinAir=0;
  //cout<<"in here too too 1b"<<endl;
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[0+i*5];
    PropagationTimeAir+=GetResultsAir[3+i*5];
    TotalGeometricPathinAir+=GetResultsAir[4+i*5];
    //cout<<"check layers "<<i<<" "<<TotalGeometricPathinAir<<" "<<GetResultsAir[MultiRayAirIceRefraction::MaxLayers*i+4]<<" "<<MultiRayAirIceRefraction::MaxLayers*i+4<<" "<<GetResultsAir[3+i*MultiRayAirIceRefraction::MaxLayers]<<endl;
  }
  //cout<<"in here too too 1c"<<endl;
  double Lvalue=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*5];  
  delete [] GetResultsAir;

  //cout<<"in here too too 2"<<endl;
  
  //auto t2b_air = std::chrono::high_resolution_clock::now();
  
  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<std::endl;
  // std::cout<<"LvalueAir for "<<Lvalue<<std::endl;
  // std::cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<std::endl;

  //auto t1b_ice = std::chrono::high_resolution_clock::now();

  double TotalHorizontalDistanceinIce=0;
  double IncidentAngleonAntenna=0;
  //double LvalueIce=0;
  double PropagationTimeIce=0;
  double TotalGeometricPathinIce=0;

  if(AntennaDepth<0){
    //cout<<"here 3 "<<endl;
    double * GetResultsIce=MultiRayAirIceRefraction::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, -AntennaDepth,Lvalue);
    TotalHorizontalDistanceinIce=GetResultsIce[0];
    IncidentAngleonAntenna=GetResultsIce[1];
    //double LvalueIce=GetResultsIce[2];
    PropagationTimeIce=GetResultsIce[3];
    TotalGeometricPathinIce=GetResultsIce[4];
    delete [] GetResultsIce;
  }
  //auto t2b_ice = std::chrono::high_resolution_clock::now();

  //cout<<"in here too too 3"<<endl;
  
  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<std::endl;
  // std::cout<<"LvalueIce "<<Lvalue<<std::endl;
  // std::cout<<"PropagationTimeIce "<<PropagationTimeIce<<" ns"<<std::endl;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;
  
  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice + Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistance "<<TotalHorizontalDistance<<" m"<<std::endl;
  // std::cout<<"TotalPropagationTime "<<TotalPropagationTime<<" ns"<<std::endl;

  // auto t2b = std::chrono::high_resolution_clock::now();
  // auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();
  // auto durationb_ice = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b_ice - t1b_ice ).count();
  // auto durationb_air = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b_air - t1b_air ).count();

  // durationb=durationb/1000;
  // durationb_ice=durationb_ice;
  // durationb_air=durationb_air;
  // std::cout<<"total time taken by the script to do solution calcuation: "<<durationb<<" ms"<<std::endl;
  // std::cout<<"total time taken by the script to do solution calcuation for Ice: "<<durationb_ice<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do solution calcuation for Air: "<<durationb_air<<" ns"<<std::endl;
  // std::cout<<" "<<std::endl;
  
  dummy[0]=AirTxHeight;
  dummy[1]=(TotalHorizontalDistance);
  dummy[2]=TotalHorizontalDistanceinAir;
  dummy[3]=TotalHorizontalDistanceinIce;
  dummy[4]=TotalPropagationTime*MultiRayAirIceRefraction::spedc;
  dummy[5]=PropagationTimeIce*MultiRayAirIceRefraction::spedc;
  dummy[6]=PropagationTimeAir*MultiRayAirIceRefraction::spedc;
  dummy[7]=TotalPropagationTime;
  dummy[8]=PropagationTimeIce;
  dummy[9]=PropagationTimeAir;
  dummy[10]=LaunchAngleAir;
  dummy[11]=IncidentAngleonIce;
  dummy[12]=MultiRayAirIceRefraction::Refl_S(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);
  dummy[13]=MultiRayAirIceRefraction::Refl_P(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);
  dummy[14]=TotalGeometricPathinAir;
  dummy[15]=TotalGeometricPathinIce;
  //std::cout<<"in raytracer "<<dummy[0]<<" "<<dummy[1]<<" "<<dummy[2]<<" "<<dummy[3]<<std::endl;
  //std::cout<<"in raytracer "<<dummy[0]<<" "<<dummy[1]<<" "<<dummy[5]<<" "<<dummy[6]<<" "<<dummy[10]<<" "<<dummy[2]<<" "<<dummy[12]<<" "<<dummy[13]<<std::endl;
  // if(AirTxHeight==5001 && StraightAngle==114.4){
  //   cout<<AirTxHeight<<" "<<StraightAngle<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<" "<<(AirTxHeight-IceLayerHeight+AntennaDepth)*tan((180-StraightAngle)*(MultiRayAirIceRefraction::pi/180.0))<<endl;  
  // cout<<dummy[0]<<" "<<dummy[1]<<" "<<dummy[5]<<" "<<dummy[6]<<" "<<dummy[10]<<" "<<dummy[2]<<" "<<dummy[12]<<" "<<dummy[13]<<endl;
  // }
}

void MultiRayAirIceRefraction::MakeTable(double IceLayerHeight, double AntennaDepth){ 

  IceLayerHeight=IceLayerHeight/100;
  AntennaDepth=AntennaDepth/100;

  std::cout<<"making the table now "<<AntennaDepth<<" "<<IceLayerHeight<<std::endl;

  ////Print out the entry number, the Tx height, ice layer height, Tx height above the icelayer height, total horizontal distance on surface, total horizontal distance in ice, RayLaunchAngle at Tx, incident angle on ice and recievd angle in ice at the antenna inside this file
  MultiRayAirIceRefraction::MakeAtmosphere();
  
  MultiRayAirIceRefraction::GridStartH=IceLayerHeight+1;
  MultiRayAirIceRefraction::GridStopH=100000;

  MultiRayAirIceRefraction::GridWidthH=MultiRayAirIceRefraction::GridStopH - MultiRayAirIceRefraction::GridStartH;
  
  MultiRayAirIceRefraction::TotalStepsH_O=(MultiRayAirIceRefraction::GridWidthH/MultiRayAirIceRefraction::GridStepSizeH_O)+1;
  MultiRayAirIceRefraction::TotalStepsTh_O=(MultiRayAirIceRefraction::GridWidthTh/MultiRayAirIceRefraction::GridStepSizeTh_O)+1;

  MultiRayAirIceRefraction::GridPoints=MultiRayAirIceRefraction::TotalStepsH_O*MultiRayAirIceRefraction::TotalStepsTh_O;

  GridPositionH.resize(MultiRayAirIceRefraction::TotalStepsH_O);
  GridPositionTh.resize(MultiRayAirIceRefraction::TotalStepsTh_O);
    
  //////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();  
  
  for(int ih=0;ih<MultiRayAirIceRefraction::TotalStepsH_O;ih++){
    for(int ith=0;ith<MultiRayAirIceRefraction::TotalStepsTh_O;ith++){
      
      double h=MultiRayAirIceRefraction::GridStartH+MultiRayAirIceRefraction::GridStepSizeH_O*ih;
      double th=MultiRayAirIceRefraction::GridStartTh+MultiRayAirIceRefraction::GridStepSizeTh_O*ith;
      
      if(ih==TotalStepsH_O-1){
	h=GridStopH;
      }
        
      if(ith==TotalStepsTh_O-1){
	th=GridStopTh;
      }      

      GridPositionH[ih]=h;
      GridPositionTh[ith]=th;
      
      double TotalHorizontalDistance=(h-IceLayerHeight+AntennaDepth)*tan((180-th)*(MultiRayAirIceRefraction::pi/180.0));
      double dummy[20];
      
      MultiRayAirIceRefraction::Air2IceRayTracing(h, TotalHorizontalDistance, IceLayerHeight, AntennaDepth,th, dummy);
      
      if((fabs(dummy[1]-TotalHorizontalDistance)/TotalHorizontalDistance<0.01 && TotalHorizontalDistance<=100) || (fabs(dummy[1]-TotalHorizontalDistance)<1 && TotalHorizontalDistance>100)){
	MultiRayAirIceRefraction::GridZValue[0].push_back(dummy[0]);//0 is AirTxHeight
	MultiRayAirIceRefraction::GridZValue[1].push_back((dummy[1]));//1 is THD
	MultiRayAirIceRefraction::GridZValue[2].push_back(dummy[5]);//2 is Optical Path in Ice
	MultiRayAirIceRefraction::GridZValue[3].push_back(dummy[6]);//3 is Optical Path in Air
	MultiRayAirIceRefraction::GridZValue[4].push_back(dummy[10]);//4 is Launch Angle in Air
	MultiRayAirIceRefraction::GridZValue[5].push_back(dummy[2]);//5 is Total Horizontal Distance in Air
	MultiRayAirIceRefraction::GridZValue[6].push_back(dummy[12]);//6 is Refl Coeff S
	MultiRayAirIceRefraction::GridZValue[7].push_back(dummy[13]);//7 is Refl Coeff P
      }else{
	MultiRayAirIceRefraction::GridZValue[0].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[1].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[2].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[3].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[4].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[5].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[6].push_back(-1000);
	MultiRayAirIceRefraction::GridZValue[7].push_back(-1000);
      }
    }
  }

  auto t2b = std::chrono::high_resolution_clock::now();
  double Duration = std::chrono::duration_cast<std::chrono::seconds>( t2b - t1b ).count();

  std::cout<<"The table took "<<Duration<<" s to make"<<std::endl;
}



double MultiRayAirIceRefraction::GetInterpolatedValue(double hR, double thR, int rtParameter){

  int MinDistBin[20];
  double MinDist[20];

  double sum1=0;
  double sum2=0;
  double NewZValue=0;
  
  double minHbin=round((hR-MultiRayAirIceRefraction::GridStartH)/MultiRayAirIceRefraction::GridStepSizeH_O);
  double minThbin=round((thR-MultiRayAirIceRefraction::GridStartTh)/MultiRayAirIceRefraction::GridStepSizeTh_O);
  
  int newHbin=minHbin;
  int newThbin=minThbin;
  
  int count=0;
  if(minHbin<=1){
    minHbin=1;
  }
  if(minThbin<=1){
    minThbin=1;
  }
  
  if(minHbin+1>MultiRayAirIceRefraction::TotalStepsH_O){
    minHbin=MultiRayAirIceRefraction::TotalStepsH_O-2;
  }
  if(minThbin+1>MultiRayAirIceRefraction::TotalStepsTh_O){
    minThbin=MultiRayAirIceRefraction::TotalStepsTh_O-2;
  }
  
  int startbinH=minHbin-1;
  int endbinH=minHbin+1;
  int startbinTh=minThbin-1;
  int endbinTh=minThbin+1;    
    
  newHbin=((minHbin-1));
  newThbin=minThbin-1;
  int newich=(minHbin-1)*TotalStepsTh_O+(minThbin-1);
  double minDist1=fabs(((hR-GridPositionH[newHbin])*(hR-GridPositionH[newHbin])+(thR-GridPositionTh[newThbin])*(thR-GridPositionTh[newThbin])));

  newHbin=((minHbin+1));
  newThbin=minThbin+1;
  newich=(minHbin+1)*TotalStepsTh_O+(minThbin+1);
  double minDist2=fabs(((hR-GridPositionH[newHbin])*(hR-GridPositionH[newHbin])+(thR-GridPositionTh[newThbin])*(thR-GridPositionTh[newThbin]))); 

  startbinH=minHbin-1;
  endbinH=minHbin+1;
  startbinTh=minThbin-1;
  endbinTh=minThbin+1;
   
  sum1=0;
  sum2=0;
  NewZValue=-1000;
 
  for(int ixn=startbinH;ixn<endbinH;ixn++){
    for(int izn=startbinTh;izn<endbinTh;izn++){

      newHbin=ixn;
      newThbin=izn;
      newich=(ixn)*TotalStepsTh_O+(izn);      
      
      if(newich>=0 && newich<MultiRayAirIceRefraction::GridPoints && ixn<MultiRayAirIceRefraction::TotalStepsH_O && izn<MultiRayAirIceRefraction::TotalStepsTh_O && ixn>=0 && izn>=0){
	MinDist[count]=fabs(((hR-MultiRayAirIceRefraction::GridPositionH[ixn])*(hR-MultiRayAirIceRefraction::GridPositionH[ixn])+(thR-MultiRayAirIceRefraction::GridPositionTh[izn])*(thR-MultiRayAirIceRefraction::GridPositionTh[izn])));
	MinDistBin[count]=newich;

	
	if(GridZValue[rtParameter][MinDistBin[count]]!=-1000){
	  sum1+=(1.0/MinDist[count])*GridZValue[rtParameter][MinDistBin[count]];
	  sum2+=(1.0/MinDist[count]);
	  NewZValue=sum1/sum2;
	}else{
	  sum1+=0;
	  sum2+=0;
	  NewZValue=-1000;
	}
	    
	if(MinDist[count]==0){
	  if(GridZValue[rtParameter][MinDistBin[count]]!=-1000){
	    NewZValue=GridZValue[rtParameter][MinDistBin[count]];
	    izn=minThbin+3;
	    ixn=minHbin+3;
	  }else{
	    NewZValue=-1000;
	    izn=minThbin+3;
	    ixn=minHbin+3;
	  }
	}
	count++;
      }
    }
  }
  
  return NewZValue;
  
}

void MultiRayAirIceRefraction::GetRayTracingSolutions(double RayLaunchAngleInAir,double AirTxHeight, double IceLayerHeight, double AntennaDepth, double dummy[20], bool &InIce){

  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MultiRayAirIceRefraction::MaxLayers;ilayer>-1;ilayer--){
    //std::cout<<ilayer<<" "<<MultiRayAirIceRefraction::ATMLAY[ilayer]/100<<" "<<MultiRayAirIceRefraction::ATMLAY[ilayer-1]/100<<std::endl;
    if(AirTxHeight<MultiRayAirIceRefraction::ATMLAY[ilayer]/100 && AirTxHeight>=MultiRayAirIceRefraction::ATMLAY[ilayer-1]/100){
      //std::cout<<"Tx Height is in this layer with a height range of "<<MultiRayAirIceRefraction::MultiRayAirIceRefraction::ATMLAY[ilayer]/100<<" m to "<<MultiRayAirIceRefraction::MultiRayAirIceRefraction::ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<std::endl;
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
  for(int ilayer=0;ilayer<MultiRayAirIceRefraction::MaxLayers;ilayer++){
    //std::cout<<ilayer<<" "<<MultiRayAirIceRefraction::ATMLAY[ilayer]/100<<" "<<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100<<std::endl;
    if(IceLayerHeight>=MultiRayAirIceRefraction::ATMLAY[ilayer]/100 && IceLayerHeight<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100){
      //std::cout<<"Ice Layer is in the layer with a height range of "<<MultiRayAirIceRefraction::ATMLAY[ilayer]/100<<" m to "<<MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<std::endl;
      ilayer=100;
    }
    if(ilayer<MultiRayAirIceRefraction::MaxLayers){
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
  double TotalGeometricPathInAir=0;
  double TimeInAir=0;
      
  ////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
  //std::cout<<"Fitting the atmosphere refrative index profile with multiple layers and propogate the ray"<<std::endl;
  for(int ilayer=MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
    //std::cout<<B_air<<std::endl;
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=MultiRayAirIceRefraction::ATMLAY[ilayer+1]/100-0.00001;
    }
	
    ////Since we have the starting height now we can find out the refactive index at that height
    Start_nh=MultiRayAirIceRefraction::Getnz_air(StartHeight);
	
    ////Set the staopping height of the ray for propogation for that layer
    if(ilayer==(SkipLayersBelow-1)+1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      StopHeight=IceLayerHeight;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      StopHeight=MultiRayAirIceRefraction::ATMLAY[ilayer]/100;
    }
	
    ////If this is the first layer then set the initial launch angle of the ray through the layers. I calculate the final launch angle by doing 180-RayLaunchAngleInAir since my raytracer only works with 0 to 90 deg. Setting an angle of 95 deg w.r.t to the vertical where 0 is up means that my raytraces takes in an launch angle of 85.
    if(ilayer==MultiRayAirIceRefraction::MaxLayers-SkipLayersAbove-1){
      StartAngle=180-RayLaunchAngleInAir;
    }
    //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<std::endl;
	
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    double* GetHitPar=MultiRayAirIceRefraction::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
    TotalHorizontalDistanceInAir+=GetHitPar[0];
    StartAngle=GetHitPar[1];
    TimeInAir+=GetHitPar[3];
    TotalGeometricPathInAir+=GetHitPar[4];
    
    ////dont forget to delete the pointer!
    delete []GetHitPar;
  }

  
  double IncidentAngleonIce=StartAngle;
  //std::cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<std::endl;

  ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
  double TotalHorizontalDistanceInIce=0;
  double RecievdAngleInIce=0;
  //double LvalueIce=GetHitPar[2];
  double TimeInIce=0;
  double TotalGeometricPathInIce=0;
  
  if(InIce==true){

    //cout<<" we are in here "<<InIce<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Section for propogating the ray through the ice
            
    ////Set the starting depth of the ray for propogation to at the ice surface
    double StartDepth=0.0;
    ////Since we have the starting height of the ice layer we can find out the refactive index of air at that height
    Start_nh=MultiRayAirIceRefraction::Getnz_air(IceLayerHeight);
	
    ////Set the stopping depth of the ray for propogation to be the depth of the antenna
    StopHeight=AntennaDepth;
    ////Set the initial launch angle or the angle of incidence of the ray
    StartAngle=IncidentAngleonIce;
    //std::cout<<"Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingDepth="<<StartDepth<<" ,StoppingDepth="<<AntennaDepth<<" ,RayLaunchAngle="<<StartAngle<<std::endl;
    
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel through ice to hit the antenna
    //// The angle of reciept/incidence at the end at the antenna
    //// The value of the L parameter for whole atmosphere fit
    double *GetHitPar=MultiRayAirIceRefraction::GetLayerHitPointPar(Start_nh, -AntennaDepth,StartDepth, StartAngle, 0);
      
    ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
    TotalHorizontalDistanceInIce=GetHitPar[0];
    RecievdAngleInIce=GetHitPar[1];
    //double LvalueIce=GetHitPar[2];
    TimeInIce=GetHitPar[3];
    TotalGeometricPathInIce=GetHitPar[4];

    delete[] GetHitPar;
  }
  //std::cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceInIce<<std::endl;
      
  //if(isnan(TotalHorizontalDistanceInAir)==false){
	
  ////define dummy/temporary variables for storing data
  for(int idum=0;idum<18;idum++){
    dummy[idum]=0;
  }
  
  dummy[0]=0;
  dummy[1]=AirTxHeight;
  dummy[2]=TotalHorizontalDistanceInAir + TotalHorizontalDistanceInIce;
  dummy[3]=TotalHorizontalDistanceInAir;
  dummy[4]=TotalHorizontalDistanceInIce;
  dummy[5]=(TimeInIce+TimeInAir)*MultiRayAirIceRefraction::spedc;
  dummy[6]=TimeInAir*MultiRayAirIceRefraction::spedc;
  dummy[7]=TimeInIce*MultiRayAirIceRefraction::spedc;
  dummy[8]=(TimeInIce+TimeInAir)*pow(10,9);
  dummy[9]=TimeInAir*pow(10,9);
  dummy[10]=TimeInIce*pow(10,9);
  dummy[11]=RayLaunchAngleInAir;
  dummy[12]=IncidentAngleonIce;
  dummy[13]=RecievdAngleInIce;
  dummy[14]=MultiRayAirIceRefraction::Refl_S(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);
  dummy[15]=MultiRayAirIceRefraction::Refl_P(IncidentAngleonIce*(MultiRayAirIceRefraction::pi/180.0), IceLayerHeight);
  dummy[16]=TotalGeometricPathInAir;
  dummy[17]=TotalGeometricPathInIce;
}

int MultiRayAirIceRefraction::MakeRayTracingTable(double AntennaDepth, double IceLayerHeight, int AntennaNumber){

  bool InIce=true;
  if(AntennaDepth<0){
    InIce=true;
  }
  if(AntennaDepth>=0){
    //cout<<"we are here "<<endl;
    InIce=false;
  }
  
  std::vector<std::vector <double>> AllTableData;
  
  ////convert cm to m
  AntennaDepth=AntennaDepth/100;
  IceLayerHeight=IceLayerHeight/100;
  
  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();
  
  ////Print out the entry number, the Tx height, ice layer height, Tx height above the icelayer height, total horizontal distance on surface, total horizontal distance in ice, RayLaunchAngle at Tx, incident angle on ice and recievd angle in ice at the antenna inside this file
  MultiRayAirIceRefraction::MakeAtmosphere();
 
  ////Define variables for the loop over Tx height and ray launch angle
  double RayLaunchAngleInAir=0;////Set zero for now and 0 deg straight up. This variable defines the initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 deg is straight up
  //double AirTxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];////Maximum height available with the refractive index data
  double AirTxHeight=100000;////Maximum height available with the refractive index data
  
  // ////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
  // double AngleStepSize=0.5;
  // double LoopStartAngle=92;
  // double LoopStopAngle=178;
  // int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize);

  ////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
  LoopStartHeight=AirTxHeight;
  LoopStopHeight=IceLayerHeight;
  if(InIce==true){
    LoopStopHeight=IceLayerHeight;
  }else{
    //cout<<"we are here "<<endl;
    LoopStopHeight=IceLayerHeight+AntennaDepth;
  }

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
  std::vector <double> temp9;
  std::vector <double> temp10;
  
  int ifileentry=0;
  //cout<<" our final angle is "<<LoopStartAngle+AngleStepSize*TotalAngleSteps<<endl;
  ////Start looping over the Tx Height and Launch angle values
  for(int ihei=0;ihei<TotalHeightSteps;ihei++){
    AirTxHeight=LoopStartHeight-HeightStepSize*ihei;
    
    if(AirTxHeight>0){
      double dummy[20];
      for(int iang=0;iang<TotalAngleSteps;iang++){
	RayLaunchAngleInAir=LoopStartAngle+AngleStepSize*iang;
	//if(RayLaunchAngleInAir>=startanglelim){	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Section for propogating the ray through the atmosphere
    
	if(AirTxHeight!=LoopStopHeight && ihei==TotalHeightSteps-1){
	  AirTxHeight=LoopStopHeight;
	}
	if(iang==TotalAngleSteps-1){
	  RayLaunchAngleInAir=LoopStopAngle;
	}

	
	GetRayTracingSolutions(RayLaunchAngleInAir,AirTxHeight,LoopStopHeight,AntennaDepth,dummy,InIce);

	// cout<<"table values are "<<RayLaunchAngleInAir<<" "<<AirTxHeight<<" "<<LoopStopHeight<<" "<<AntennaDepth<<" "<<InIce<<" dummies are "<<dummy[1]<<" "<<dummy[2]<<" "<<dummy[7]<<" "<<dummy[6]<<" "<<dummy[11]<<" "<<dummy[3]<<" "<<dummy[14]<<endl;
	
	temp1.push_back(dummy[1]);///AirTx Height
	temp2.push_back(dummy[2]);///THDTotal   0       
	temp3.push_back(dummy[7]);///OpticalPathIce 1
	temp4.push_back(dummy[6]);///OpticalPathAir 2
	temp5.push_back(dummy[11]);///LaunchAngleAir 3         
	temp6.push_back(dummy[3]);///THDAir 4
	temp7.push_back(dummy[14]);///ReflectionCoefficientS 5
	temp8.push_back(dummy[15]);///ReflectionCoefficientP 6
	temp9.push_back(dummy[16]);///GeometricPathAir 7
	temp10.push_back(dummy[17]);///GeometricPathIce 8
	
	//cout<<RayLaunchAngleInAir<<" "<<AirTxHeight<<" "<<LoopStopHeight<<" "<<AntennaDepth<<" "<<InIce<<" dummies are "<<dummy[16]<<" "<<dummy[17]<<endl;
	
	ifileentry++;
     
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
  AllTableData.push_back(temp9);
  AllTableData.push_back(temp10);

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
  auto durationb = std::chrono::duration_cast<std::chrono::seconds>( t2b - t1b ).count();

  //durationb=durationb/1000000;
  std::cout<<"total time taken by the script to generate the table: "<<durationb<<" s"<<std::endl;
  return 0;
  
}
