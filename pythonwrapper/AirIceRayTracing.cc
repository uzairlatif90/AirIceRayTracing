#include "AirIceRayTracing.h"

////This Function reads in the values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file. The a,b and c values are mass overburden values and are not required in this code.
int AirIceRayTracing::readATMpar(std::string atmosFileName){
  
  ////Open the file
  std::ifstream ain(atmosFileName);
  
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
	for (int i=0; i<5; i++){ AirIceRayTracing::ATMLAY[i]=dummya[i]; }
      }    
      if(n1==1){
	for (int i=0; i<5; i++){ AirIceRayTracing::abc[i][0]=dummya[i]; }
      }
      if(n1==2){
	for (int i=0; i<5; i++){ AirIceRayTracing::abc[i][1]=dummya[i]; }
      }
      if(n1==3){
	for (int i=0; i<5; i++){ AirIceRayTracing::abc[i][2]=dummya[i]; }
      }
      n1++;
    }////end the while loop
    
    ain.close();
  }////if condition to check if file is open

  AirIceRayTracing::abc[4][0]=AirIceRayTracing::abc[3][0];
  AirIceRayTracing::abc[4][1]=AirIceRayTracing::abc[3][1];
  AirIceRayTracing::abc[4][2]=AirIceRayTracing::abc[3][2];

  AirIceRayTracing::ATMLAY[4]=150000*100;////set max possible height that can be used in cm

  //std::cout<<"abc values are "<<AirIceRayTracing::abc[4][0]<<" "<<AirIceRayTracing::abc[4][1]<<" "<<AirIceRayTracing::abc[4][2]<<std::endl;
  
  return 0;
}

int AirIceRayTracing::readnhFromFile(std::string atmosFileName){

  AirIceRayTracing::nh_data.clear();
  AirIceRayTracing::lognh_data.clear();
  AirIceRayTracing::h_data.clear();
  
  ////Open the file
  std::ifstream ain(atmosFileName);
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
	
	if(dummy1*100>=AirIceRayTracing::ATMLAY[layer]){////change the layer once the data of all the heights of that layer has been read in
	  if(layer>0){////now since the layer has finished and the temporary vectors have been filled in. Now we push the vectors in the main 2d height and refractice index vectors
	    AirIceRayTracing::h_data.push_back(temp1);
	    AirIceRayTracing::nh_data.push_back(temp2);
	    AirIceRayTracing::lognh_data.push_back(temp3);

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
      AirIceRayTracing::h_data.push_back(temp1);
      AirIceRayTracing::nh_data.push_back(temp2);
      AirIceRayTracing::lognh_data.push_back(temp3);
      ////clear the vectors now for storing the next layer
      temp1.clear();
      temp2.clear();
      temp3.clear();
    }
    layer++;
    
    ain.close();
  }////if condition to check if file is open

  ////The file reading condition "while (getline(ain,line))" reads the last the datapoint of the file twice. This is to to remove the last repeat data point in all the data arrays
  AirIceRayTracing::h_data[AirIceRayTracing::h_data.size()-1].erase(AirIceRayTracing::h_data[AirIceRayTracing::h_data.size()-1].end() - 1);
  AirIceRayTracing::nh_data[AirIceRayTracing::nh_data.size()-1].erase(AirIceRayTracing::nh_data[AirIceRayTracing::nh_data.size()-1].end() - 1);
  AirIceRayTracing::lognh_data[AirIceRayTracing::lognh_data.size()-1].erase(AirIceRayTracing::lognh_data[AirIceRayTracing::lognh_data.size()-1].end() - 1);

  AirIceRayTracing::MaxLayers=AirIceRayTracing::h_data.size()+1;////store the total number of layers present in the data

  //cout<<"max layers are "<<MaxLayers<<endl;
  
  return 0;
}

////Get the value of the B parameter for the ice refractive index model
double AirIceRayTracing::GetB_ice(double z){
  double zabs=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

////Get the value of the C parameter for the ice refractive index model
double AirIceRayTracing::GetC_ice(double z){
  double zabs=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth in ice
double AirIceRayTracing::Getnz_ice(double z){
  z=fabs(z);
  return AirIceRayTracing::A_ice+AirIceRayTracing::GetB_ice(z)*exp(-AirIceRayTracing::GetC_ice(z)*z);
}

int AirIceRayTracing::FillInAirRefractiveIndex(){
  
  double N0=0;
  for(int ilayer=0;ilayer<5;ilayer++){
    double hlow=AirIceRayTracing::ATMLAY[ilayer]/100;
    AirIceRayTracing::C_air[ilayer]=1.0/(AirIceRayTracing::abc[ilayer][2]/100);
    if(ilayer>0){
      N0=AirIceRayTracing::A_air+AirIceRayTracing::B_air[ilayer-1]*exp(-hlow*AirIceRayTracing::C_air[ilayer-1]);
    }
    if(ilayer==0){
      N0=gsl_spline_eval(spline, 0, accelerator);
    }
    AirIceRayTracing::B_air[ilayer]=((N0-1)/exp(-hlow*AirIceRayTracing::C_air[ilayer]));
  }
 
  return 0;   
}

////Get the value of the B parameter for the air refractive index model
double AirIceRayTracing::GetB_air(double z){
  double zabs=fabs(z);
  double B=0;
  int whichlayer=0;

  if(AirIceRayTracing::UseConstantRefractiveIndex==false){
    
    for(int ilayer=0;ilayer<AirIceRayTracing::MaxLayers-1;ilayer++){

      if(zabs<AirIceRayTracing::ATMLAY[ilayer+1]/100 && zabs>=AirIceRayTracing::ATMLAY[ilayer]/100){
	whichlayer=ilayer;
	ilayer=100;
      }  
    }
    if(zabs>=AirIceRayTracing::ATMLAY[AirIceRayTracing::MaxLayers-1]/100){
      whichlayer=AirIceRayTracing::MaxLayers-1;
    }
 
    B=AirIceRayTracing::B_air[whichlayer];
  }
  
  if(AirIceRayTracing::UseConstantRefractiveIndex==true){
    B=0;
  }
  return B;
}

////Get the value of the C parameter for the air refractive index model
double AirIceRayTracing::GetC_air(double z){
  double zabs=fabs(z);
  double C=0;
  int whichlayer=0;

  if(AirIceRayTracing::UseConstantRefractiveIndex==false){
    
    for(int ilayer=0;ilayer<AirIceRayTracing::MaxLayers-1;ilayer++){
      if(zabs<AirIceRayTracing::ATMLAY[ilayer+1]/100 && zabs>=AirIceRayTracing::ATMLAY[ilayer]/100){
	whichlayer=ilayer;
	ilayer=100;
      }
    }
  
    if(zabs>=AirIceRayTracing::ATMLAY[AirIceRayTracing::MaxLayers-1]/100){
      whichlayer=AirIceRayTracing::MaxLayers-1;
    }
    C=AirIceRayTracing::C_air[whichlayer];

  }

  if(AirIceRayTracing::UseConstantRefractiveIndex==true){
    C=1e-9;
  }
  return C;
}

////Get the value of refractive index model for a given depth in air
double AirIceRayTracing::Getnz_air(double z){
  double zabs=fabs(z);
  double output=0;
  if(AirIceRayTracing::UseConstantRefractiveIndex==false){
    output=AirIceRayTracing::A_air+AirIceRayTracing::GetB_air(zabs)*exp(-AirIceRayTracing::GetC_air(zabs)*zabs);
  }
  if(AirIceRayTracing::UseConstantRefractiveIndex==true){
    output=AirIceRayTracing::A_const;
  }
  return output;
}

/* E-feild Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflection coefficient. The transmission coefficient is t=1+r */
double AirIceRayTracing::Refl_S(double thetai, double IceLayerHeight){

  double Nair=AirIceRayTracing::Getnz_air(IceLayerHeight);
  double Nice=AirIceRayTracing::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;
  
  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double rS=(num/den);

  if(isnan(rS)){
    rS=1;
  }
  return (rS);
}

double AirIceRayTracing::Trans_S(double thetai, double IceLayerHeight){

  double Nair=AirIceRayTracing::Getnz_air(IceLayerHeight);
  double Nice=AirIceRayTracing::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;
  
  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double tS=1+(num/den);

  if(isnan(tS)){
    tS=0;
  }
  return (tS);
}

/* E-feild Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflection coeffient. The transmission coefficient is t=(n_1/n_2)*(1+r) */
double AirIceRayTracing::Refl_P(double thetai, double IceLayerHeight){
   
  double Nair=AirIceRayTracing::Getnz_air(IceLayerHeight);
  double Nice=AirIceRayTracing::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;

  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double rP=-(num)/(den);
  if(isnan(rP)){
    rP=1;
  }
  return (rP);
}

double AirIceRayTracing::Trans_P(double thetai, double IceLayerHeight){
   
  double Nair=AirIceRayTracing::Getnz_air(IceLayerHeight);
  double Nice=AirIceRayTracing::Getnz_ice(0); 
  double n1=Nair;
  double n2=Nice;

  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double tP=(1-(num/den))*(n1/n2);

  if(isnan(tP)){
    tP=0;
  }
  return (tP);
}

////Use GSL minimiser which uses Brent's Method to find root for a given function
double AirIceRayTracing::FindFunctionRoot(gsl_function F,double x_lo, double x_hi,const gsl_root_fsolver_type *T,double tolerance, int iterations)
{
  int status;
  int iter = 0, max_iter = iterations;
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
double AirIceRayTracing::fDnfR(double x,void *params){
  
  struct AirIceRayTracing::fDnfR_params *p= (struct AirIceRayTracing::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*(A+B*exp(C*x))-L*L+sqrt(A*A-L*L)*sqrt(pow(A+B*exp(C*x),2)-L*L)));;
}

// ////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
// double AirIceRayTracing::fdxdz(double x,void *params){
  
//   struct AirIceRayTracing::fdxdz_params *p= (struct AirIceRayTracing::fdxdz_params *) params;
//   double Lang = p->lang;
//   double Z0 = p->z0;
//   double Z1 = p->z1;
//   int AirOrIce = p->airorice;

//   double output=0,dumx=0;
//   if(AirOrIce==0){
//     dumx=(AirIceRayTracing::Getnz_ice(Z0)*sin(x))/AirIceRayTracing::Getnz_ice(Z1);
//   }
//   if(AirOrIce==1){
//     dumx=(AirIceRayTracing::Getnz_air(Z0)*sin(x))/AirIceRayTracing::Getnz_air(Z1);
//   }
//   //output=((dumx/sqrt(1-dumx*dumx)) - tan(Lang));
//   //cout<<"output is "<<output<<" "<<x<<endl;
//   output=dumx - sin(Lang);
  
//   return output;
// }

////The function used to calculate ray propogation time in ice
double AirIceRayTracing::ftimeD(double x,void *params){

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
double AirIceRayTracing::fpathD(double x,void *params){

  struct AirIceRayTracing::ftimeD_params *p= (struct AirIceRayTracing::ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;

  //integral sec(sin^(-1)(L/(A + B e^(C x)))) dx = (log((A + B e^(C x)) (sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + 1)) - (A log(A sqrt(A^2 - L^2) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + B sqrt(A^2 - L^2) e^(C x) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + A^2 + A B e^(C x) - L^2))/sqrt(A^2 - L^2) + (A C x)/sqrt(A^2 - L^2))/C;
  
  return (log((A + B*exp(C*x))*(sqrt((A*A + 2*A*B*exp(C*x) + B*B*exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x))) ) + 1)) - (A*log(A*sqrt(A*A - L*L)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/(( A + B*exp(C*x))*(A + B*exp(C*x)))) + B*sqrt(A*A - L*L)*exp(C*x)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x)))) + A*A + A*B*exp(C*x) - L*L))/sqrt(A*A - L*L) + (A*C*x)/sqrt(A*A - L*L))/C ;

}

double AirIceRayTracing::GetRayHorizontalPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct AirIceRayTracing::fDnfR_params params2a;
  struct AirIceRayTracing::fDnfR_params params2b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params2a = {A, AirIceRayTracing::GetB_ice(RxDepth), -AirIceRayTracing::GetC_ice(RxDepth), Lvalue};
    params2b = {A, AirIceRayTracing::GetB_ice(TxDepth), -AirIceRayTracing::GetC_ice(TxDepth), Lvalue};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params2a = {A, AirIceRayTracing::GetB_air(RxDepth), -AirIceRayTracing::GetC_air(RxDepth), Lvalue};
    params2b = {A, AirIceRayTracing::GetB_air(TxDepth), -AirIceRayTracing::GetC_air(TxDepth), Lvalue};
  }
  double x1=+AirIceRayTracing::fDnfR(RxDepth,&params2a)-AirIceRayTracing::fDnfR(TxDepth,&params2b);
  if(AirOrIce==1){
    x1*=-1;
  }
  
  return x1;
}

double AirIceRayTracing::GetRayPropagationTime(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct AirIceRayTracing::ftimeD_params params3a;
  struct AirIceRayTracing::ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, AirIceRayTracing::GetB_ice(RxDepth), -AirIceRayTracing::GetC_ice(RxDepth), AirIceRayTracing::spedc, Lvalue,0};
    params3b = {A, AirIceRayTracing::GetB_ice(TxDepth), -AirIceRayTracing::GetC_ice(TxDepth), AirIceRayTracing::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, AirIceRayTracing::GetB_air(RxDepth), -AirIceRayTracing::GetC_air(RxDepth), AirIceRayTracing::spedc, Lvalue,1};
    params3b = {A, AirIceRayTracing::GetB_air(TxDepth), -AirIceRayTracing::GetC_air(TxDepth), AirIceRayTracing::spedc, Lvalue,1};
  }
  double RayTimeIn2ndLayer=+AirIceRayTracing::ftimeD(RxDepth,&params3a)-AirIceRayTracing::ftimeD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  
  return RayTimeIn2ndLayer;
}

double AirIceRayTracing::GetRayGeometricPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce){
  
  struct AirIceRayTracing::ftimeD_params params3a;
  struct AirIceRayTracing::ftimeD_params params3b;
  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    params3a = {A, AirIceRayTracing::GetB_ice(RxDepth), -AirIceRayTracing::GetC_ice(RxDepth), AirIceRayTracing::spedc, Lvalue,0};
    params3b = {A, AirIceRayTracing::GetB_ice(TxDepth), -AirIceRayTracing::GetC_ice(TxDepth), AirIceRayTracing::spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    params3a = {A, AirIceRayTracing::GetB_air(RxDepth), -AirIceRayTracing::GetC_air(RxDepth), AirIceRayTracing::spedc, Lvalue,1};
    params3b = {A, AirIceRayTracing::GetB_air(TxDepth), -AirIceRayTracing::GetC_air(TxDepth), AirIceRayTracing::spedc, Lvalue,1};
  }
  double RayGeometricPath=AirIceRayTracing::fpathD(RxDepth,&params3a)-AirIceRayTracing::fpathD(TxDepth,&params3b);//0
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
double *AirIceRayTracing::GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce){

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

  double SurfaceRayIncidentAngle=IncidentAng*(AirIceRayTracing::pi/180.0);////Angle at which the ray is incident on the second layer
  double RayAngleInside2ndLayer=0;////Use Snell's Law to find the angle of transmission in the 2ndlayer

  double A=0;
  double nzRx=0;
  double nzTx=0;
  double GSLFnLimit=0;

  if(AirOrIce==0){
    //std::cout<<"in ice"<<std::endl;
    A=AirIceRayTracing::A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //std::cout<<"in air"<<std::endl;
    A=AirIceRayTracing::A_air;
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
  // struct AirIceRayTracing::fdxdz_params params1 = {RayAngleInside2ndLayer, RxDepth, TxDepth, AirOrIce};
  // F1.function = &fdxdz;
  // F1.params = &params1;
  // //cout<<"limits are "<<RayAngleInside2ndLayer*(AirIceRayTracing::pi/180)<<" "<<GSLFnLimit*(180.0/AirIceRayTracing::pi)<<endl;
  // ReceiveAngle=AirIceRayTracing::FindFunctionRoot(F1,0.0*(AirIceRayTracing::pi/180),GSLFnLimit, gsl_root_fsolver_brent,0.00000001);

  double Lang = RayAngleInside2ndLayer;
  double Z0 = RxDepth;
  double Z1 = TxDepth;

  if(AirOrIce==0){
    ReceiveAngle= asin((AirIceRayTracing::Getnz_ice(Z1)*sin(Lang))/AirIceRayTracing::Getnz_ice(Z0));
  }
  if(AirOrIce==1){
    ReceiveAngle= asin((AirIceRayTracing::Getnz_air(Z1)*sin(Lang))/AirIceRayTracing::Getnz_air(Z0));   
  }
  
  //std::cout<<"The angle from vertical at which the target recieves the ray is "<<ReceiveAngle*(180/AirIceRayTracing::pi)<<" deg"<<std::endl;
  
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
  // F4.function = &AirIceRayTracing::fDnfR;
  // F4.params = &params2b;
  // gsl_deriv_central (&F4, TxDepth, 1e-8, &result, &abserr);
  // AngleOfEntryIn2ndLayer=atan(result)*(180.0/AirIceRayTracing::pi);
  // if(TxDepth==RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=180-ReceiveAngle;
  // }
  // if(TxDepth!=RxDepth && TMath::IsNaN(AngleOfEntryIn2ndLayer)==true){
  //   AngleOfEntryIn2ndLayer=90;
  // }
  //std::cout<<"AngleOfEntryIn2ndLayer= "<<AngleOfEntryIn2ndLayer<<" ,RayAngleInside2ndLayer="<<RayAngleInside2ndLayer*(180/AirIceRayTracing::pi)<<std::endl;

  // auto durationa = std::chrono::duration_cast<std::chrono::nanoseconds>( t2a - t1a ).count();
  // auto durationb = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  // auto durationc = std::chrono::duration_cast<std::chrono::nanoseconds>( t2c - t1c ).count();
  // auto durationd = std::chrono::duration_cast<std::chrono::nanoseconds>( t2d - t1d ).count();

  // std::cout<<"total time taken by the script to do a: "<<durationa<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do b: "<<durationb<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do c: "<<durationc<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do d: "<<durationd<<" ns"<<std::endl;  
  
  output[0]=x1;
  output[1]=ReceiveAngle*(180/AirIceRayTracing::pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  output[4]=x1_Geo;
  
  return output;
}

////This function flattens out 2d std::vectors into 1d std::vectors
std::vector<double> AirIceRayTracing::flatten(const std::vector<std::vector<double> >& v) {
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
double * AirIceRayTracing::GetAirPropagationPar(double LaunchAngleAir, double AirTxHeight, double IceLayerHeight){
  double *output=new double[5*AirIceRayTracing::MaxLayers+2];

  //auto t1a = std::chrono::high_resolution_clock::now();  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=AirIceRayTracing::MaxLayers;ilayer>-1;ilayer--){
    if(AirTxHeight<AirIceRayTracing::ATMLAY[ilayer]/100 && AirTxHeight>=AirIceRayTracing::ATMLAY[ilayer-1]/100){
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
  for(int ilayer=0;ilayer<AirIceRayTracing::MaxLayers;ilayer++){
    if(IceLayerHeight>=AirIceRayTracing::ATMLAY[ilayer]/100 && IceLayerHeight<AirIceRayTracing::ATMLAY[ilayer+1]/100){
      //cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
      ilayer=100;
    }
    if(ilayer<AirIceRayTracing::MaxLayers){
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
  for(int ilayer=AirIceRayTracing::MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
    
    ////Set the starting height of the ray for propogation for that layer
    if(ilayer==AirIceRayTracing::MaxLayers-SkipLayersAbove-1){
      ////If this is the first layer then set the start height to be the height of the source
      StartHeight=AirTxHeight;
    }else{
      ////If this is any layer after the first layer then set the start height to be the starting height of the layer
      StartHeight=AirIceRayTracing::ATMLAY[ilayer+1]/100-0.00001;
    }
    
    ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
    Start_nh=AirIceRayTracing::Getnz_air(StartHeight);//gsl_spline_eval(spline, StartHeight, accelerator);
    
    ////Set the stopping height of the ray for propogation for that layer
    if(ilayer==(SkipLayersBelow-1)+1){
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      StopHeight=IceLayerHeight;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      StopHeight=AirIceRayTracing::ATMLAY[ilayer]/100;
    }
    
    ////If this is the first layer then set the initial launch angle of the ray through the layers
    if(ilayer==AirIceRayTracing::MaxLayers-SkipLayersAbove-1){
      StartAngle=180-LaunchAngleAir;
    }
    //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A_air<<" ,B="<<B_air[ilayer]<<" ,C="<<C_air[ilayer]<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<" , UserLaunchAngle "<<LaunchAngleAir<<std::endl;
    
    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    if(ilayer==AirIceRayTracing::MaxLayers-SkipLayersAbove-1){ 
      //auto t1c = std::chrono::high_resolution_clock::now();  
      //cout<<"in layer "<<ilayer<<endl;
      double* GetHitPar=AirIceRayTracing::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
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
    if(ilayer<AirIceRayTracing::MaxLayers-SkipLayersAbove-1){
      Lvalue.push_back(Lvalue[0]);
      double nzStopHeight=Getnz_air(StopHeight);
      double RecAng=asin(Lvalue[0]/nzStopHeight);
      RecAng=RecAng*(180/AirIceRayTracing::pi);
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
  
  for(int i=0;i<Lvalue.size();i++){
    output[0+i*5]=TotalHorizontalDistance[i];
    output[1+i*5]=ReceiveAngle[i];
    output[2+i*5]=Lvalue[i];
    output[3+i*5]=PropagationTime[i];
    output[4+i*5]=TotalGeometricPath[i];
    //output[4+i*5]=0;

    //cout<<"lval size "<<i<<" "<<output[AirIceRayTracing::MaxLayers*i+4]<<" "<<AirIceRayTracing::MaxLayers*i+4<<" "<<output[AirIceRayTracing::MaxLayers*i+3]<<endl;
  }
  output[5*AirIceRayTracing::MaxLayers+1]=Lvalue.size();
  //auto t2b = std::chrono::high_resolution_clock::now();

  // for(int i=0;i<5*AirIceRayTracing::MaxLayers+2;i++){
  //   cout<<"check array A "<<i<<" "<<output[i]<<endl;
  // }
  
  // auto durationa = std::chrono::duration_cast<std::chrono::nanoseconds>( t2a - t1a ).count();
  // auto durationb = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  
  // std::cout<<"total time taken by the script to do a: "<<durationa<<" ns"<<std::endl;
  // std::cout<<"total time taken by the script to do b: "<<durationb<<" ns"<<std::endl;

  return output;
}

////Get Propogation parameters for ray propagating in ice
double * AirIceRayTracing::GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth, double Lvalue){
  double *output=new double[5];

  double StartAngle=IncidentAngleonIce;
  double StartDepth=0.0;
  double StopDepth=AntennaDepth;
  double nzStopDepth=Getnz_ice(StopDepth);
  
  double TotalHorizontalDistance=GetRayHorizontalPath(A_ice, StopDepth, StartDepth, Lvalue, 0);
  double ReceiveAngle=asin(Lvalue/nzStopDepth)*(180/AirIceRayTracing::pi);
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
double AirIceRayTracing::MinimizeforLaunchAngle(double x, void *params){

  struct AirIceRayTracing::MinforLAng_params *p= (struct AirIceRayTracing::MinforLAng_params *) params;
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
  int FilledLayers=GetResultsAir[5*AirIceRayTracing::MaxLayers+1];
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
int AirIceRayTracing::MakeAtmosphere(std::string atmosFileName){
   
  ////Fill in the n(h) and h arrays and ATMLAY and a,b and c (these 3 are the mass overburden parameters) from the data file
  AirIceRayTracing::readATMpar(atmosFileName);
  AirIceRayTracing::readnhFromFile(atmosFileName);
  
  ////Flatten out the height and the refractive index std::vectors to be used for setting the up the spline interpolation.
  std::vector <double> flattened_h_data=flatten(AirIceRayTracing::h_data);
  std::vector <double> flattened_nh_data=flatten(AirIceRayTracing::nh_data);

  ////Set up the GSL cubic spline interpolation. This used for interpolating values of refractive index at different heights.
  AirIceRayTracing::accelerator =  gsl_interp_accel_alloc();
  AirIceRayTracing::spline = gsl_spline_alloc (gsl_interp_cspline,flattened_h_data.size());
  gsl_spline_init(AirIceRayTracing::spline, flattened_h_data.data(), flattened_nh_data.data(), flattened_h_data.size());
 
  AirIceRayTracing::FillInAirRefractiveIndex();

  // flattened_h_data.clear();
  // flattened_nh_data.clear();
  
  return 0;
}

////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.
bool AirIceRayTracing::GetRayTracingSolution(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& geometricalPathLengthInIce, double& geometricalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& AngleOfIncidenceOnIce, double &RecievedAngleInIce){
  
  double AirTxHeight=SrcHeightASL;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx;////Horizontal distance
  IceLayerHeight=IceLayerHeight;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary;////Depth of antenna in the ice

  double thR=0;
  if(AntennaDepth<0){
    thR=180-(atan( HorizontalDistance/(AirTxHeight-IceLayerHeight-AntennaDepth) )*(180.0/AirIceRayTracing::pi) );
  }
  if(AntennaDepth>=0){
    thR=180-(atan( HorizontalDistance/(AirTxHeight-(IceLayerHeight+AntennaDepth)) )*(180.0/AirIceRayTracing::pi) );
  }

  double dummy[20];
  AirIceRayTracing::Air2IceRayTracing(AirTxHeight, HorizontalDistance, IceLayerHeight, AntennaDepth,thR, dummy);
  
  opticalPathLengthInIce=dummy[5];
  opticalPathLengthInAir=dummy[6];
  geometricalPathLengthInIce=dummy[14];
  geometricalPathLengthInAir=dummy[13];
  
  launchAngle=dummy[10];
  horizontalDistanceToIntersectionPoint=dummy[2];
  AngleOfIncidenceOnIce=dummy[11];
  RecievedAngleInIce=dummy[12];
  
  bool CheckSolution=false;
  double checkminimisation=dummy[1]-HorizontalDistance;

  //cout<<"raytrace arguments are "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<AntennaDepth<<" "<<IceLayerHeight<<" "<<thR<<endl;
  if((fabs(dummy[1]-HorizontalDistance)/HorizontalDistance<0.01 && HorizontalDistance<=100) || (fabs(dummy[1]-HorizontalDistance)<1 && HorizontalDistance>100)){
    CheckSolution=true;
  }
  if(dummy[1]<0){
    CheckSolution=false;
  }

  //cout<<"raytrace results are "<<dummy[1]<<" "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<transmissionCoefficientS<<" "<<transmissionCoefficientP<<" "<<CheckSolution<<endl;
  
  return CheckSolution;

}

void AirIceRayTracing::Air2IceRayTracing(double AirTxHeight, double HorizontalDistance, double IceLayerHeight,double AntennaDepth, double StraightAngle, double dummy[20]){
    
  //std::cout<<"parameters are "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<std::endl;
  ////For recording how much time the process took
  //auto t1b = std::chrono::high_resolution_clock::now();  

  gsl_function F1;
  struct AirIceRayTracing::MinforLAng_params params1;
  if(AntennaDepth>=0){
    IceLayerHeight=AntennaDepth+IceLayerHeight;
    AntennaDepth=0;
    params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  }
  if(AntennaDepth<0){
    params1 = { AirTxHeight, IceLayerHeight, -AntennaDepth, HorizontalDistance};
  }
  F1.function = & AirIceRayTracing::MinimizeforLaunchAngle;
  F1.params = &params1;
 
  ////Set the initial angle limits for the minimisation
  double startanglelim=90;
  double endanglelim=180;

  startanglelim=StraightAngle-16;
  endanglelim=StraightAngle;

  if(AirIceRayTracing::UseConstantRefractiveIndex==false){  
  
    if(startanglelim<90.001){
      startanglelim=90.001;
      ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
      bool checknan=false;
      double TotalHorizontalDistanceinAirt=0;
      int FilledLayerst=0;
      while(checknan==false && startanglelim>89.9){
	double *GetResultsAirTest1= AirIceRayTracing::GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
	TotalHorizontalDistanceinAirt=0;
	FilledLayerst=GetResultsAirTest1[5*AirIceRayTracing::MaxLayers+1];
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

  }else{
    startanglelim=90;
  }
  
  if(endanglelim<90.001 && endanglelim>90.00){
    //startanglelim=90.05;
    endanglelim=90.05;
  }

  
  //cout<<"angles are "<<startanglelim<<" "<<endanglelim<<" "<<StraightAngle<<endl;
  
  //auto t1b_air = std::chrono::high_resolution_clock::now();
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir= AirIceRayTracing::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_bisection,0.000000001,40);
  //std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  
  double * GetResultsAir= AirIceRayTracing::GetAirPropagationPar(LaunchAngleAir,AirTxHeight,IceLayerHeight);

  int FilledLayers=GetResultsAir[5*AirIceRayTracing::MaxLayers+1];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  double TotalGeometricPathinAir=0;
  
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[0+i*5];
    PropagationTimeAir+=GetResultsAir[3+i*5];
    TotalGeometricPathinAir+=GetResultsAir[4+i*5];
    //cout<<"check layers "<<i<<" "<<TotalGeometricPathinAir<<" "<<GetResultsAir[AirIceRayTracing::MaxLayers*i+4]<<" "<<AirIceRayTracing::MaxLayers*i+4<<" "<<GetResultsAir[3+i*AirIceRayTracing::MaxLayers]<<endl;
  }
  double Lvalue=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*5];  
  delete [] GetResultsAir;
  
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
    double * GetResultsIce=AirIceRayTracing::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, -AntennaDepth,Lvalue);
    TotalHorizontalDistanceinIce=GetResultsIce[0];
    IncidentAngleonAntenna=GetResultsIce[1];
    //double LvalueIce=GetResultsIce[2];
    PropagationTimeIce=GetResultsIce[3];
    TotalGeometricPathinIce=GetResultsIce[4];
    delete [] GetResultsIce;
  }
  //auto t2b_ice = std::chrono::high_resolution_clock::now();
  
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
  dummy[4]=TotalPropagationTime*AirIceRayTracing::spedc;
  dummy[5]=PropagationTimeIce*AirIceRayTracing::spedc;
  dummy[6]=PropagationTimeAir*AirIceRayTracing::spedc;
  dummy[7]=TotalPropagationTime;
  dummy[8]=PropagationTimeIce;
  dummy[9]=PropagationTimeAir;
  dummy[10]=LaunchAngleAir;
  dummy[11]=asin((AirIceRayTracing::Getnz_air(IceLayerHeight)/AirIceRayTracing::Getnz_ice(0))*sin(IncidentAngleonIce*(AirIceRayTracing::pi/180)))*(180./AirIceRayTracing::pi);//IncidentAngleonIce;
  dummy[12]=IncidentAngleonAntenna;
  dummy[13]=TotalGeometricPathinAir;
  dummy[14]=TotalGeometricPathinIce;
  
}
