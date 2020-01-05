#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multimin.h>
#include <sys/time.h>

const double pi=4.0*atan(1.0); /**< Gives back value of Pi */
const double spedc=299792458.0; /**< Speed of Light in m/s */

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

////Get the value of the B parameter for the air refractive index model
double GetB_air(double z){
  double zabs=fabs(z);
  double B=0;
  
  B=0.000382116;
  return B;
}

////Get the value of the C parameter for the air refractive index model
double GetC_air(double z){
  double zabs=fabs(z);
  double C=0;

  C=0.000156803;
  return C;
}

////Get the value of refractive index model for a given depth in air
double Getnz_air(double z){
  z=fabs(z);
  return A_air+GetB_air(z)*exp(-GetC_air(z)*z);
}

////Use GSL minimiser which uses Brent's Method to find root for a given function. This will be used to find roots wherever it is needed in my code.
double FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  //cout<<x_lo<<" "<<x_hi<<" "<<endl;
  // printf ("error: %s\n", gsl_strerror (status));  
  // printf ("using %s method\n", gsl_root_fsolver_name (s));
  // printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.0001);
      
      //printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      //if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      //}
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}


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
struct fdxdz_params { double lang, z0,z1; int airorice; };
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
  double RayAngleInside2ndLayer=0;

  double A=0;
  double nzRx=0;
  double nzTx=0;
  double GSLFnLimit=0;
  
  if(AirOrIce==0){
    //cout<<"in ice"<<endl;
    A=A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //cout<<"in air"<<endl;
    A=A_air;
    nzRx=Getnz_air(RxDepth);
    nzTx=Getnz_air(TxDepth);
  }

    if(AirOrIce==0){
    //cout<<"in ice"<<endl;
    A=A_ice;
    nzRx=Getnz_ice(RxDepth);
    nzTx=Getnz_ice(TxDepth);
  }
  if(AirOrIce==1){
    //cout<<"in air"<<endl;
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
  ReceiveAngle=FindFunctionRoot(F1,1*(pi/180),GSLFnLimit);
  
  ////calculate the distance of the point of incidence on the 2ndLayer surface and also the value of the L parameter of the solution
  Lvalue=nzRx*sin(ReceiveAngle);
  struct fDnfR_params params2a;
  struct fDnfR_params params2b;
  if(AirOrIce==0){
    //cout<<"in ice"<<endl;
    params2a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), Lvalue};
    params2b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), Lvalue};
  }
  if(AirOrIce==1){
    //cout<<"in air"<<endl;
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
    //cout<<"in ice"<<endl;
    params3a = {A, GetB_ice(RxDepth), -GetC_ice(RxDepth), spedc, Lvalue,0};
    params3b = {A, GetB_ice(TxDepth), -GetC_ice(TxDepth), spedc, Lvalue,0};
  }
  if(AirOrIce==1){
    //cout<<"in air"<<endl;
    params3a = {A, GetB_air(RxDepth), -GetC_air(RxDepth), spedc, Lvalue,1};
    params3b = {A, GetB_air(TxDepth), -GetC_air(TxDepth), spedc, Lvalue,1};
  }
  RayTimeIn2ndLayer=+ftimeD(RxDepth,&params3a)-ftimeD(TxDepth,&params3b);
  if(AirOrIce==1){
    RayTimeIn2ndLayer*=-1;
  }
  
  output[0]=x1;
  output[1]=ReceiveAngle*(180/pi);
  output[2]=Lvalue;
  output[3]=RayTimeIn2ndLayer;
  
  return output;
}

////Get Propogation parameters for ray propagating in air
double * GetAirPropagationPar(double LaunchAngle, double AirTxHeight, double IceLayerHeight){
  double *output=new double[4];

  double StartAngle=180-LaunchAngle;
  double StartHeight=AirTxHeight;
  double Start_nh=Getnz_air(StartHeight);
  double StopHeight=IceLayerHeight;
  
  double* GetHitPar=GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
  
  double TotalHorizontalDistance=GetHitPar[0];
  double ReceiveAngle=GetHitPar[1];
  double Lvalue=GetHitPar[2];
  double PropagationTime=GetHitPar[3];

  delete []GetHitPar;

  output[0]=TotalHorizontalDistance;
  output[1]=ReceiveAngle;
  output[2]=Lvalue;
  output[3]=PropagationTime;

  return output;
}

////Get Propogation parameters for ray propagating in ice
double * GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth){
  double *output=new double[4];

  double StartAngle=IncidentAngleonIce;
  double StartDepth=0.0;
  double Start_nh=Getnz_air(IceLayerHeight);
  double StopHeight=AntennaDepth;

  double* GetHitPar=GetLayerHitPointPar(Start_nh, AntennaDepth,StartDepth, StartAngle, 0);
 
  double TotalHorizontalDistance=GetHitPar[0];
  double ReceiveAngle=GetHitPar[1];
  double Lvalue=GetHitPar[2];
  double PropagationTime=GetHitPar[3];

  delete []GetHitPar;

  output[0]=TotalHorizontalDistance;
  output[1]=ReceiveAngle;
  output[2]=Lvalue;
  output[3]=PropagationTime;

  return output;
}

////This function will be minimised to get the launch angle in air. It looks at the total horizontal distance and tries to find a ray in air and ice which can cover that horizontal distance perfectly. It takes into account the discontinuity in refractive index at the air ice boundary
struct MinforLAng_params { double airtxheight, icelayerheight, antennadepth, horizontaldistance; };
double MinimizeforLaunchAngle(double x, void *params){

  struct MinforLAng_params *p= (struct MinforLAng_params *) params;
  double AirTxHeight = p->airtxheight;
  double IceLayerHeight = p->icelayerheight;
  double AntennaDepth = p->antennadepth;
  double HorizontalDistance = p->horizontaldistance;
  
  double * GetResultsAir=GetAirPropagationPar(x,AirTxHeight,IceLayerHeight);
  double TotalHorizontalDistanceinAir=GetResultsAir[0];
  double IncidentAngleonIce=GetResultsAir[1];
  delete [] GetResultsAir;

  double * GetResultsIce=GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  delete [] GetResultsIce;
 
  return TotalHorizontalDistanceinIce + TotalHorizontalDistanceinAir - HorizontalDistance;
  
}

////This function is used to measure the amount of time the code takes to run
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp (){
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void Air2IceRayTracing_wROOTplot(double AirSrcHeight, double HorizontalDistance, double IceLayerHeight, double AntennaDepth){

  ////For recording how much time the process took
  timestamp_t t0 = get_timestamp();
  
  // double AirSrcHeight=5000;////Height of the source
  // double HorizontalDistance=1000;////Horizontal distance
  // double IceLayerHeight=3000;////Height where the ice layer starts off
  // double AntennaDepth=200;////Depth of antenna in the ice  
  bool PlotRayPath=false;
  
  if(AirSrcHeight<IceLayerHeight){
    cout<<"WARNING: AirSrcHeight is less than IceLayerHeight."<<endl;
    cout<<"Setting AirSrcHeight to be 100 m above the IceLayerHeight"<<endl;
    AirSrcHeight=IceLayerHeight+100;
  }
  
  gsl_function F1;
  struct MinforLAng_params params1 = { AirSrcHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  F1.function = &MinimizeforLaunchAngle;
  F1.params = &params1;

  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir=FindFunctionRoot(F1,90,180);
  cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<endl;
  
  double * GetResultsAir=GetAirPropagationPar(LaunchAngleAir,AirSrcHeight,IceLayerHeight);
  double TotalHorizontalDistanceinAir=GetResultsAir[0];
  double IncidentAngleonIce=GetResultsAir[1];
  double LvalueAir=GetResultsAir[2];
  double PropagationTimeAir=GetResultsAir[3]*pow(10,9);
  delete [] GetResultsAir;

  cout<<"***********Results for Air************"<<endl;
  cout<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<endl;
  cout<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<endl;
  cout<<"LvalueAir "<<LvalueAir<<endl;
  cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<endl;
  
  double * GetResultsIce=GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  double IncidentAngleonAntenna=GetResultsAir[1];
  double LvalueIce=GetResultsAir[2];
  double PropagationTimeIce=GetResultsAir[3]*pow(10,9);
  delete [] GetResultsIce;

  cout<<" "<<endl;
  cout<<"***********Results for Ice************"<<endl;
  cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<endl;
  cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<endl;
  cout<<"LvalueIce "<<LvalueIce<<endl;
  cout<<"PropagationTimeIce "<<PropagationTimeIce<<" ns"<<endl;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;

  cout<<" "<<endl;
  cout<<"***********Results for Ice + Air************"<<endl;
  cout<<"TotalHorizontalDistance "<<TotalHorizontalDistance<<" m"<<endl;
  cout<<"TotalPropagationTime "<<TotalPropagationTime<<" ns"<<endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice

  if(PlotRayPath==true){
    ////Print out the ray path x and y values in a file
    ofstream aout("RayPathinAirnIce.txt");
  
    TMultiGraph *mg=new TMultiGraph();
    TMultiGraph *mgB=new TMultiGraph();
  
    TGraph *grRefracted=new TGraph();
    TGraph *grStraight=new TGraph();
    TGraph *grResidual=new TGraph();
    TGraph *grAirIce=new TGraph();
    TGraph *grIceLayer=new TGraph();

    ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
    double StraightLine_slope=tan(pi/2-LaunchAngleAir*(pi/180));
    double StraightLine_y_intercept=AirSrcHeight;
    //cout<<"StraightLine slope "<<StraightLine_slope<<" ,StraightLine_y_intercept "<<StraightLine_y_intercept<<endl;
  
    ////Define ray variables for plotting and/or storing ray path as it comes down from the atmosphere
    double Refracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray
    double LastRefracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray in the previous iteration
    double LastHeight=0;////Y coordinate variable which stores the Height of the ray in the previous iteration
    double LayerHoriOffset=0;////X axes layer offset that needs to be calculated to align the layers with each other
    double LayerStartHeight=0;////The starting height for the propagation in the layer
    double LayerStopHeight=0;////The stopping height for the propagation in the layer
    int ipoints=0;////variable for counting the total number of samples that make up the ray path

    ////Get and Set the A,B,C and L parameters for the layer
    struct fDnfR_params params2;

    ////Start looping over the layers to trace out the ray
    for(int il=0;il<1;il++){
    
      ////If this is the first layer then set the start height to be the height of the source
      LayerStartHeight=AirSrcHeight;
    
      ////If this is the last layer then set the stopping height to be the height of the ice layer
      LayerStopHeight=IceLayerHeight;
    
      //cout<<il<<" A="<<layerAs[il]<<" ,B="<<layerBs[il]<<" ,C="<<layerCs[il]<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<" ,LayerStartHeight="<<LayerStartHeight<<" ,LayerStopHeight="<<LayerStopHeight<<endl;

      ////Start tracing out the ray as it propagates through the layer
      for(double i=LayerStartHeight;i>LayerStopHeight;i--){
	params2 = {A_air, GetB_air(-i), GetC_air(-i), LvalueAir};
      
	////Calculate the x (distance) value of the ray for given y (height) value
	Refracted_x=fDnfR(-i,&params2)-fDnfR(-AirSrcHeight,&params2);

	// ////If the ray just started off in a new layer we need to offset the x values of the new layer so that the ray aligns with the previous layer.
	// if(ipoints>0 && fabs(i-LayerStartHeight)<0.0001){
	// 	LayerHoriOffset=(Refracted_x-LastRefracted_x);
	// 	//cout<<il<<" layer offset is "<<LayerHoriOffset<<endl;
	// }
	// Refracted_x=Refracted_x-LayerHoriOffset;

	////Caclulate the y value of the straight line
	double StraightLine_y=StraightLine_slope*Refracted_x+StraightLine_y_intercept;

	grAirIce->SetPoint(ipoints,Refracted_x,i);
	grRefracted->SetPoint(ipoints,Refracted_x,i);
	grStraight->SetPoint(ipoints,Refracted_x,StraightLine_y);
	grIceLayer->SetPoint(ipoints,Refracted_x,IceLayerHeight);
      
	////To make sure that the residual in air is only caculated above the ice surface
	if(StraightLine_y>=IceLayerHeight){
	  grResidual->SetPoint(ipoints,Refracted_x,i-StraightLine_y);
	}

	//cout<<ipoints<<" "<<Refracted_x<<" "<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<i-StraightLine_y<<endl;
	aout<<ipoints<<" "<<Refracted_x<<" "<<i<<endl;

	////If you want to check the the transition between different layers uncomment these lines
	// if(fabs(i-LayerStopHeight)<10){
	// 	cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<endl;;
	// }
      
	// if(fabs(i-LayerStartHeight)<10){
	// 	cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<endl;;
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

      double refractedpath=TotalHorizontalDistanceinAir-fDnfR((double)i,&params3a)+fDnfR(0,&params3b);
      grAirIce->SetPoint(ipoints,refractedpath,(double)i+IceLayerHeight);
      grIceLayer->SetPoint(ipoints,refractedpath,IceLayerHeight);
      aout<<ipoints<<" "<<refractedpath<<" "<<(double)i+IceLayerHeight<<endl;
      ipoints++;
    }

    grRefracted->SetLineColor(2);//red
    grStraight->SetLineColor(4);//blue
    grRefracted->SetLineWidth(5);
    grStraight->SetLineWidth(2);
  
    grStraight->SetTitle("Straight Line with same launch angle");
    grRefracted->SetTitle("Refracted Ray");

    TString mgtitle="Launch Angle=";
    mgtitle+=180-LaunchAngleAir;
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

    mgB->SetTitle("RayPath through Air and Ice;Distance (m); Height (m)");
    grAirIce->SetMarkerStyle(20);
    grIceLayer->SetMarkerStyle(20);
    grIceLayer->SetMarkerColor(kBlue);
    mgB->Add(grAirIce);
    mgB->Add(grIceLayer);
    TCanvas *cb=new TCanvas("cb","cb");
    cb->cd();
    mgB->Draw("ALP");
  }
  
  timestamp_t t1 = get_timestamp();
  
  double secs = (t1 - t0) / 1000000.0L;
  cout<<"total time taken by the script: "<<secs<<" s"<<endl;

}
