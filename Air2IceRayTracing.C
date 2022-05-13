#include "RayTracingFunctions.cc"

int main(int argc, char **argv){
    
  if(argc==1){
    std::cout<<"No Extra Command Line Argument Passed Other Than Program Name"<<std::endl;
    std::cout<<"Example run command: ./Air2IceRayTracing 5000 1000 3000 200"<<std::endl;
    std::cout<<"Here 5000 m is Tx Height in air in m, 1000 is the horizontal distance btw Tx in air and Rx in ice in m, 3000 is Ice Layer Height in m and 200 is the Antenna Depth in ice in m"<<std::endl;
    return 0;
  }
  if(argc<5){
    std::cout<<"More Arguments needed!"<<std::endl;
    std::cout<<"Example run command: ./Air2IceRayTracing 5000 1000 3000 200"<<std::endl;
    std::cout<<"Here 5000 m is Tx Height in air in m, 1000 is the horizontal distance btw Tx in air and Rx in ice in m, 3000 is Ice Layer Height in m and 200 is the Antenna Depth in ice in m"<<std::endl;
    return 0;
  }
  if(argc==5){
    std::cout<<"Tx Height in air is set at "<<atof(argv[1])<<" m, the horizontal distance btw Tx in air and Rx in ice is set at "<<atof(argv[2])<<" m, Ice Layer Height is set at "<<atof(argv[3]) <<" m, Antenna Depth is set at "<<atof(argv[4])<<" m"<<std::endl;
    if(atof(argv[1])<atof(argv[3])){
      std::cout<<"WARNING: AirTxHeight is less than IceLayerHeight."<<std::endl;
      std::cout<<"Please set the AirTxHeight to be above the IceLayerHeight"<<std::endl;
      return 0;
    }
  } 
  if(argc>5){
    std::cout<<"More Arguments than needed!"<<std::endl;
    std::cout<<"Example run command: ./Air2IceRayTracing 5000 1000 3000 200"<<std::endl;
    std::cout<<"Here 5000 m is Tx Height in air in m, 1000 is the horizontal distance btw Tx in air and Rx in ice in m, 3000 is Ice Layer Height in m and 200 is the Antenna Depth in ice in m"<<std::endl;
    return 0;
  }

  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();

  auto t1b_atm = std::chrono::high_resolution_clock::now();
  RayTracingFunctions::MakeAtmosphere();
  auto t2b_atm = std::chrono::high_resolution_clock::now();
  
  double AirTxHeight=atof(argv[1]);////Height of the source
  double HorizontalDistance=atof(argv[2]);////Horizontal distance
  double IceLayerHeight=atof(argv[3]);////Height where the ice layer starts off
  double AntennaDepth=atof(argv[4]);////Depth of antenna in the ice
  
  // // double AirTxHeight=5000;////Height of the source
  // // double HorizontalDistance=1000;////Horizontal distance
  // // double IceLayerHeight=3000;////Height where the ice layer starts off
  // // double AntennaDepth=200;////Depth of antenna in the ice  

  bool StoreRayPath=false;
  
  gsl_function F1;
  struct  RayTracingFunctions::MinforLAng_params params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  F1.function = & RayTracingFunctions::MinimizeforLaunchAngle;
  F1.params = &params1;
  
  ////Set the initial angle limits for the minimisation
  double startanglelim=90;
  double endanglelim=180;

  double StraightAngle=180-(atan( HorizontalDistance/(AirTxHeight-IceLayerHeight+AntennaDepth) )*(180.0/RayTracingFunctions::pi) );

  ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
  // bool checknan=false;
  // double TotalHorizontalDistanceinAirt=0;
  // int FilledLayerst=0;
  // while(checknan==false && startanglelim>89.9){
  //   double *GetResultsAirTest1= RayTracingFunctions::GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
  //   TotalHorizontalDistanceinAirt=0;
  //   FilledLayerst=GetResultsAirTest1[4*RayTracingFunctions::MaxLayers];
  //   for(int i=0;i<FilledLayerst;i++){
  //     TotalHorizontalDistanceinAirt+=GetResultsAirTest1[i*4];
  //   }
  //   delete []GetResultsAirTest1;

  //   //std::cout<<"startanglelim is "<<startanglelim<<" "<<TotalHorizontalDistanceinAirt<<std::endl;
  //   if((isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0) || startanglelim>180){
  //     checknan=true;
  //   }else{
  //     startanglelim=startanglelim+0.05;
  //   }
  // }
  
  // checknan=false;
  // while(checknan==false && endanglelim<180.1){
  //   double *GetResultsAirTest2= RayTracingFunctions::GetAirPropagationPar(endanglelim,AirTxHeight,IceLayerHeight);
  //   TotalHorizontalDistanceinAirt=0;
  //   FilledLayerst=GetResultsAirTest2[4*RayTracingFunctions::MaxLayers];
  //   for(int i=0;i<FilledLayerst;i++){
  //     TotalHorizontalDistanceinAirt+=GetResultsAirTest2[i*4];
  //   }
  //   delete []GetResultsAirTest2;

  //   //std::cout<<"endanglelim is "<<endanglelim<<" "<<TotalHorizontalDistanceinAirt<<std::endl;
  //   if((isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0) || endanglelim<90){
  //     checknan=true;
  //   }else{
  //     endanglelim=endanglelim-0.05;
  //   }
  // }
  
  startanglelim=StraightAngle-16;
  endanglelim=StraightAngle;

  if(startanglelim<90.00){
    startanglelim=90.05;
    ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
    bool checknan=false;
    double TotalHorizontalDistanceinAirt=0;
    int FilledLayerst=0;
    while(checknan==false && startanglelim>89.9){
      double *GetResultsAirTest1= RayTracingFunctions::GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
      TotalHorizontalDistanceinAirt=0;
      FilledLayerst=GetResultsAirTest1[4*RayTracingFunctions::MaxLayers];
      for(int i=0;i<FilledLayerst;i++){
	TotalHorizontalDistanceinAirt+=GetResultsAirTest1[i*4];
      }
      delete []GetResultsAirTest1;
    
      if((isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0) || startanglelim>endanglelim-1){    
	checknan=true;
      }else{
	startanglelim=startanglelim+0.05;
      }
    }
  }
  
  std::cout<<"Launch Angle search range is:  Startangle "<<startanglelim<<" ,Endangle "<<endanglelim<<std::endl;
  // //std::cout<<" "<<std::endl;

  auto t1b_air = std::chrono::high_resolution_clock::now();
  
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir= RayTracingFunctions::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  //std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  
  double * GetResultsAir= RayTracingFunctions::GetAirPropagationPar(LaunchAngleAir,AirTxHeight,IceLayerHeight);
  int FilledLayers=GetResultsAir[4*RayTracingFunctions::MaxLayers];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[i*4];
    PropagationTimeAir+=GetResultsAir[3+i*4]*pow(10,9);
  }
  double Lvalue=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*4];  
  delete [] GetResultsAir;

  auto t2b_air = std::chrono::high_resolution_clock::now();
  
  std::cout<<" "<<std::endl;
  std::cout<<"***********Results for Air************"<<std::endl;
  std::cout<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<std::endl;
  std::cout<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<std::endl;
  std::cout<<"LvalueAir for "<<Lvalue<<std::endl;
  std::cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<std::endl;

  auto t1b_ice = std::chrono::high_resolution_clock::now();
  
  double * GetResultsIce=RayTracingFunctions::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth,Lvalue);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  double IncidentAngleonAntenna=GetResultsIce[1];
  //double LvalueIce=GetResultsIce[2];
  double PropagationTimeIce=GetResultsIce[3]*pow(10,9);
  delete [] GetResultsIce;

  auto t2b_ice = std::chrono::high_resolution_clock::now();
  
  std::cout<<" "<<std::endl;
  std::cout<<"***********Results for Ice************"<<std::endl;
  std::cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<std::endl;
  std::cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<std::endl;
  std::cout<<"LvalueIce "<<Lvalue<<std::endl;
  std::cout<<"PropagationTimeIce "<<PropagationTimeIce<<" ns"<<std::endl;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;
  
  std::cout<<" "<<std::endl;
  std::cout<<"***********Results for Ice + Air************"<<std::endl;
  std::cout<<"TotalHorizontalDistance "<<TotalHorizontalDistance<<" m"<<std::endl;
  std::cout<<"TotalPropagationTime "<<TotalPropagationTime<<" ns"<<std::endl;

  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();
  auto durationb_ice = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b_ice - t1b_ice ).count();
  auto durationb_air = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b_air - t1b_air ).count();
  auto durationb_atm = std::chrono::duration_cast<std::chrono::microseconds>( t2b_atm - t1b_atm ).count();

  durationb=durationb/1000;
  durationb_ice=durationb_ice;
  durationb_air=durationb_air;
  durationb_atm=durationb_atm/1000;
  std::cout<<"total time taken by the script to do solution calcuation: "<<durationb<<" ms"<<std::endl;
  std::cout<<"total time taken by the script to do solution calcuation for Ice: "<<durationb_ice<<" ns"<<std::endl;
  std::cout<<"total time taken by the script to do solution calcuation for Air: "<<durationb_air<<" ns"<<std::endl;
  std::cout<<"total time taken by the script to do solution calcuation for Atm: "<<durationb_atm<<" ms"<<std::endl;
  std::cout<<" "<<std::endl;
  
  // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice

  ////For recording how much time the process took
  t1b = std::chrono::high_resolution_clock::now();
  
  if(StoreRayPath==true){
    ////Print out the ray path x and y values in a file
    std::ofstream aout("RayPathinAirnIce.txt");

    ////Find out how many atmosphere layers are above the source or Tx which we do not need
    int skiplayer=0;
    for(int ilayer=RayTracingFunctions::MaxLayers;ilayer>-1;ilayer--){
      //std::cout<<ilayer<<" "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" "<<RayTracingFunctions::ATMLAY[ilayer-1]/100<<std::endl;
      if(AirTxHeight<RayTracingFunctions::ATMLAY[ilayer]/100 && AirTxHeight>=RayTracingFunctions::ATMLAY[ilayer-1]/100){
	//std::cout<<"Tx Height is in this layer with a height range of "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" m to "<<RayTracingFunctions::ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<std::endl;
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
    for(int ilayer=0;ilayer<RayTracingFunctions::MaxLayers;ilayer++){
      //std::cout<<ilayer<<" "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" "<<RayTracingFunctions::ATMLAY[ilayer+1]/100<<std::endl;
      if(IceLayerHeight>=RayTracingFunctions::ATMLAY[ilayer]/100 && IceLayerHeight<RayTracingFunctions::ATMLAY[ilayer+1]/100){
	//std::cout<<"Ice Layer is in the layer with a height range of "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" m to "<<RayTracingFunctions::ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<std::endl;
	ilayer=100;
      }
      if(ilayer<RayTracingFunctions::MaxLayers){
	skiplayer++;
      }
    }
    int SkipLayersBelow=skiplayer;
    //std::cout<<"The total number of layers that need to be skipped from below is "<<skiplayer<<std::endl;
    //std::cout<<"max layers "<<RayTracingFunctions::MaxLayers<<std::endl;
  
    double Start_nh=0;
    double StartHeight=0;
    double StopHeight=0;
    double StartAngle=0;
    double TotalHorizontalDistance=0;
    std::vector <double> layerLs;////std::vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them

    double RecieveAngle=0;
    double Lvalue=0;
    
    for(int ilayer=RayTracingFunctions::MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
 
      ////Set the starting height of the ray for propogation for that layer
      if(ilayer==RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	////If this is the first layer then set the start height to be the height of the source
	StartHeight=AirTxHeight;
      }else{
	////If this is any layer after the first layer then set the start height to be the starting height of the layer
	StartHeight=RayTracingFunctions::ATMLAY[ilayer+1]/100-0.00001;
      }

      ////Since we have the starting height now we can find out the refactive index at that height from data using spline interpolation
      Start_nh=gsl_spline_eval(RayTracingFunctions::spline, StartHeight, RayTracingFunctions::accelerator);
      
      ////Set the stopping height of the ray for propogation for that layer
      if(ilayer==(SkipLayersBelow-1)+1){
	////If this is the last layer then set the stopping height to be the height of the ice layer
	StopHeight=IceLayerHeight;
      }else{
	////If this is NOT the last layer then set the stopping height to be the end height of the layer
	StopHeight=RayTracingFunctions::ATMLAY[ilayer]/100;
      }
      
      ////If this is the first layer then set the initial launch angle of the ray through the layers
      if(ilayer==RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	StartAngle=180-LaunchAngleAir;
      }
      //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<std::endl;
      
      ////Get the hit parameters from the function. The output is:
      //// How much horizontal distance did the ray travel in the layer
      //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
      //// The value of the L parameter for that layer
      if(ilayer==RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	double* GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
	TotalHorizontalDistance+=GetHitPar[0];
	RecieveAngle=GetHitPar[1];
	StartAngle=GetHitPar[1];
	////Store in the values of A,B,C and L for tha layer
	Lvalue=GetHitPar[2];
	layerLs.push_back(GetHitPar[2]);
	delete []GetHitPar;  
      }
      if(ilayer<RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	double nzStopHeight=RayTracingFunctions::Getnz_air(StopHeight);
	double RecAng=asin(Lvalue/nzStopHeight);
	RecAng=RecAng*(180/RayTracingFunctions::pi);
	double THD=RayTracingFunctions::GetRayOpticalPath(RayTracingFunctions::A_air, StopHeight, StartHeight, Lvalue, 1);
	TotalHorizontalDistance+=THD;
	StartAngle=RecAng;
      ////Store in the values of A,B,C and L for tha layer
	layerLs.push_back(Lvalue);
      }
      //std::cout<<"Total horizontal distance is "<<TotalHorizontalDistance<<" "<<GetHitPar[0]<<std::endl;
    }

    //std::cout<<"Total horizontal distance is "<<TotalHorizontalDistance<<std::endl;
    
    ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
    double StraightLine_slope=tan(RayTracingFunctions::pi/2-LaunchAngleAir*(RayTracingFunctions::pi/180));
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
    struct RayTracingFunctions::fDnfR_params params2a;
    struct RayTracingFunctions::fDnfR_params params2b;

    ////Start looping over the layers to trace out the ray
    for(int il=0;il<RayTracingFunctions::MaxLayers-SkipLayersAbove-SkipLayersBelow;il++){
    
      if(il==0){
	////If this is the first layer then set the start height to be the height of the source
	LayerStartHeight=AirTxHeight;
      }else{
	////If this is any layer after the first layer then set the start height to be the starting height of the next layer or the end height of the previous layer
	LayerStartHeight=LastHeight-0.00001;
      }

      if(il==RayTracingFunctions::MaxLayers-SkipLayersAbove-SkipLayersBelow-1){
	////If this is the last layer then set the stopping height to be the height of the ice layer
	LayerStopHeight=IceLayerHeight;
      }else{
	////If this is NOT the last layer then set the stopping height to be the end height of the layer
	LayerStopHeight=(RayTracingFunctions::ATMLAY[RayTracingFunctions::MaxLayers-SkipLayersAbove-SkipLayersBelow-il-1]/100);
      }
    
      //std::cout<<il<<" A="<<layerAs[il]<<" ,B="<<layerBs[il]<<" ,C="<<layerCs[il]<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<" ,LayerStartHeight="<<LayerStartHeight<<" ,LayerStopHeight="<<LayerStopHeight<<std::endl;
      //std::cout<<" new layer "<<std::endl;
      ////Start tracing out the ray as it propagates through the layer
      for(double i=LayerStartHeight;i>LayerStopHeight-1;i=i-1){

	if(i<LayerStopHeight){
	  i=LayerStopHeight;
	}
	
	////Get and Set the A,B,C and L parameters for the layer
	params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-i), RayTracingFunctions::GetC_air(-i), layerLs[il]};
	params2b = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-(i)), RayTracingFunctions::GetC_air(-(i)), layerLs[il]};
      
	////Calculate the x (distance) value of the ray for given y (height) value
	Refracted_x=RayTracingFunctions::fDnfR(-i,&params2a)-RayTracingFunctions::fDnfR(-(LayerStartHeight),&params2b)+LastRefracted_x;
	
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
    struct RayTracingFunctions::fDnfR_params params3a;
    struct RayTracingFunctions::fDnfR_params params3b;
    
    for(int i=0;i>-(AntennaDepth+1);i=i-1){
      params3a = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(i), RayTracingFunctions::GetC_ice(i), Lvalue};
      params3b = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(0), RayTracingFunctions::GetC_ice(0), Lvalue};

      double refractedpath=LastRefracted_x-RayTracingFunctions::fDnfR((double)i,&params3a)+RayTracingFunctions::fDnfR(0,&params3b);
      aout<<ipoints<<" "<<refractedpath<<" "<<(double)i+IceLayerHeight<<std::endl;
      ipoints++;
    }

  } 
 
  t2b = std::chrono::high_resolution_clock::now();
  durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  durationb=durationb/1000;
  std::cout<<"total time taken by the script to store rays: "<<durationb<<" ms"<<std::endl;

  return 0;
}
