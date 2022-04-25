#include "RayTracingFunctions.cc"

void Air2IceRayTracing_wROOTplot(double AirTxHeight, double HorizontalDistance, double IceLayerHeight, double AntennaDepth){

  RayTracingFunctions::MakeAtmosphere();
  
  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();

  //cout<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<endl;
  
  // double AirTxHeight=5000;////Height of the source
  // double HorizontalDistance=1000;////Horizontal distance
  // double IceLayerHeight=3000;////Height where the ice layer starts off
  // double AntennaDepth=200;////Depth of antenna in the ice  
  bool PlotRayPath=true;
  
  if(AirTxHeight<IceLayerHeight){
    cout<<"WARNING: AirTxHeight is less than IceLayerHeight."<<endl;
    cout<<"Setting AirTxHeight to be 100 m above the IceLayerHeight"<<endl;
    AirTxHeight=IceLayerHeight+100;
  }
  
  gsl_function F1;
  struct RayTracingFunctions::MinforLAng_params params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  F1.function = &RayTracingFunctions::MinimizeforLaunchAngle;
  F1.params = &params1;

  ////Set the initial angle limits for the minimisation

  ////Set the initial angle limits for the minimisation
  double startanglelim=90;
  double endanglelim=180;

  //cout<<"max layers are "<<RayTracingFunctions::MaxLayers<<endl;
  
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
    
    if(isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0){
      checknan=true;
    }else{
      startanglelim=startanglelim+0.05;
    }
  }
  
  checknan=false;
  while(checknan==false && endanglelim<180.1){
    double *GetResultsAirTest2= RayTracingFunctions::GetAirPropagationPar(endanglelim,AirTxHeight,IceLayerHeight);
    TotalHorizontalDistanceinAirt=0;
    FilledLayerst=GetResultsAirTest2[4*RayTracingFunctions::MaxLayers];
    for(int i=0;i<FilledLayerst;i++){
      TotalHorizontalDistanceinAirt+=GetResultsAirTest2[i*4];
    }
    delete []GetResultsAirTest2;
    
    if(isnan(TotalHorizontalDistanceinAirt)==false && (TotalHorizontalDistanceinAirt)>0){
      checknan=true;
    }else{
      endanglelim=endanglelim-0.05;
    }
  }

  std::cout<<"Launch Angle search range is:  Startangle "<<startanglelim<<" ,Endangle "<<endanglelim<<std::endl;
  // //std::cout<<" "<<std::endl;
  
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleTx=RayTracingFunctions::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleTx<<" deg"<<endl;
 
  double * GetResultsAir=RayTracingFunctions::GetAirPropagationPar(LaunchAngleTx,AirTxHeight,IceLayerHeight);
 
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

  cout<<"***********Results for Air************"<<endl;
  cout<<setprecision(10)<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<endl;
  cout<<setprecision(10)<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<endl;
  cout<<"LvalueAir for "<<Lvalue<<endl;
  cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<endl;
  
  double * GetResultsIce=RayTracingFunctions::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth,Lvalue);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  double IncidentAngleonAntenna=GetResultsIce[1];
  //double LvalueIce=GetResultsIce[2];
  double PropagationTimeIce=GetResultsIce[3]*pow(10,9);
  delete [] GetResultsIce;

  cout<<" "<<endl;
  cout<<"***********Results for Ice************"<<endl;
  cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<endl;
  cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<endl;
  cout<<"LvalueIce "<<Lvalue<<endl;
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
    //ofstream aout("RayPathinAirnIce.txt");

    ////Find out how many atmosphere layers are above the source or Tx which we do not need
    int skiplayer=0;
    for(int ilayer=RayTracingFunctions::MaxLayers;ilayer>-1;ilayer--){
      //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<endl;
      if(AirTxHeight<RayTracingFunctions::ATMLAY[ilayer]/100 && AirTxHeight>=RayTracingFunctions::ATMLAY[ilayer-1]/100){
	//cout<<"Tx Height is in this layer with a height range of "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" m to "<<RayTracingFunctions::ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<endl;
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
    for(int ilayer=0;ilayer<RayTracingFunctions::MaxLayers;ilayer++){
      //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer+1]/100<<endl;
      if(IceLayerHeight>=RayTracingFunctions::ATMLAY[ilayer]/100 && IceLayerHeight<RayTracingFunctions::ATMLAY[ilayer+1]/100){
	//cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
	ilayer=100;
      }
      if(ilayer<RayTracingFunctions::MaxLayers){
	skiplayer++;
      }
    }
    int SkipLayersBelow=skiplayer;
    //cout<<"The total number of layers that need to be skipped from below is "<<skiplayer<<endl;
    
    TMultiGraph *mg=new TMultiGraph();
    TMultiGraph *mgB=new TMultiGraph();
  
    TGraph *grRefracted=new TGraph();
    TGraph *grStraight=new TGraph();
    TGraph *grResidual=new TGraph();
    TGraph *grIceLayer=new TGraph();
    TGraph *grAirIce=new TGraph();
    
    double Start_nh=0;
    double StartHeight=0;
    double StopHeight=0;
    double StartAngle=0;
    double TotalHorizontalDistance=0;
    vector <double> layerLs;////vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them
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
	StartAngle=180-LaunchAngleTx;
      }
      //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;

      //cout<<ilayer<<" StartingHeight="<<StartHeight<<" "<<gsl_spline_eval(spline, StartHeight, accelerator)<<" ,StoppingHeight= "<<StopHeight<<" "<<gsl_spline_eval(spline, StopHeight, accelerator)<<endl;
      
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
	//Lvalue=0.427404936;//GetHitPar[2];
	//GetHitPar[2]=Lvalue;
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
      //cout<<"start angle is "<<StartAngle<<" "<<GetHitPar[2]<<endl;
      
    }

    cout<<"Total horizontal distance is "<<TotalHorizontalDistance<<endl;
    
    ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
    double StraightLine_slope=tan(RayTracingFunctions::pi/2-LaunchAngleTx*(RayTracingFunctions::pi/180));
    double StraightLine_y_intercept=AirTxHeight;
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
    
      //cout<<il<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<" ,LayerStartHeight="<<LayerStartHeight<<" ,LayerStopHeight="<<LayerStopHeight<<endl;
      //cout<<" new layer "<<endl;
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
	//   cout<<il<<" layer offset is "<<LayerHoriOffset<<endl;
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
	//aout<<ipoints<<" "<<Refracted_x<<" "<<i<<endl;

	// ////If you want to check the the transition between different layers uncomment these lines
	// if(fabs(i-LayerStopHeight)<10){
	//   cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<" "<<GetB_air(-i)<<" "<<GetC_air(-i)<<" "<<layerLs[il]<<" "<<fDnfR(-i,&params2a)<<" "<<-fDnfR(-(LayerStartHeight),&params2b)<<" "<<LayerStartHeight<<" "<<LastRefracted_x<<" "<<GetB_air(-LayerStartHeight)<<" "<<GetC_air(-LayerStartHeight)<<endl;
	// }
      
	// if(fabs(i-LayerStartHeight)<10){
	//   cout<<"ipoints= "<<ipoints<<" ,ref_x="<<Refracted_x<<" ,i="<<i<<" "<<Refracted_x<<" "<<StraightLine_y<<" "<<StraightLine_y-i<<" "<<GetB_air(-i)<<" "<<GetC_air(-i)<<" "<<layerLs[il]<<" "<<fDnfR(-i,&params2a)<<" "<<-fDnfR(-(LayerStartHeight),&params2b)<<" "<<LayerStartHeight<<" "<<LastRefracted_x<<" "<<GetB_air(-LayerStartHeight)<<" "<<GetC_air(-LayerStartHeight)<<endl;
	// }
      
	ipoints++;
	LastHeight=i;
      }
      LastRefracted_x=Refracted_x;
    }    

    //cout<<"last ref values are "<<LastRefracted_x<<" "<<LastHeight<<endl;
    
    ////Print out the ray path in ice too  
    struct RayTracingFunctions::fDnfR_params params3a;
    struct RayTracingFunctions::fDnfR_params params3b;

    //for(int i=0;i>-(0.1);i=i-1){
    for(int i=0;i>-(AntennaDepth+1);i=i-1){
      params3a = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(i), RayTracingFunctions::GetC_ice(i), Lvalue};
      params3b = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(0), RayTracingFunctions::GetC_ice(0), Lvalue};

      double refractedpath=LastRefracted_x-RayTracingFunctions::fDnfR((double)i,&params3a)+RayTracingFunctions::fDnfR(0,&params3b);
      grAirIce->SetPoint(ipoints,refractedpath,(double)i+IceLayerHeight);
      grIceLayer->SetPoint(ipoints,refractedpath,IceLayerHeight);
      //cout<<ipoints<<" "<<refractedpath<<" "<<(double)i+IceLayerHeight<<" "<<i<<endl;
      ipoints++;
    }

    grRefracted->SetLineColor(2);//red
    grStraight->SetLineColor(4);//blue
    grRefracted->SetLineWidth(5);
    grStraight->SetLineWidth(2);
  
    grStraight->SetTitle("Straight Line with same launch angle");
    grRefracted->SetTitle("Refracted Ray");

    TString mgtitle="Launch Angle=";
    mgtitle+=LaunchAngleTx;
    mgtitle+=" deg where 0 deg is straight down;Distance (m); Height (m)";
    mg->SetTitle(mgtitle);
    grResidual->SetTitle("Difference of Red with Blue;Distance (m); Height (m)");
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

    // delete mg;
    // delete mgB;

    // delete grAirIce;
    // delete grRefracted;
    // delete grStraight;
    // delete grResidual;
    // delete grIceLayer;
    
  }
  
  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  durationb=durationb/1000;
  std::cout<<"total time taken by the script: "<<durationb<<" ms"<<std::endl;
}
