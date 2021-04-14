#include "RayTracingFunctions.cc"

void SingleRayAirIceRefraction_wROOTGr(double AntennaDepth, double RayLaunchAngle, double AirTxHeight, double IceLayerHeight){
  //void SingleRayAirIceRefraction_wROOTGr(){

  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();
  
  // double AntennaDepth=200;////Depth of antenna in the ice
  // double RayLaunchAngle=124.765;////Initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 is vertically down
  // double AirTxHeight=5000;////Height of the source
  // double IceLayerHeight=3000;////Height where the ice layer starts off

  ////Fill in the n(h) and h arrays and ATMLAY and a,b and c (these 3 are the mass overburden parameters) from the data file
  RayTracingFunctions::MakeAtmosphere();
  
  if(AirTxHeight>h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]){
    ////Maximum height available with the refractive index data
    cout<<"Tx Height is set higher than maximum available height for atmospheric refractive index which is "<<h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1]<<endl;
    AirTxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];
    cout<<"Setting Tx Height to be the maximum available height"<<endl;
  }

  if(RayLaunchAngle<=90){
    cout<<"RayLaunchAngle has been set at "<<RayLaunchAngle<<" deg which is outside of the allowed range of 90 deg <RayLaunchAngle< 180 deg"<<endl;
    RayLaunchAngle=135;
    cout<<"Setting RayLaunchAngle at"<<RayLaunchAngle<<endl;
  }
  
  ////Print out the ray path x and y values in a file
  ofstream aout("RayPathinAirnIce.txt");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Section for propogating the ray through the atmosphere
  
  ////Find out how many atmosphere layers are above the source or Tx which we do not need
  int skiplayer=0;
  for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
    //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<endl;
    if(AirTxHeight<ATMLAY[ilayer]/100 && AirTxHeight>=ATMLAY[ilayer-1]/100){
      cout<<"Tx Height is in this layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<endl;
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
    if(IceLayerHeight>=ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
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
  double TotalHorizontalDistance=0;
  vector <double> layerLs;////vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them

  double RecieveAngle=0;
  double Lvalue=0;
  
  ////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
  cout<<"Fitting the atmosphere refrative index profile with multiple layers and propogate the ray"<<endl;
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
      StartAngle=180-RayLaunchAngle;
    }
    //cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<endl;

    ////Get the hit parameters from the function. The output is:
    //// How much horizontal distance did the ray travel in the layer
    //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
    //// The value of the L parameter for that layer
    if(ilayer==MaxLayers-SkipLayersAbove-1){
      double* GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
      TotalHorizontalDistance+=GetHitPar[0];
      RecieveAngle=GetHitPar[1];
      StartAngle=GetHitPar[1];
      ////Store in the values of A,B,C and L for tha layer
      Lvalue=GetHitPar[2];
      layerLs.push_back(GetHitPar[2]);
      delete []GetHitPar;  
    }
    if(ilayer<MaxLayers-SkipLayersAbove-1){
      double nzStopHeight=RayTracingFunctions::Getnz_air(StopHeight);
      double RecAng=asin(Lvalue/nzStopHeight);
      RecAng=RecAng*(180/RayTracingFunctions::pi);
      double THD=RayTracingFunctions::GetRayOpticalPath(RayTracingFunctions::A_air, StopHeight, StartHeight, Lvalue, 1);
      TotalHorizontalDistance+=THD;
      StartAngle=RecAng;
      ////Store in the values of A,B,C and L for tha layer
      layerLs.push_back(Lvalue);
    }
    
  }

  double IncidentAngleonIce=StartAngle;
  cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<endl;
  cout<<"Now propagating the ray through the ice towards the antenna"<<endl;
  
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
  //cout<<"Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingDepth="<<StartDepth<<" ,StoppingDepth="<<AntennaDepth<<" ,RayLaunchAngle="<<StartAngle<<endl;
  
  ////Get the hit parameters from the function. The output is:
  //// How much horizontal distance did the ray travel through ice to hit the antenna
  //// The angle of reciept/incidence at the end at the antenna
  //// The value of the L parameter for whole atmosphere fit
  double *GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, AntennaDepth,StartDepth, StartAngle, 0);

  ////SLF here stands for Single Layer Fitting. These variables store the hit parameters
  double TotalHorizontalDistanceinIce=GetHitPar[0];
  double RecievdAngleinIce=GetHitPar[1];
  double LvalueIce=GetHitPar[2];

  delete[] GetHitPar;
  
  cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceinIce<<endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice
  
  TMultiGraph *mg=new TMultiGraph();
  TMultiGraph *mgB=new TMultiGraph();
  
  TGraph *grRefracted=new TGraph();
  TGraph *grStraight=new TGraph();
  TGraph *grResidual=new TGraph();
  TGraph *grAirIce=new TGraph();
  TGraph *grIceLayer=new TGraph();

  ////Make a straight line at the same launch angle as the refracted ray in air to calculate the residual
  double StraightLine_slope=tan(RayTracingFunctions::pi/2-RayLaunchAngle*(RayTracingFunctions::pi/180));
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
      LayerStopHeight=IceLayerHeight;
    }else{
      ////If this is NOT the last layer then set the stopping height to be the end height of the layer
      LayerStopHeight=(ATMLAY[MaxLayers-SkipLayersAbove-SkipLayersBelow-il-1]/100)-1;
    }
    
    //std::cout<<il<<" A="<<layerAs[il]<<" ,B="<<layerBs[il]<<" ,C="<<layerCs[il]<<" ,L="<<layerLs[il]<<" , StartHeight="<<StartHeight<<" ,StopHeight="<<StopHeight<<" ,LayerStartHeight="<<LayerStartHeight<<" ,LayerStopHeight="<<LayerStopHeight<<std::endl;
    //std::cout<<" new layer "<<std::endl;
    ////Start tracing out the ray as it propagates through the layer
    for(double i=LayerStartHeight;i>LayerStopHeight-0.01;i=i-0.01){
      
      ////Get and Set the A,B,C and L parameters for the layer
      params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-i), RayTracingFunctions::GetC_air(-i), layerLs[il]};
      params2b = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-(i)), RayTracingFunctions::GetC_air(-(i)), layerLs[il]};
      
      ////Calculate the x (distance) value of the ray for given y (height) value
      Refracted_x=RayTracingFunctions::fDnfR(-i,&params2a)-fDnfR(-(LayerStartHeight),&params2b)+LastRefracted_x;
      
      ////If the ray just started off in a new layer we might need to offset the x values of the new layer so that the ray aligns with the previous layer.
      // if(ipoints>0 && fabs(i-LayerStartHeight-1)<2){
      //   LayerHoriOffset=(Refracted_x-LastRefracted_x);
      //   std::cout<<il<<" layer offset is "<<LayerHoriOffset<<std::endl;
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
  
  for(int i=0;i>-(AntennaDepth+1);i--){
    params3a = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(i), RayTracingFunctions::GetC_ice(i), LvalueIce};
    params3b = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(0), RayTracingFunctions::GetC_ice(0), LvalueIce};
    
    double refractedpath=LastRefracted_x-RayTracingFunctions::fDnfR((double)i,&params3a)+RayTracingFunctions::fDnfR(0,&params3b);
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

  mgB->SetTitle("RayPath through Air and Ice;Distance (m); Height (m)");
  grAirIce->SetMarkerStyle(20);
  grIceLayer->SetMarkerStyle(20);
  grIceLayer->SetMarkerColor(kBlue);
  mgB->Add(grAirIce);
  mgB->Add(grIceLayer);
  TCanvas *cb=new TCanvas("cb","cb");
  cb->cd();
  mgB->Draw("ALP");

  // delete accelerator;
  // delete spline;
  // flattened_h_data.clear();
  // flattened_nh_data.clear();

  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  durationb=durationb/1000;
  std::cout<<"total time taken by the script: "<<durationb<<" ms"<<std::endl;
}
