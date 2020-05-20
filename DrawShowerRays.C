#include "RayTracingFunctions.cc"

int Air2IceRayTracing_gr(double AirTxHeight, double HorizontalDistance, double IceLayerHeight, double AntennaDepth, TGraph *grAirIce, Double_t displacement){
  
  ////For recording how much time the process took
  timestamp_t t0 = get_timestamp();
  
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
  double startanglelim=91.5;
  double endanglelim=178.5;

  ////Start opening up the angle limit range until the air minimisation function becomes undefined or gives out a nan. Then set the limits within that range.
  bool checknan=false;
  double TotalHorizontalDistanceinAirt=0;
  int FilledLayerst=0;
  while(checknan==false){
    double *GetResultsAirTest1=RayTracingFunctions::GetAirPropagationPar(startanglelim,AirTxHeight,IceLayerHeight);
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
    }
  }
  checknan=false;
  while(checknan==false){
    double *GetResultsAirTest2=RayTracingFunctions::GetAirPropagationPar(endanglelim,AirTxHeight,IceLayerHeight);
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
  
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleTx=RayTracingFunctions::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  //cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleTx<<" deg"<<endl;
  
  double * GetResultsAir=RayTracingFunctions::GetAirPropagationPar(LaunchAngleTx,AirTxHeight,IceLayerHeight);
  int FilledLayers=GetResultsAir[4*MaxLayers];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[i*4];
    PropagationTimeAir+=GetResultsAir[3+i*4]*pow(10,9);
  }
  double Lvalue=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*4];  
  
  delete [] GetResultsAir;
  
  double * GetResultsIce=RayTracingFunctions::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth,Lvalue);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  double IncidentAngleonAntenna=GetResultsIce[1];
  //double LvalueIce=GetResultsIce[2];
  double PropagationTimeIce=GetResultsIce[3]*pow(10,9);
  delete [] GetResultsIce;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;
  
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////This section now is for storing and plotting the ray. We calculate the x and y coordinates of the ray as it travels through the air and ice
  if(PlotRayPath==true){
    ////Print out the ray path x and y values in a file
    //ofstream aout("RayPathinAirnIce.txt");

    ////Find out how many atmosphere layers are above the source or Tx which we do not need
    int skiplayer=0;
    for(int ilayer=MaxLayers;ilayer>-1;ilayer--){
      //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer-1]/100<<endl;
      if(AirTxHeight<ATMLAY[ilayer]/100 && AirTxHeight>=ATMLAY[ilayer-1]/100){
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
      //cout<<ilayer<<" "<<ATMLAY[ilayer]/100<<" "<<ATMLAY[ilayer+1]/100<<endl;
      if(IceLayerHeight>=ATMLAY[ilayer]/100 && IceLayerHeight<ATMLAY[ilayer+1]/100){
  	//cout<<"Ice Layer is in the layer with a height range of "<<ATMLAY[ilayer]/100<<" m to "<<ATMLAY[ilayer+1]/100<<" m and is at a height of "<<IceLayerHeight<<" m"<<endl;
  	ilayer=100;
      }
      if(ilayer<MaxLayers){
  	skiplayer++;
      }
    }
    int SkipLayersBelow=skiplayer;
    //cout<<"The total number of layers that need to be skipped from below is "<<skiplayer<<endl;
    
    double Start_nh=0;
    double StartHeight=0;
    double StopHeight=0;
    double StartAngle=0;
    double TotalHorizontalDistance=0;
    vector <double> layerLs;////vector for storing the A,B,C and L values of each of the atmosphere layers as the ray passes through them
    double RecieveAngle=0;
    double Lvalue=0;
    
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
  	StartAngle=180-LaunchAngleTx;
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
      //cout<<"start angle is "<<StartAngle<<" "<<GetHitPar[2]<<endl;
      
    }

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
    
      ////Start tracing out the ray as it propagates through the layer
      for(double i=LayerStartHeight;i>LayerStopHeight-0.01;i=i-0.01){
	
  	////Get and Set the A,B,C and L parameters for the layer
  	params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-i), RayTracingFunctions::GetC_air(-i), layerLs[il]};
  	params2b = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(-(i)), RayTracingFunctions::GetC_air(-(i)), layerLs[il]};
      
  	////Calculate the x (distance) value of the ray for given y (height) value
  	Refracted_x=RayTracingFunctions::fDnfR(-i,&params2a)-RayTracingFunctions::fDnfR(-(LayerStartHeight),&params2b)+LastRefracted_x;
	
  	////Caclulate the y value of the straight line
  	double StraightLine_y=StraightLine_slope*Refracted_x+StraightLine_y_intercept;

  	grAirIce->SetPoint(ipoints,displacement+Refracted_x,i);
      
  	ipoints++;
  	LastHeight=i;
      }
      LastRefracted_x=Refracted_x;
    }
    
    ////Print out the ray path in ice too  
    struct RayTracingFunctions::fDnfR_params params3a;
    struct RayTracingFunctions::fDnfR_params params3b;
    
    for(int i=0;i>-(AntennaDepth+1);i--){
      params3a = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(i), RayTracingFunctions::GetC_ice(i), Lvalue};
      params3b = {RayTracingFunctions::A_ice, RayTracingFunctions::GetB_ice(0), RayTracingFunctions::GetC_ice(0), Lvalue};

      double refractedpath=LastRefracted_x-RayTracingFunctions::fDnfR((double)i,&params3a)+RayTracingFunctions::fDnfR(0,&params3b);
      grAirIce->SetPoint(ipoints,displacement+refractedpath,(double)i+IceLayerHeight);
      ipoints++;
    }
    
  }
  
  timestamp_t t1 = get_timestamp();  
  double secs = (t1 - t0) / 1000000.0L;
  //cout<<"total time taken by the script: "<<secs<<" s"<<endl;

  return 0;
}

void DrawShowerRays(){
  RayTracingFunctions::MakeAtmosphere();
  
  double pi=4*atan(1.0);
  
  double AntennaDepth=180;////Depth of antenna in the ice
  double IceLayerHeight=2800;////Height where the ice layer starts off
  
  double initialheight=500+IceLayerHeight;
  double showerangle=30;
  double diststep=10;

  double initialhdist=(initialheight-IceLayerHeight)*tan(showerangle*(pi/180));
  initialhdist=initialhdist+100;
  cout<<"initial h distance is "<<initialhdist<<", 100 m from impact point"<<endl; 
    
  int totalpositions=((initialheight-IceLayerHeight)/cos(showerangle*(pi/180)))/diststep;
  totalpositions=100;
  //cout<<"total positions are "<<totalpositions<<endl;

  TMultiGraph *mg=new TMultiGraph();
  TGraph *gr1=new TGraph();
  TGraph *gr2=new TGraph();
  TGraph *gr3[totalpositions];
  for(int ipos=0;ipos<totalpositions;ipos++){
    gr3[ipos]=new TGraph();
  }
  
  Double_t layerdist=0;
  Int_t ilayer=0;
  while(layerdist<initialhdist){
    gr1->SetPoint(ilayer,layerdist,IceLayerHeight);
    layerdist+=50;
    ilayer++;
  }
  gr1->SetPoint(ilayer,layerdist,IceLayerHeight);
  
  int ipos=0;
    
  double rotx=0;
  double roty=IceLayerHeight+11;
  double firstx=0;
  
  double hdist=0;
  double x=0;
  double y=0;

  while(roty>=IceLayerHeight+10){
    x=0;
    y=initialheight-diststep*ipos;

    TVector3 v1(x,y-(((initialheight-IceLayerHeight)/2)+IceLayerHeight),0);
    v1.RotateZ(showerangle*(pi/180));
    rotx=v1.X()-firstx;
    roty=v1.Y()+(((initialheight-IceLayerHeight)/2)+IceLayerHeight);

    if(ipos==0){
      firstx=rotx;
      rotx=0;
      initialhdist=(roty-IceLayerHeight)*tan(showerangle*(pi/180));
      initialhdist=initialhdist+100;
      hdist=initialhdist;
      //cout<<"initial h distance is "<<initialhdist<<endl; 
    }
    
    if(roty>=10+IceLayerHeight){
      // if(y>=IceLayerHeight){
      // 	gr1->SetPoint(ipos,x,y);
      // }
      gr2->SetPoint(ipos,rotx,roty);
      //cout<<ipos<<" "<<rotx<<" "<<hdist-rotx<<" "<<roty<<endl;
      Air2IceRayTracing_gr(roty, hdist-rotx,IceLayerHeight, AntennaDepth,gr3[ipos],rotx);
      gr3[ipos]->SetMarkerColor(ipos+30);
      gr3[ipos]->SetLineColor(ipos+30);
      ipos++;
    }
  }
  
  //gr1->SetMarkerStyle(20);
  gr1->SetLineWidth(2);
  gr2->SetMarkerStyle(20);

  gr1->SetMarkerColor(kRed);
  gr2->SetMarkerColor(kBlue);
  
  mg->Add(gr1);
  mg->Add(gr2);

  cout<<"total positions are "<<totalpositions<<" "<<ipos<<endl; 
  for(int ip=0;ip<ipos;ip++){
    mg->Add(gr3[ip]);
  }
  mg->SetTitle("Shower Zenith=30 deg, Antenna 100 m away from IP, 180 m deep ;Distance (m) ;Altitude (m)");
  TCanvas *c1=new TCanvas("c1","c1");
  c1->cd();
  // mg->GetYaxis()->SetRangeUser(IceLayerHeight-AntennaDepth,IceLayerHeight+AntennaDepth);
  // mg->GetXaxis()->SetRangeUser(hdist-10,hdist+2);  
  // mg->GetYaxis()->SetRangeUser(IceLayerHeight-AntennaDepth,IceLayerHeight+AntennaDepth);
  // mg->GetXaxis()->SetRangeUser(hdist-150,hdist+2);  

  mg->Draw("ALP");  
}
