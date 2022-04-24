#include "IceRayTracing.cc"
#include "RayTracingFunctions.cc"

void MakeMultiRayPlot(Double_t z0, Double_t LaunchInterval){

  ////Make the atmosphere
  RayTracingFunctions::MakeAtmosphere();
  
  ////Plot the ray solutions
  //Double_t z0=-180;//Tx z position in meter
  Double_t ZLowerLimit=z0;//Set the lower limit for the rays
  //Double_t LaunchInterval=0.25;//launch angle step size in degrees
  Int_t TotalRays=90.0/LaunchInterval;//total number of rays
  Double_t MaxRayHeightInAir=50;//max ray height in air in meter
  
  ////store ray paths
  TGraph *grIce[TotalRays];
  TGraph *grAir[TotalRays];
  TMultiGraph *mg=new TMultiGraph();

  Double_t DummyVariable=0;
  Double_t LaunchAngle=0;//variable defined for the for-loop
  
  for(Int_t iang=0;iang<TotalRays;iang++){
    LaunchAngle=iang*LaunchInterval;
    
    Double_t LvalueR=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(LaunchAngle*(IceRayTracing::pi/180));
    Double_t LvalueRa=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(LaunchAngle*(IceRayTracing::pi/180));
    Double_t zn=ZLowerLimit;

    /* This function returns the zmax for Reflected/Refracted ray path in a TGraph*/
    Double_t zmax=IceRayTracing::GetZmax(IceRayTracing::A_ice,LvalueRa)+0.0000001;
    
    if(zmax>1e-5){
      zn=ZLowerLimit;
      /* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
      grIce[iang]=IceRayTracing::GetFullRefractedRayPath(z0,DummyVariable,ZLowerLimit,zmax,LvalueRa,0);
    }else{
      zn=ZLowerLimit;
      /* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
      grIce[iang]=IceRayTracing::GetFullReflectedRayPath(z0,DummyVariable,ZLowerLimit,LvalueR);
    }
    grIce[iang]->SetLineColor(kBlue);

    struct IceRayTracing::fDnfR_L_params params1b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), -z0};
    struct IceRayTracing::fDnfR_L_params params1c = {IceRayTracing::A_ice, IceRayTracing::GetB(0.0000001), -IceRayTracing::GetC(0.0000001), 0.0000001};
    double TotalHoriDist=fDnfR_L(LvalueR,&params1b) - fDnfR_L(LvalueR,&params1c);

    ////If the in-ice ray got reflected then start tracing the ray in air
    if(zmax<1e-5){
      double Refracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray
      double LastRefracted_x=TotalHoriDist;////X coordinate variable which stores the horizontal distance of the refracted ray in the previous iteration
      double LayerStartHeight=0;////The starting height for the propagation in the layer
      double LayerStopHeight=0;////The stopping height for the propagation in the layer
      int ipoints=0;////variable for counting the total number of samples that make up the ray path
      Double_t LastB=0;
      Double_t LastC=0;
      
      ////Get and Set the A,B,C and L parameters for the layer
      struct RayTracingFunctions::fDnfR_params params2a;
      struct RayTracingFunctions::fDnfR_params params2b;
      params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(0), RayTracingFunctions::GetC_air(0), LvalueR};
      Refracted_x=RayTracingFunctions::fDnfR(0,&params2a);  
      if(isnan(Refracted_x)==false){
	grAir[iang]=new TGraph();
	grAir[iang]->SetLineColor(kRed);
	////Define ray variables for plotting and/or storing ray path as it comes down from the atmosphere

	LayerStartHeight=0;
	LayerStopHeight=MaxRayHeightInAir;
      
	////Start tracing out the ray as it propagates through the layer
	for(double i=0;i<LayerStopHeight+0.1;i=i+0.1){
	  ////Get and Set the A,B,C and L parameters for the layer
	  params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(i), RayTracingFunctions::GetC_air(i), LvalueR};
	  params2b = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(i), RayTracingFunctions::GetC_air(i), LvalueR};
	  
	  if(RayTracingFunctions::GetB_air(i)!=LastB && i>LayerStartHeight){
	    LayerStartHeight=i;
	    LastRefracted_x=Refracted_x;	
	  } 
	  ////Calculate the x (distance) value of the ray for given y (height) value
	  Refracted_x=RayTracingFunctions::fDnfR(i,&params2a)-RayTracingFunctions::fDnfR((LayerStartHeight),&params2b)+LastRefracted_x;  
	  grAir[iang]->SetPoint(ipoints,Refracted_x,i);
	  
	  ipoints++;
	  
	  LastB=RayTracingFunctions::GetB_air(i);
	  LastC=RayTracingFunctions::GetC_air(i);
	}
	mg->Add(grAir[iang]);

	/* Setup the function that will be used to calculate the angle of reception for all the rays */
	gsl_function F5;
	double result, abserr;
 
	/* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
	struct IceRayTracing::fDnfR_params paramsIAngB = {IceRayTracing::A_ice, IceRayTracing::GetB(0.0000001), IceRayTracing::GetC(0.0000001), LvalueR};
	F5.function = &IceRayTracing::fDnfR; 
	F5.params = &paramsIAngB;
	gsl_deriv_central (&F5, -0.0000001, 1e-8, &result, &abserr);
	double IncidenceAngleInIce=atan(result)*(180.0/IceRayTracing::pi);
	
      }
    }
    
    mg->Add(grIce[iang]);
  }//iang loop

  TString title="Depth vs Distance, Tx Depth=";
  title+=z0;
  title+=" m; Distance (m);Depth (m)";
  mg->SetTitle(title); 
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->cd();
  mg->Draw("AL");
  mg->GetXaxis()->SetNdivisions(20);
  c1->SetGridx();
  
}
