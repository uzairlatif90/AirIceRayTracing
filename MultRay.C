#include "IceRayTracing.cc"
#include "RayTracingFunctions.cc"

void MultRay(){

  RayTracingFunctions::MakeAtmosphere();
  
  //Plot the ray solutions
  Double_t x0=0;//Tx x positions
  Double_t z0=-180;//Tx z position
  Double_t x1=1000;//Tx z position
  //Double_t z1=z0-100;//Set the lower limit for the rays
  Double_t z1=z0;//Set the lower limit for the rays
  Double_t setzrange=-250;//Set the lower limit for plotting
  Double_t setxrange=800;//Set the distance limit for plotting
  Double_t launchangle=0;//variable defined for the for loop
  
  Int_t totray=90*2*2;//total number of rays
  Double_t launchinterval=0.25;//launch angle step size in degrees

  ////store ray paths
  TGraph *grR[totray];
  TGraph *grAir[totray];
  TMultiGraph *mg=new TMultiGraph();
  
  //for(Int_t iang=201;iang<202;iang++){
  for(Int_t iang=0;iang<totray;iang++){
    launchangle=iang*launchinterval;
    
    Double_t lvalueR=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t lvalueRa=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t zn=z1;

    /* This function returns the zmax for Reflected/Refracted ray path in a TGraph*/
    Double_t zmax=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa)+0.0000001;
    
    if(zmax>1e-5){
      zn=z1;
      /* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
      grR[iang]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,zmax,lvalueRa);
    }else{
      zn=z1;
      /* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
      grR[iang]=IceRayTracing::GetFullReflectedRayPath(z0,x1,z1,lvalueR);
    }
    grR[iang]->SetLineColor(kBlue);

    struct IceRayTracing::fDnfR_L_params params1b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), -z0};
    struct IceRayTracing::fDnfR_L_params params1c = {IceRayTracing::A_ice, IceRayTracing::GetB(0.0000001), -IceRayTracing::GetC(0.0000001), 0.0000001};
    double TotalHoriDist=fDnfR_L(lvalueR,&params1b) - fDnfR_L(lvalueR,&params1c);
    
    if(zmax<1e-5){
      double Refracted_x=0;////X coordinate variable which stores the horizontal distance of the refracted ray
      double LastRefracted_x=TotalHoriDist;////X coordinate variable which stores the horizontal distance of the refracted ray in the previous iteration
      double LayerStartHeight=0;////The starting height for the propagation in the layer
      double LayerStopHeight=0;////The stopping height for the propagation in the layer
      int ipoints=0;////variable for counting the total number of samples that make up the ray path
      Double_t lastB=0;
      Double_t lastC=0;
      
      ////Get and Set the A,B,C and L parameters for the layer
      struct RayTracingFunctions::fDnfR_params params2a;
      struct RayTracingFunctions::fDnfR_params params2b;
      params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(0), RayTracingFunctions::GetC_air(0), lvalueR};
      Refracted_x=RayTracingFunctions::fDnfR(0,&params2a);  
      if(isnan(Refracted_x)==false){
	//cout<<iang<<" "<<zmax<<" we are in "<<launchangle<<endl;
	grAir[iang]=new TGraph();
	grAir[iang]->SetLineColor(kRed);
	////Define ray variables for plotting and/or storing ray path as it comes down from the atmosphere

	LayerStartHeight=0;
	LayerStopHeight=50;
      
	////Start tracing out the ray as it propagates through the layer
	for(double i=0;i<LayerStopHeight+0.1;i=i+0.1){
	  ////Get and Set the A,B,C and L parameters for the layer
	  params2a = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(i), RayTracingFunctions::GetC_air(i), lvalueR};
	  params2b = {RayTracingFunctions::A_air, RayTracingFunctions::GetB_air(i), RayTracingFunctions::GetC_air(i), lvalueR};
	  
	  if(RayTracingFunctions::GetB_air(i)!=lastB && i>LayerStartHeight){
	    LayerStartHeight=i;
	    LastRefracted_x=Refracted_x;	
	  } 
	  ////Calculate the x (distance) value of the ray for given y (height) value
	  Refracted_x=RayTracingFunctions::fDnfR(i,&params2a)-RayTracingFunctions::fDnfR((LayerStartHeight),&params2b)+LastRefracted_x;  
	  grAir[iang]->SetPoint(ipoints,Refracted_x,i);
	  //cout<<i<<" "<<Refracted_x<<" "<<RayTracingFunctions::fDnfR(i,&params2a)<<" "<<RayTracingFunctions::fDnfR((LayerStartHeight),&params2b)<<endl;
	  ipoints++;
	  
	  lastB=RayTracingFunctions::GetB_air(i);
	  lastC=RayTracingFunctions::GetC_air(i);
	}
	mg->Add(grAir[iang]);
      }
    }
    
    mg->Add(grR[iang]);
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
