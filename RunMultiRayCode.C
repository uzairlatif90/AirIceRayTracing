#include "MultiRayAirIceRefraction.cc"

void RunMultiRayCode(){
  
  ////All variables are in m here
  double AntennaDepth=-200;////Depth of antenna in the ice
  double IceLayerHeight=3000;////Height where the ice layer starts off
  double AntennaNumber=0;
  double AirTxHeight=5000;
  double HorizontalDistance=1000;
  bool UseTable=true;
  
  double opticalPathLengthInIce;
  double opticalPathLengthInAir;
  double geometricalPathLengthInIce;
  double geometricalPathLengthInAir;
  double launchAngle;
  double horidist2interpnt;
  double transmissionCoefficientS;
  double transmissionCoefficientP;
  double recieveAngleinIce;
  
  bool CheckSol=false;////check if solution exists or not
  
  if(UseTable==true){
    
    //MultiRayAirIceRefraction::MakeTable(IceLayerHeight,AntennaDepth);
    MultiRayAirIceRefraction::MakeRayTracingTable(AntennaDepth*100,IceLayerHeight*100,AntennaNumber);

    CheckSol=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint_Table(AirTxHeight*100, HorizontalDistance*100 ,AntennaDepth*100, IceLayerHeight*100,0, opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt,transmissionCoefficientS,transmissionCoefficientP,recieveAngleinIce);
  }else{
    
    MultiRayAirIceRefraction::MakeAtmosphere();
    CheckSol=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(AirTxHeight*100, HorizontalDistance*100 ,AntennaDepth*100, IceLayerHeight*100,opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt,transmissionCoefficientS,transmissionCoefficientP,recieveAngleinIce);
    
  }

  if(CheckSol==true){
    cout<<" We have a solution!!!"<<endl;
    cout<<"AntennaNumber: "<< AntennaNumber<<endl;
    cout<<"AirTxHeight: "<<AirTxHeight<<endl;
    cout<<"HorizontalDistance: "<<HorizontalDistance<<endl;
    cout<<"opticalPathLengthInIce: "<<opticalPathLengthInIce/100<<endl;
    cout<<"opticalPathLengthInAir: "<<opticalPathLengthInAir/100<<endl;
    cout<<"launchAngle: "<<launchAngle*(180/MultiRayAirIceRefraction::pi)<<endl;
    cout<<"horidist2interpnt: "<<horidist2interpnt/100<<endl;
    cout<<"transmissionCoefficientS: "<<transmissionCoefficientS<<endl;
    cout<<"transmissionCoefficientP: "<<transmissionCoefficientP<<endl;
    cout<<"recieveAngleinIce: "<<recieveAngleinIce*(180/MultiRayAirIceRefraction::pi)<<endl;
  }else{
    cout<<" We do NOT have a solution!!!"<<endl;
  }
  
}
