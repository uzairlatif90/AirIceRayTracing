#include "MultiRayAirIceRefraction.cc"

std::vector <double> AntennaDepths;
std::vector <int> AntennaTableAlreadyMade;

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

    AntennaDepths.push_back(AntennaDepth*100);
    
    for(int i=0;i<AntennaDepths.size();i++){
    if(i==0){
      cout<<"\n Making table for Antenna "<<i<<" at "<< AntennaDepths[i]<<" cm "<<endl;
      MultiRayAirIceRefraction::MakeRayTracingTable(AntennaDepths[i], IceLayerHeight*100, i); 
      AntennaTableAlreadyMade.push_back(i);
    }

    bool MakeTable=true;
    for(int j=0;j<AntennaTableAlreadyMade.size();j++){
      if(AntennaDepths[i]==AntennaDepths[AntennaTableAlreadyMade[j]]){
	MakeTable=false;
      }
    }
     
    if(MakeTable==true){
      cout<<"\n Making table for Antenna "<<i<<" at "<< AntennaDepths[i] <<" cm "<<endl;
      MultiRayAirIceRefraction::MakeRayTracingTable(AntennaDepths[i], IceLayerHeight*100, i); 
      AntennaTableAlreadyMade.push_back(i);
    }
  }
    
   
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
