#include "AirIceRayTracing.cc"
/* 
g++ -O -g -Wall -I. -pthread -std=c++17 -m64 -fPIC -shared -o libAirIceRayTracing.so TraceIceToAir.C -lgsl -lgslcblas
 */
void TraceIceToAir(  double AntennaDepth, double IceLayerHeight,double AirTxHeight, double HorizontalDistance, double ArrayParameters[10]){
  
  // ////All variables are in m here
  // double AntennaDepth=-100;////Depth of antenna in the ice
  // double IceLayerHeight=3000;////Height where the ice layer starts off
  // double AntennaNumber=0;
  // double AirTxHeight=3000+200;
  // double HorizontalDistance=100;
  
  double opticalPathLengthInIce;
  double opticalPathLengthInAir;
  double geometricalPathLengthInIce;
  double geometricalPathLengthInAir;
  double launchAngle;
  double horidist2interpnt;
  double AngleOfIncidenceOnIce;
  double receivedAngle;
  
  bool CheckSol=false;////check if solution exists or not
  
  AirIceRayTracing::MakeAtmosphere("Atmosphere.dat");

  // AirIceRayTracing::A_const=AirIceRayTracing::Getnz_air(IceLayerHeight);
  // AirIceRayTracing::UseConstantRefractiveIndex=true;
  // AirIceRayTracing::A_air=AirIceRayTracing::A_const;
  
  CheckSol=AirIceRayTracing::GetRayTracingSolution(AirTxHeight, HorizontalDistance ,AntennaDepth, IceLayerHeight,opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt, AngleOfIncidenceOnIce, receivedAngle);

  std::swap(launchAngle,receivedAngle);
  receivedAngle=180-receivedAngle;
  if(CheckSol==true){
    std::cout<<" We have a solution!!!"<<std::endl;
    std::cout<<"AirTxHeight: "<<AirTxHeight<<std::endl;
    std::cout<<"HorizontalDistance: "<<HorizontalDistance<<std::endl;
    std::cout<<"geometricalPathLengthInIce: "<<geometricalPathLengthInIce<<std::endl;
    std::cout<<"geometricalPathLengthInAir: "<<geometricalPathLengthInAir<<std::endl;
    std::cout<<"launchAngle: "<<launchAngle<<std::endl;
    std::cout<<"RecievedAngle: "<<receivedAngle<<std::endl;
    std::cout<<"horidist2interpnt: "<<horidist2interpnt<<std::endl;
    std::cout<<"AngleOfIncidenceOnIce: "<<AngleOfIncidenceOnIce<<std::endl;

    ArrayParameters[0]=AirTxHeight;
    ArrayParameters[1]=HorizontalDistance;
    ArrayParameters[2]=geometricalPathLengthInIce;
    ArrayParameters[3]=geometricalPathLengthInAir;
    ArrayParameters[4]=launchAngle;
    ArrayParameters[5]=receivedAngle;
    ArrayParameters[6]=horidist2interpnt;
    ArrayParameters[7]=AngleOfIncidenceOnIce;
    ArrayParameters[8]=0;
    ArrayParameters[9]=0;

  }else{
    std::cout<<" We do NOT have a solution!!!"<<std::endl;
    ArrayParameters[0]=-1000;
    ArrayParameters[1]=-1000;
    ArrayParameters[2]=-1000;
    ArrayParameters[3]=-1000;
    ArrayParameters[4]=-1000;
    ArrayParameters[5]=-1000;
    ArrayParameters[6]=-1000;
    ArrayParameters[7]=-1000;
    ArrayParameters[8]=-1000;
    ArrayParameters[9]=-1000;

  }

  
}

extern "C" {
  void Py_TraceIceToAir(  double AntennaDepth, double IceLayerHeight,double AirTxHeight, double HorizontalDistance, double ArrayParameters[10]){
    return TraceIceToAir( AntennaDepth, IceLayerHeight, AirTxHeight, HorizontalDistance, ArrayParameters);
    }
}
