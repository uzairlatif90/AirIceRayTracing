#include "MultiRayAirIceRefraction.h"
#include "RayTracingFunctions.cc"

////Store the maximum possible height allowed by GDAS tables
double MaxAirTxHeight=0;
////Store the minimum possible height allowed antenna tables
double MinAirTxHeight=0;

////Multidimensional vector to store data read in from the antenna tables for interpolation puprposes
std::vector<std::vector<std::vector <double>>> AllTableAllAntData;

////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
double AngleStepSize=0.5;
double LoopStartAngle=92;
double LoopStopAngle=180;
int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize)+1;

////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
double HeightStepSize=20;
double LoopStartHeight=0;
double LoopStopHeight=0;
int TotalHeightSteps=0;

////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){
  
  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice
  // double AirTxHeight=5000;////Height of the source
  // double HorizontalDistance=1000;////Horizontal distance
  // double IceLayerHeight=3000;////Height where the ice layer starts off
  // double AntennaDepth=200;////Depth of antenna in the ice  

  //cout<<"Minimser values "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<endl;
  
  gsl_function F1;
  struct RayTracingFunctions::MinforLAng_params params1 = { AirTxHeight, IceLayerHeight, AntennaDepth, HorizontalDistance};
  F1.function = &RayTracingFunctions::MinimizeforLaunchAngle;
  F1.params = &params1;

  ////Set the initial angle limits for the minimisation
  double startanglelim=90;
  double endanglelim=180;

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
  
  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double LaunchAngleAir=RayTracingFunctions::FindFunctionRoot(F1,startanglelim,endanglelim,gsl_root_fsolver_brent,0.000000001);
  
  //std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  double * GetResultsAir=RayTracingFunctions::GetAirPropagationPar(LaunchAngleAir,AirTxHeight,IceLayerHeight);
  int FilledLayers=GetResultsAir[4*RayTracingFunctions::MaxLayers];
  double TotalHorizontalDistanceinAir=0;
  double PropagationTimeAir=0;
  double LvalueAir;
  for(int i=0;i<FilledLayers;i++){
    TotalHorizontalDistanceinAir+=GetResultsAir[i*4];
    PropagationTimeAir+=GetResultsAir[3+i*4];
  }
  LvalueAir=GetResultsAir[2];
  double IncidentAngleonIce=GetResultsAir[1+(FilledLayers-1)*4];    
  delete [] GetResultsAir;
  
  // std::cout<<"Result from the minimization: Air Launch Angle: "<<LaunchAngleAir<<" deg"<<std::endl;
  // std::cout<<"***********Results for Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinAir "<<TotalHorizontalDistanceinAir<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonIce "<<IncidentAngleonIce<<" deg"<<std::endl;
  //   std::cout<<"LvalueAir "<<LvalueAir<<std::endl;
  // std::cout<<"PropagationTimeAir "<<PropagationTimeAir<<" ns"<<std::endl;


  double * GetResultsIce=RayTracingFunctions::GetIcePropagationPar(IncidentAngleonIce, IceLayerHeight, AntennaDepth,LvalueAir);
  double TotalHorizontalDistanceinIce=GetResultsIce[0];
  //double IncidentAngleonAntenna=GetResultsIce[1];
  //double LvalueIce=GetResultsIce[2];
  double PropagationTimeIce=GetResultsIce[3];
  delete [] GetResultsIce;

  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice************"<<std::endl;
  // std::cout<<"TotalHorizontalDistanceinIce "<<TotalHorizontalDistanceinIce<<" m"<<std::endl;
  // std::cout<<"IncidentAngleonAntenna "<<IncidentAngleonAntenna<<" deg"<<std::endl;
  // std::cout<<"LvalueIce "<<LvalueIce<<std::endl;
  // std::cout<<"PropagationTimeIce "<<PropagationTimeIce*pow(10,9)<<" ns"<<std::endl;

  double TotalHorizontalDistance=TotalHorizontalDistanceinIce+TotalHorizontalDistanceinAir;
  //double TotalPropagationTime=PropagationTimeIce+PropagationTimeAir;

  // std::cout<<" "<<std::endl;
  // std::cout<<"***********Results for Ice + Air************"<<std::endl;
  // std::cout<<"TotalHorizontalDistance "<<TotalHorizontalDistance<<" m"<<std::endl;
  // std::cout<<"TotalPropagationTime "<<TotalPropagationTime*pow(10,9)<<" ns"<<std::endl;

  opticalPathLengthInIce=(PropagationTimeIce*RayTracingFunctions::spedc)*100;
  opticalPathLengthInAir=(PropagationTimeAir*RayTracingFunctions::spedc)*100;
  launchAngle=(LaunchAngleAir)*(RayTracingFunctions::pi/180);
  horizontalDistanceToIntersectionPoint=TotalHorizontalDistanceinAir*100;
  
  reflectionCoefficientS=RayTracingFunctions::Refl_S(IncidentAngleonIce*(RayTracingFunctions::pi/180),IceLayerHeight);
  reflectionCoefficientP=RayTracingFunctions::Refl_P(IncidentAngleonIce*(RayTracingFunctions::pi/180),IceLayerHeight);
  
  //cout<<" in minimser "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<endl;
  
  bool CheckSolution=false;
  double checkminimisation=TotalHorizontalDistance-HorizontalDistance;
  if(fabs(checkminimisation)<1){
    CheckSolution=true;
  }
  
  return CheckSolution;

}

////Just a simple function for interpolating in 1 dimension between two points (xa,ya) and (xb,yb)
double MultiRayAirIceRefraction::oneDLinearInterpolation(double x, double xa, double ya, double xb, double yb){
  double y=ya+(yb-ya)*((x-xa)/(xb-xa));
  return y;
}

////Just a simple function for extrapolating in 1 dimension outside two points
double MultiRayAirIceRefraction::Extrapolate(int Par, int index, double TotalHorizontalDistance, int AntennaNumber){
  double x1,x2,y1,y2,x3,m,c,y3;

  x1=AllTableAllAntData[AntennaNumber][1][index];
  x2=AllTableAllAntData[AntennaNumber][1][index+1];
  y1=AllTableAllAntData[AntennaNumber][Par][index];
  y2=AllTableAllAntData[AntennaNumber][Par][index+1];
  x3=TotalHorizontalDistance;
  
  m=(y2-y1)/(x2-x1);
  c=y1-m*x1;
  y3=m*x3+c;
  
  return y3;
}

////Find the limit to extrapolate for all the parameters using the Air Launch Angle parameter
double MultiRayAirIceRefraction::FindExtrapolationLimit( int index, double TotalHorizontalDistance, int AntennaNumber){
  double x1,x2,y1,y2,m,c,xlimit;
  //double y3,x3;
  
  x1=(AllTableAllAntData[AntennaNumber][1][index]);
  x2=(AllTableAllAntData[AntennaNumber][1][index+1]);
  y1=AllTableAllAntData[AntennaNumber][4][index];
  y2=AllTableAllAntData[AntennaNumber][4][index+1];
  //x3=(TotalHorizontalDistance);
  
  m=(y1-y2)/(x1-x2);
  c=y1-m*x1;
  //y3=m*x3+c;
  xlimit=(90-c)/m;
  
  return xlimit;
}

int MultiRayAirIceRefraction::FindClosestAirTxHeight(double IceLayerHeight, double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber){

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;  
  int CurrentHeightStep=floor((ParValue-LoopStopHeight)/HeightStepSize);

  int mainindex=(TotalHeightSteps-CurrentHeightStep-1);
  int index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
  index1=index1+mainindex;
  if(AllTableAllAntData[AntennaNumber][0][index1]<ParValue){
    mainindex=mainindex-1;
    index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
    index1=index1+mainindex;
  }
  int index2=((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))+index1;

  if(index2==TotalTableEntries-1){
   index1=index1-((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))-1;
   index2=index2-((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))-1;
  }
  
  if(ParValue==(AllTableAllAntData[AntennaNumber][0][index1]+AllTableAllAntData[AntennaNumber][0][index1-1])/2 && index1>0){
    mainindex=(TotalHeightSteps-CurrentHeightStep-1)-1;
    index1=(mainindex*((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)));
    index1=index1+mainindex;
    index2=((LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize))+index1;
  }
  
  RStartIndex=index1;
  REndIndex=index2;
  ClosestVal=fabs(AllTableAllAntData[AntennaNumber][0][index1]-ParValue);

  return 0;
}

int MultiRayAirIceRefraction::FindClosestTHD(double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber){

  double minimum=100000000000;
  int index2=0;
  double minval=0;
  //double lastminval=0;
  for(int ipnt=StartIndex;ipnt<EndIndex+1;ipnt++){
    minval=fabs(AllTableAllAntData[AntennaNumber][1][ipnt]-ParValue);
    if(minval<minimum && AllTableAllAntData[AntennaNumber][1][ipnt]>ParValue){
      minimum=minval;
    }else{
      index2=ipnt; 
      break;
    }
    //lastminval=minval;
  }
  double index1=index2-1;

  minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index2]);
  if(minimum>fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1])){
    minimum=fabs(ParValue-AllTableAllAntData[AntennaNumber][1][index1]);
  }
  
  RStartIndex=index1;
  REndIndex=index2;
  ClosestVal=minimum;
  
  //cout<<"start index is "<<index1<<" last index is "<<index2<<" "<<minimum<<endl;
  return 0;
}

////Interpolate the value of the given parameter for a given TxHeight and THD
int MultiRayAirIceRefraction::GetParValues(double AntennaNumber, double AirTxHeight, double TotalHorizontalDistance, double IceLayerHeight, double &AirTxHeight1, double Par1[6],double &AirTxHeight2, double Par2[6]){

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
  MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
  MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];
  
  int startindex[6];
  int endindex[6];
  double closestvalue[6];
  double x1,x2,y1,y2;
  double MaxTotalHorizontalDistance=0;
  double MinTotalHorizontalDistance=0;
  
  ////Find out the AirTxheight1 bin range for the given Tx height
  FindClosestAirTxHeight(IceLayerHeight,AirTxHeight,0, TotalTableEntries, startindex[0], endindex[0], closestvalue[0], AntennaNumber);
  ////Set the first Tx height
  AirTxHeight1=AllTableAllAntData[AntennaNumber][0][startindex[0]];

  ////Find the maximum and minimum possible values of THD for AirTxHeight1. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
  MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[0]];
  MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[0]];
  //cout<<" minimum distance is 1 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[0]+1]<<endl;
  if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){
    ////Find out the THD bin range for the given THD and AirTxHeight1
    FindClosestTHD(TotalHorizontalDistance, startindex[0], endindex[0], startindex[1], endindex[1], closestvalue[1], AntennaNumber);
    
    if(closestvalue[1]!=0){////If the given THD does not match a bin value
      x1=AllTableAllAntData[AntennaNumber][1][startindex[1]];
      x2=AllTableAllAntData[AntennaNumber][1][endindex[1]];
      for(int ipar=0;ipar<6;ipar++){
	y1=AllTableAllAntData[AntennaNumber][2+ipar][startindex[1]];
	y2=AllTableAllAntData[AntennaNumber][2+ipar][endindex[1]];
	Par1[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
      }
    }
    if(closestvalue[1]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the same and there is no interpolation needed
      startindex[1]=startindex[1]+1;
      endindex[1]=startindex[1];
      for(int ipar=0;ipar<6;ipar++){
	Par1[ipar]=AllTableAllAntData[AntennaNumber][2+ipar][startindex[1]];
      }
    }
    
  }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
    ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
    double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);

    if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
      double THDLimit=FindExtrapolationLimit(startindex[0],TotalHorizontalDistance,AntennaNumber);
      double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
      double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
      if(DistancePercentage<DistancePercentageLimit){
	for(int ipar=0;ipar<6;ipar++){
	  Par1[ipar]=Extrapolate(2+ipar, startindex[0], TotalHorizontalDistance, AntennaNumber);
	}
      }else{
	for(int ipar=0;ipar<6;ipar++){
	  Par1[ipar]=-pow(10,9);
	}
      }
    }
    
  }
  
  ////Find out the Txheight2 bin range for the given Tx height
  startindex[2]=startindex[0]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize)+1;
  endindex[2]=startindex[2]+(LoopStopAngle-LoopStartAngle)*(1.0/AngleStepSize);
    
  if(closestvalue[0]!=0 && AirTxHeight>MinAirTxHeight && startindex[2]<TotalTableEntries){////If the AirTxHeight does not exactly match a table value and AirTxHeight1 is not equal to AirTxHeight2 AND AirTxHeight1 is not the minimum possible height in the table

    ////Set the second Tx height
    AirTxHeight2=AllTableAllAntData[AntennaNumber][0][startindex[2]];
    ////Find the maximum and minimum possible values of THD for AirTxHeight2. If given THD from the user is in this range then do interpolation otherwise we have to do extrapolation
    MaxTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][startindex[2]];
    MinTotalHorizontalDistance=AllTableAllAntData[AntennaNumber][1][endindex[2]];
    //cout<<" minimum distance is 2 :"<<MinTotalHorizontalDistance<<" "<<AllTableAllAntData[AntennaNumber][1][endindex[2]+1]<<endl;
    if(TotalHorizontalDistance<=MaxTotalHorizontalDistance){

      FindClosestTHD(TotalHorizontalDistance, startindex[2], endindex[2], startindex[3], endindex[3], closestvalue[2], AntennaNumber);
       
      if(closestvalue[2]!=0){////If the given THD does not match a bin value
	x1=AllTableAllAntData[AntennaNumber][1][startindex[3]];
	x2=AllTableAllAntData[AntennaNumber][1][endindex[3]];
	for(int ipar=0;ipar<6;ipar++){
	  y1=AllTableAllAntData[AntennaNumber][2+ipar][startindex[3]];
	  y2=AllTableAllAntData[AntennaNumber][2+ipar][endindex[3]];
	  Par2[ipar]=oneDLinearInterpolation(TotalHorizontalDistance,x1,y1,x2,y2);
	}
      }
      if(closestvalue[2]==0){////if the value of THD matches exactly the value of one the values in the table then go into this if condition and set the the start and stop indexes to be the sohame and there is no interpolation needed
	startindex[3]=startindex[3]+1;
	endindex[3]=startindex[3];
	for(int ipar=0;ipar<6;ipar++){
	  Par2[ipar]=AllTableAllAntData[AntennaNumber][2+ipar][startindex[3]];
	}
      }

    }else{///Do extrapolation as the given THD from the user does NOT lie in the range of values
      ////check if the total horizontal distance is outside 100% of the whole total horizontal distance range for that Tx height. If that is the case then do not extrapolate and return a no solution number
      double DistanceRange=fabs(MaxTotalHorizontalDistance-MinTotalHorizontalDistance);
      if(TotalHorizontalDistance>MaxTotalHorizontalDistance){
	double THDLimit=FindExtrapolationLimit(startindex[2],TotalHorizontalDistance,AntennaNumber);
	double DistancePercentageLimit=fabs(THDLimit-MaxTotalHorizontalDistance)/DistanceRange;
	double DistancePercentage=fabs(TotalHorizontalDistance-MaxTotalHorizontalDistance)/DistanceRange;
	if(DistancePercentage<DistancePercentageLimit){
	  for(int ipar=0;ipar<6;ipar++){
	    Par2[ipar]=Extrapolate(2+ipar, startindex[2], TotalHorizontalDistance, AntennaNumber);
	  }
	}else{
	  for(int ipar=0;ipar<6;ipar++){
	    Par2[ipar]=-pow(10,9);
	  }
	}
      }
      
    }
    
  }else{////If AirTxHeight is exactly equal to a bin value then AirTxHeight1 and AirTxHeight2 are equal and we dont need to do interpolation for AirTxHeight2 as we already have done it for AirTxHeight1
    ////Also if AirTxHeight1 is equal to the minimum possible height then stop right there and dont the interpolation for AirTxHeight2
    AirTxHeight2=AirTxHeight1;
    for(int ipar=0;ipar<6;ipar++){
      Par2[ipar]=Par1[ipar];
    }
  }
  
  return 0; 
}

////This functions reads in the antenna tables and interpolates (or extrapolates) from the table to provide output value for raytracing
bool MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint_Table(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary, double IceLayerHeight, int AntennaNumber, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP){

  double AirTxHeight=SrcHeightASL/100;////Height of the source
  double HorizontalDistance=HorizontalDistanceToRx/100;////Horizontal distance
  IceLayerHeight=IceLayerHeight/100;////Height where the ice layer starts off
  //double AntennaDepth=RxDepthBelowIceBoundary/100;////Depth of antenna in the ice

  //cout<<"Table values "<<AirTxHeight<<" "<<HorizontalDistance<<" "<<IceLayerHeight<<" "<<AntennaDepth<<endl;
  // double PropagationTimeIce=0;
  // double PropagationTimeAir=0;
  // double LaunchAngleAir=0;
  // double TotalHorizontalDistanceinAir=0;
  //double TotalHorizontalDistance=0;

  int TotalTableEntries=AllTableAllAntData[AntennaNumber][0].size()-1;
  MaxAirTxHeight=AllTableAllAntData[AntennaNumber][0][0];
  MinAirTxHeight=AllTableAllAntData[AntennaNumber][0][TotalTableEntries];
  
  double x1=0,x2=0,y1=0,y2=0;
  double AirTxHeight1,AirTxHeight2;
  double Par1[6];
  double Par2[6];
  double ParInterpolatedValues[6];
  ///0 is OpticalPathIce
  ///1 is OpticalPathAir
  ///2 is LaunchAngleAir
  ///3 is THDAir

  // cout<<"in here "<<AirTxHeight<<" "<<HorizontalDistance<<endl;
  if(AirTxHeight<=MaxAirTxHeight && AirTxHeight>=MinAirTxHeight && AirTxHeight>0){
    //cout<<"in here "<<endl;
   
    GetParValues(AntennaNumber,AirTxHeight,HorizontalDistance,IceLayerHeight,AirTxHeight1,Par1,AirTxHeight2,Par2);
    x1=AirTxHeight1;
    x2=AirTxHeight2;
    for(int ipar=0;ipar<6;ipar++){
      y1=Par1[ipar];
      y2=Par2[ipar];
      double GetInterpolatedParValue=0;
      int checkval=0;
      if(y1==-pow(10,9) || y2==-pow(10,9)){
	checkval=1;
      }
      if(x1!=x2 && checkval==0){
	GetInterpolatedParValue=oneDLinearInterpolation(AirTxHeight,x1,y1,x2,y2);
      }else{////If AirTxHeight is exactly equal to a bin value then AirTxHeight1 and AirTxHeight2 are equal and we dont need to do interpolation for AirTxHeight2 as we already have done it for AirTxHeight1
	if(x1==x2 && y1==y2){
	  GetInterpolatedParValue=Par1[ipar];
	}
	if(y2==-pow(10,9) && y1==-pow(10,9)){
	  ipar=5;
	}

      }
      ParInterpolatedValues[ipar]=GetInterpolatedParValue;
    }
    
  }
  
  opticalPathLengthInIce=ParInterpolatedValues[0]*100;
  opticalPathLengthInAir=ParInterpolatedValues[1]*100;
  launchAngle=ParInterpolatedValues[2]*(RayTracingFunctions::pi/180);
  horizontalDistanceToIntersectionPoint=ParInterpolatedValues[3]*100;
  reflectionCoefficientS=ParInterpolatedValues[4];
  reflectionCoefficientP=ParInterpolatedValues[5];

  bool CheckSolBool=false;
  if( (y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9)) ){
    CheckSolBool=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint(SrcHeightASL, HorizontalDistanceToRx, RxDepthBelowIceBoundary, IceLayerHeight*100, opticalPathLengthInIce, opticalPathLengthInAir, launchAngle, horizontalDistanceToIntersectionPoint,reflectionCoefficientS,reflectionCoefficientP);
    //cout<<"in here 2 "<<opticalPathLengthInIce<<" "<<opticalPathLengthInAir<<" "<<launchAngle<<" "<<horizontalDistanceToIntersectionPoint<<" "<<y1<<" "<<y2<<endl;
  }
  
  bool CheckSolution=true;

  if(y2==-pow(10,9) && y1==-pow(10,9)){
    // std::cout<<"The given Total Horizontal Distance bigger than the maximum percentage of the Total Horizontal Distance range for the given AirTxHeight! Cannot extrapolate! "<<HorizontalDistance<<" "<<AirTxHeight<<" "<<IceLayerHeight<<std::endl;
    CheckSolution=false;
  }
  if( ((y1==-pow(10,9) && y2!=-pow(10,9)) || (y2==-pow(10,9) && y1!=-pow(10,9))) && CheckSolBool==false ){
    // std::cout<<"The minimiser failed for the weird extrapolation case"<<std::endl;
    CheckSolution=false;    
  }
  if(AirTxHeight>MaxAirTxHeight){
//    std::cout<<"AirTx Height is greater than the maximum possible height in GDAS!"<<std::endl;
    CheckSolution=false;
  }
  if(AirTxHeight<MinAirTxHeight){
//    std::cout<<"AirTx Height is less than the minimum possible height in the table!"<<std::endl;
    CheckSolution=false;
  }
  if(AirTxHeight<0){
//    std::cout<<"AirTx Height is less than the zero!"<<std::endl;
    CheckSolution=false;
  }

  if(CheckSolution==false){
    opticalPathLengthInIce=0;
    opticalPathLengthInAir=0;
    launchAngle=0;
    horizontalDistanceToIntersectionPoint=0;
  }
  
  return CheckSolution;
}


int MultiRayAirIceRefraction::MakeRayTracingTable(double AntennaDepth, double IceLayerHeight, int AntennaNumber){

  std::vector<std::vector <double>> AllTableData;
  
  ////convert cm to m
  AntennaDepth=AntennaDepth/100;
  IceLayerHeight=IceLayerHeight/100;
  
  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();
  
  ////Print out the entry number, the Tx height, ice layer height, Tx height above the icelayer height, total horizontal distance on surface, total horizontal distance in ice, RayLaunchAngle at Tx, incident angle on ice and recievd angle in ice at the antenna inside this file
  RayTracingFunctions::MakeAtmosphere();
 
  ////Define variables for the loop over Tx height and ray launch angle
  double RayLaunchAngleInAir=0;////Set zero for now and 0 deg straight up. This variable defines the initial launch angle of the ray w.r.t to the vertical in the atmosphere. 0 deg is straight up
  //double AirTxHeight=h_data[h_data.size()-1][h_data[h_data.size()-1].size()-1];////Maximum height available with the refractive index data
  double AirTxHeight=100000;////Maximum height available with the refractive index data
  
  // ////Set the variables for the for loop that will loop over the launch angle values. All values are in degrees
  // double AngleStepSize=0.5;
  // double LoopStartAngle=92;
  // double LoopStopAngle=178;
  // int TotalAngleSteps=floor((LoopStopAngle-LoopStartAngle)/AngleStepSize);

  ////Set the variables for the for loop that will loop over the Tx height values above the ice layer. All values are in degrees
  LoopStartHeight=AirTxHeight;
  LoopStopHeight=IceLayerHeight;
  TotalHeightSteps=floor((LoopStartHeight-LoopStopHeight)/HeightStepSize)+1;

  ////Temporary vectors
  std::vector <double> temp1;
  std::vector <double> temp2;
  std::vector <double> temp3;
  std::vector <double> temp4;
  std::vector <double> temp5;
  std::vector <double> temp6;
  std::vector <double> temp7;
  std::vector <double> temp8;
  
  int ifileentry=0;
  //cout<<" our final angle is "<<LoopStartAngle+AngleStepSize*TotalAngleSteps<<endl;
  ////Start looping over the Tx Height and Launch angle values
  for(int ihei=0;ihei<TotalHeightSteps;ihei++){
    AirTxHeight=LoopStartHeight-HeightStepSize*ihei;
    if(AirTxHeight!=IceLayerHeight && ihei==TotalHeightSteps-1){
      AirTxHeight=IceLayerHeight;
    }
    if(AirTxHeight>0){
      for(int iang=0;iang<TotalAngleSteps;iang++){
	RayLaunchAngleInAir=LoopStartAngle+AngleStepSize*iang;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Section for propogating the ray through the atmosphere
	
	////Find out how many atmosphere layers are above the source or Tx which we do not need
	int skiplayer=0;
	for(int ilayer=RayTracingFunctions::MaxLayers;ilayer>-1;ilayer--){
	  //std::cout<<ilayer<<" "<<RayTracingFunctions::ATMLAY[ilayer]/100<<" "<<RayTracingFunctions::ATMLAY[ilayer-1]/100<<std::endl;
	  if(AirTxHeight<RayTracingFunctions::ATMLAY[ilayer]/100 && AirTxHeight>=RayTracingFunctions::ATMLAY[ilayer-1]/100){
	    //std::cout<<"Tx Height is in this layer with a height range of "<<MultiRayAirIceRefraction::RayTracingFunctions::ATMLAY[ilayer]/100<<" m to "<<MultiRayAirIceRefraction::RayTracingFunctions::ATMLAY[ilayer-1]/100<<" m and is at a height of "<<AirTxHeight<<" m"<<std::endl;
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
      
	////Define variables for ray propogation through mutliple layers in the atmosphere
	double Start_nh=0;
	double StartHeight=0;
	double StopHeight=0;
	double StartAngle=0;
	double TotalHorizontalDistanceInAir=0;
	double TimeInAir=0;
      
	////Start loop over the atmosphere layers and analyticaly propagate the ray through the atmosphere
	//std::cout<<"Fitting the atmosphere refrative index profile with multiple layers and propogate the ray"<<std::endl;
	for(int ilayer=RayTracingFunctions::MaxLayers-SkipLayersAbove-1;ilayer>SkipLayersBelow-1;ilayer--){
	  //std::cout<<B_air<<std::endl;
	  ////Set the starting height of the ray for propogation for that layer
	  if(ilayer==RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	    ////If this is the first layer then set the start height to be the height of the source
	    StartHeight=AirTxHeight;
	  }else{
	    ////If this is any layer after the first layer then set the start height to be the starting height of the layer
	    StartHeight=RayTracingFunctions::ATMLAY[ilayer+1]/100-0.00001;
	  }
	
	  ////Since we have the starting height now we can find out the refactive index at that height
	  Start_nh=RayTracingFunctions::Getnz_air(StartHeight);
	
	  ////Set the staopping height of the ray for propogation for that layer
	  if(ilayer==(SkipLayersBelow-1)+1){
	    ////If this is the last layer then set the stopping height to be the height of the ice layer
	    StopHeight=IceLayerHeight;
	  }else{
	    ////If this is NOT the last layer then set the stopping height to be the end height of the layer
	    StopHeight=RayTracingFunctions::ATMLAY[ilayer]/100;
	  }
	
	  ////If this is the first layer then set the initial launch angle of the ray through the layers. I calculate the final launch angle by doing 180-RayLaunchAngleInAir since my raytracer only works with 0 to 90 deg. Setting an angle of 95 deg w.r.t to the vertical where 0 is up means that my raytraces takes in an launch angle of 85.
	  if(ilayer==RayTracingFunctions::MaxLayers-SkipLayersAbove-1){
	    StartAngle=180-RayLaunchAngleInAir;
	  }
	  //std::cout<<ilayer<<" Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingHeight="<<StartHeight<<" ,StoppingHeight="<<StopHeight<<" ,RayLaunchAngle"<<StartAngle<<std::endl;
	
	  ////Get the hit parameters from the function. The output is:
	  //// How much horizontal distance did the ray travel in the layer
	  //// The angle of reciept/incidence at the end or the starting angle for propogation through the next layer
	  //// The value of the L parameter for that layer
	  double* GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, StopHeight, StartHeight, StartAngle, 1);
	  TotalHorizontalDistanceInAir+=GetHitPar[0];
	  StartAngle=GetHitPar[1];
	  TimeInAir+=GetHitPar[3];
	
	  ////dont forget to delete the pointer!
	  delete []GetHitPar;
	}
      
	double IncidentAngleonIce=StartAngle;
	//std::cout<<"Total horizontal distance travelled by the ray using Multiple Layer fitting is "<<TotalHorizontalDistance<<std::endl;
      
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Section for propogating the ray through the ice
            
	////Set the starting depth of the ray for propogation to at the ice surface
	double StartDepth=0.0;
	////Since we have the starting height of the ice layer we can find out the refactive index of air at that height
	Start_nh=RayTracingFunctions::Getnz_air(IceLayerHeight);
	
	////Set the stopping depth of the ray for propogation to be the depth of the antenna
	StopHeight=AntennaDepth;
	////Set the initial launch angle or the angle of incidence of the ray
	StartAngle=IncidentAngleonIce;
	//std::cout<<"Starting n(h)="<<Start_nh<<" ,A="<<A<<" ,B="<<B<<" ,C="<<C<<" StartingDepth="<<StartDepth<<" ,StoppingDepth="<<AntennaDepth<<" ,RayLaunchAngle="<<StartAngle<<std::endl;
    
	////Get the hit parameters from the function. The output is:
	//// How much horizontal distance did the ray travel through ice to hit the antenna
	//// The angle of reciept/incidence at the end at the antenna
	//// The value of the L parameter for whole atmosphere fit
	double *GetHitPar=RayTracingFunctions::GetLayerHitPointPar(Start_nh, AntennaDepth,StartDepth, StartAngle, 0);
      
	////SLF here stands for Single Layer Fitting. These variables store the hit parameters
	double TotalHorizontalDistanceInIce=GetHitPar[0];
	double RecievdAngleInIce=GetHitPar[1];
	//double LvalueIce=GetHitPar[2];
	double TimeInIce=GetHitPar[3];
      
	//std::cout<<"Total horizontal distance travelled by the ray in ice is  "<<TotalHorizontalDistanceInIce<<std::endl;
      
	//if(isnan(TotalHorizontalDistanceInAir)==false){
	
	////define dummy/temporary variables for storing data
	double dummy[16];
        for(int idum=0;idum<16;idum++){
          dummy[idum]=0;
        }

        dummy[0]=ifileentry;
        dummy[1]=AirTxHeight;
        dummy[2]=TotalHorizontalDistanceInAir + TotalHorizontalDistanceInIce;
        dummy[3]=TotalHorizontalDistanceInAir;
        dummy[4]=TotalHorizontalDistanceInIce;
        dummy[5]=(TimeInIce+TimeInAir)*RayTracingFunctions::spedc;
        dummy[6]=TimeInAir*RayTracingFunctions::spedc;
        dummy[7]=TimeInIce*RayTracingFunctions::spedc;
        dummy[8]=(TimeInIce+TimeInAir)*pow(10,9);
        dummy[9]=TimeInAir*pow(10,9);
        dummy[10]=TimeInIce*pow(10,9);
        dummy[11]=RayLaunchAngleInAir;
        dummy[12]=IncidentAngleonIce;
        dummy[13]=RecievdAngleInIce;
        dummy[14]=RayTracingFunctions::Refl_S(IncidentAngleonIce*(RayTracingFunctions::pi/180.0), IceLayerHeight);
        dummy[15]=RayTracingFunctions::Refl_P(IncidentAngleonIce*(RayTracingFunctions::pi/180.0), IceLayerHeight);

        temp1.push_back(dummy[1]);///AirTx Height
        temp2.push_back(dummy[2]);///THDTotal          
        temp3.push_back(dummy[7]);///OpticalPathIce
        temp4.push_back(dummy[6]);///OpticalPathAir
        temp5.push_back(dummy[11]);///LaunchAngleAir          
        temp6.push_back(dummy[3]);///THDAir
        temp7.push_back(dummy[14]);///ReflectionCoefficientS
        temp8.push_back(dummy[15]);///ReflectionCoefficientP
	ifileentry++;
	//}
	delete[] GetHitPar;
     
      }//// end of iang loop
    }////if condition to make sure AirTxHeight does not go below zero
  }//// end of ihei loop

  
  AllTableData.push_back(temp1);
  AllTableData.push_back(temp2);
  AllTableData.push_back(temp3);
  AllTableData.push_back(temp4);
  AllTableData.push_back(temp5);
  AllTableData.push_back(temp6);
  AllTableData.push_back(temp7);
  AllTableData.push_back(temp8);

  AllTableAllAntData.push_back(AllTableData);
  
  AllTableData.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();
  temp4.clear();
  temp5.clear();
  temp6.clear();
  temp7.clear();
  temp8.clear();
  
  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  durationb=durationb/1000000;
  std::cout<<"total time taken by the script to generate the table: "<<durationb<<" s"<<std::endl;
  return 0;
  
}
