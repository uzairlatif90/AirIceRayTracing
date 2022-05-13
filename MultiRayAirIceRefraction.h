
#ifndef _INCLUDE_MULTIRAYAIRICEREFRACTION_H_
#define _INCLUDE_MULTIRAYAIRICEREFRACTION_H_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstdlib>

#include <sys/time.h>

#include "RayTracingFunctions.h"

namespace MultiRayAirIceRefraction{

  /********Stuff for Interpolation**********/
  static std::vector <double> GridPositionH;
  static std::vector <double> GridPositionTh;
  static std::vector <double> GridZValue[8];

  static double GridStartTh=90.1;
  static double GridStopTh=179.9;
  
  static double GridStepSizeH_O=100;
  static double GridStepSizeTh_O=0.1;
  static double GridWidthH=1000;
  static double GridWidthTh=GridStopTh-GridStartTh;

  static int GridPoints=100;////just set a non-zeronumber for now
  static int TotalStepsH_O=100;////just set a non-zeronumber for now
  static int TotalStepsTh_O=100;////just set a non-zeronumber for now
  static double GridStartH=1000;////just set a non-zeronumber for now
  static double GridStopH=100000;////just set a non-zeronumber for now  

  ////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.  
  bool GetHorizontalDistanceToIntersectionPoint(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary,double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint,  double& reflectionCoefficientS, double& reflectionCoefficientP);
  
  ////Just a simple function for interpolating in 1 dimension between two points (xa,ya) and (xb,yb)
  double oneDLinearInterpolation(double x, double xa, double ya, double xb, double yb);

  ////Just a simple function for extrapolating in 1 dimension outside two points
  double Extrapolate(int Par, int index, double TotalHorizontalDistance, int AntennaNumber);

  ////Find the limit to extrapolate for all the parameters using the Air Launch Angle parameter
  double FindExtrapolationLimit(int index, double TotalHorizontalDistance, int AntennaNumber);

  int FindClosestAirTxHeight(double IceLayerHeight, double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber);
  
  int FindClosestTHD(double ParValue, int StartIndex, int EndIndex, int &RStartIndex, int &REndIndex, double &ClosestVal, int AntennaNumber);
  
  ////Interpolate the value of the given parameter for a given TxHeight and THD
  int GetParValues(double AntennaNumber, double AirTxHeight, double TotalHorizontalDistance, double IceLayerHeight, double &AirTxHeight1, double Par1[6],double &AirTxHeight2, double Par2[6]);

  ////This functions reads in the antenna tables and interpolates (or extrapolates) from the table to provide output value for raytracing
  bool GetHorizontalDistanceToIntersectionPoint_Table(double SrcHeightASL, double HorizontalDistanceToRx ,double RxDepthBelowIceBoundary, double IceLayerHeight, int AntennaNumber, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint, double& reflectionCoefficientS, double& reflectionCoefficientP);

  void Air2IceRayTracing(double AirTxHeight, double HorizontalDistance, double IceLayerHeight,double AntennaDepth, double StraightAngle, double dummy[14]);
  
  void MakeTable(double IceLayerHeight,double AntennaDepth);

  double GetInterpolatedValue(double hR, double thR, int rtParameter);
  
  
  void GetRayTracingSolutions(double RayLaunchAngleInAir,double AirTxHeight, double IceLayerHeight, double AntennaDepth, double dummy[16]);
  
  ////This is the main function which will make the RayTracing Table that will be used for interpolation by COREAS
  ////The arguments are:
  ////1. AntennaDepth is depth of antenna in the ice and is given in m and is positive
  ////2. IceLayerHeight is the height in m in a.s.l where the ice layer starts off
  int MakeRayTracingTable(double AntennaDepth, double IceLayerHeight, int AntennaNumber);
  
}
#endif
