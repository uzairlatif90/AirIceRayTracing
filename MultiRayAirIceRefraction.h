
#ifndef _INCLUDE_MULTIRAYAIRICEREFRACTION_H_
#define _INCLUDE_MULTIRAYAIRICEREFRACTION_H_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstdlib>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>

#include <sys/time.h>

namespace MultiRayAirIceRefraction{


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
  
  ////This is the main function which will make the RayTracing Table that will be used for interpolation by COREAS
  ////The arguments are:
  ////1. AntennaDepth is depth of antenna in the ice and is given in m and is positive
  ////2. IceLayerHeight is the height in m in a.s.l where the ice layer starts off
  int MakeRayTracingTable(double AntennaDepth, double IceLayerHeight, int AntennaNumber);

}
#endif
