from AirIceRayTracing import *

AntennaDepth=-10
IceLayerHeight=3000
AirTxHeight=8050
HorizontalDistance=10000
ii=ctypes.c_double*10
ArrayParameters=ii(1,1,1,1,1,1,1,1,1,1)
Py_TraceIceToAir(  AntennaDepth, IceLayerHeight, AirTxHeight, HorizontalDistance, ArrayParameters);
for x in ArrayParameters:
    print(x)
