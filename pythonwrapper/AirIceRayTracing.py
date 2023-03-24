import ctypes
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
handle = ctypes.CDLL(dir_path + "/libAirIceRayTracing.so")     

handle.Py_TraceIceToAir.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double*10] 
  
def Py_TraceIceToAir( AntennaDepth, IceLayerHeight, AirTxHeight, HorizontalDistance, ArrayParameters):
    return handle.Py_TraceIceToAir(AntennaDepth, IceLayerHeight, AirTxHeight, HorizontalDistance, ArrayParameters)
