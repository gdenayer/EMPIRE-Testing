import os
from ctypes import *
import time

# Load the library
print("Using EMPIRE_API:", os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
print("This is the pong!")

# Define the argument and return types
libempire_api.EMPIRE_API_Connect.argtypes = [c_char_p]
libempire_api.EMPIRE_API_Connect.restype = None

libempire_api.EMPIRE_API_Disconnect.argtypes = []
libempire_api.EMPIRE_API_Disconnect.restype = None

libempire_api.EMPIRE_API_recvSignal_double.argtypes = [c_char_p, c_int, POINTER(c_double)]
libempire_api.EMPIRE_API_recvSignal_double.restype = None

libempire_api.EMPIRE_API_sendSignal_double.argtypes = [c_char_p, c_int, POINTER(c_double)]
libempire_api.EMPIRE_API_sendSignal_double.restype = None

# Connect to the API
libempire_api.EMPIRE_API_Connect(b"pong.xml")

# Receive signal (signal1)
signal1 = (c_double * 5)(1, 2, 3, 4, 5)
size1 = c_int(5)
libempire_api.EMPIRE_API_recvSignal_double(b"signal1", size1, signal1)
print("Pong: signal1 received")

# Send signal (signal2)
signal2 = (c_double * 10)(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
size2 = c_int(10)
libempire_api.EMPIRE_API_sendSignal_double(b"signal2", size2, signal2)
print("Ping: signal2 sent")

# Disconnect
libempire_api.EMPIRE_API_Disconnect()
   
