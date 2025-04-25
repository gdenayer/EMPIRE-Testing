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
libempire_api.EMPIRE_API_Connect(b"pong.xml");

# Receive signal (signal1)
signal1 = (c_double * 1)(-1);
size1 = c_int(1);
print('Receving ...')
libempire_api.EMPIRE_API_recvSignal_double(b"signal1", size1, signal1);
print('Received: {}'.format(signal1[0]));

# Send signal (signal2)
signal2 = (c_double * 1)(10000);
size2 = c_int(1);
print('Sending ...')
libempire_api.EMPIRE_API_sendSignal_double(b"signal2", size2, signal2);
print('Sent: {}'.format(signal2[0]));

# Disconnect
libempire_api.EMPIRE_API_Disconnect();
