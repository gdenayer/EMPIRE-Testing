import os
from ctypes import *
import time

# Load the shared library
print("Using EMPIRE_API:", os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])

# Define C function signatures
libempire_api.EMPIRE_API_Connect.argtypes = [c_char_p]
libempire_api.EMPIRE_API_Connect.restype = c_int

libempire_api.EMPIRE_API_Disconnect.argtypes = []
libempire_api.EMPIRE_API_Disconnect.restype = c_int

libempire_api.EMPIRE_API_sendSignal_double.argtypes = [c_char_p, c_int, POINTER(c_double)]
libempire_api.EMPIRE_API_sendSignal_double.restype = c_int

libempire_api.EMPIRE_API_recvSignal_double.argtypes = [c_char_p, c_int, POINTER(c_double)]
libempire_api.EMPIRE_API_recvSignal_double.restype = c_int

EMPIRE_API_getUserDefinedText = libempire_api.EMPIRE_API_getUserDefinedText
EMPIRE_API_getUserDefinedText.argtypes = [c_char_p]
EMPIRE_API_getUserDefinedText.restype = c_char_p

# Begin API interaction
print("This is the ping!")

libempire_api.EMPIRE_API_Connect(b"./ping.xml")

# Get user-defined text
print("User Message 1:", EMPIRE_API_getUserDefinedText(b"myMessage").decode())

time.sleep(10)  # Edit the XML during this time if needed

# Receive signal
signal2 = (c_double * 10)(*range(1, 11))
size2 = c_int(10)
libempire_api.EMPIRE_API_recvSignal_double(b"signal2", size2, signal2)
print("Ping: signal2 received")

# Send signal
signal1 = (c_double * 5)(1, 2, 3, 4, 5)
size1 = c_int(5)
libempire_api.EMPIRE_API_sendSignal_double(b"signal1", size1, signal1)
print("Ping: signal1 sent")

# Get updated user-defined text
print("User Message 2:", EMPIRE_API_getUserDefinedText(b"myMessage").decode())

# Disconnect
libempire_api.EMPIRE_API_Disconnect()
