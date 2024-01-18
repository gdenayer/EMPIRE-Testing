# import necessary modules
import os, sys
from ctypes import cdll
from ctypes import *

# add EMPIRE libraries path
EMPIRELibPath = os.path.abspath(os.environ['EMPIRE_LD_LIBRARY_PATH'])
sys.path.append(EMPIRELibPath)

# import EMPIRE_API library
libempire_mapper = cdll.LoadLibrary(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'])
