"""
Implementation of binary file I/O for C++ STL and GMM classes.
"""

# $Source: /usr/data0/leipzig_work/tmat_cvs/src/binaryIO.py,v $ #

# requires struct, string, Exception
from numpy import *
from struct import *
from string import *


__all__ = ['set_IO_sizeBits', 
  'loadSTDVectorLong', 'loadSTDVectorDouble', 'loadSTDVectorComplexDouble',
  'loadSTDVectorSTDVectorLong', 'loadSTDVectorSTDVectorDouble', 'loadSTDVectorSTDVectorComplexDouble',   
  'loadSTDString', 
  'loadGMMDenseMatrixLong', 'loadGMMDenseMatrixDouble', 'loadGMMDenseMatrixComplexDouble', 
  'loadntuple2Long', 'loadntuple2_ulong', 'loadntuple2_size_t',
  'loadntuple2Double', 'loadntuple2ComplexDouble', 
  'loadntuple3Long', 'loadntuple3_ulong', 'loadntuple3_size_t', 
  'loadntuple3Double', 'loadntuple3ComplexDouble', 
  'load_size_t', 'load_ulong', 'load_long', 'load_double', 'load_complex_double', 'load_int32_t',
  'loadMappedBinary', 'load_simple_object', 'interp_data']
  
class Error(Exception):
    """Base class for exceptions in this module.
		   
			 Attributes:
			   message: explaination of the error
		"""

    def __init__(self,message):
	    self.message = message
			
class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

global s_size_t_type, s_ulong_type, s_long_type
s_size_t_type = '=Q'
s_ulong_type = '=Q'
s_long_type = '=q'

def set_IO_sizeBits(nBits):
  """
  set number of bits in imported representation of size_t:
  """
  global s_size_t_type, s_ulong_type, s_long_type
  if nBits == 32:
      s_size_t_type = '=L' # 32-bit unsigned long  
      s_ulong_type = '=L'  # 32-bit unsigned long
      s_long_type = '=l'   # 32-bit long
  elif (nBits == 64):
      s_size_t_type = '=Q' # 64-bit unsigned long
      s_ulong_type  = '=Q' # 64-bit unsigned long
      s_long_type   = '=q' # 64-bit long
  else:
      raise Error, 'no local representation for %d bit size_t' % nBits	 

set_IO_sizeBits(64) # default case


def loadSTDVectorLong(fp):
  """ V = loadSTDVectorLong(fp)
  %
  % load a std::vector<long> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % V = loadSTDVectorLong(fp);
  % fclose(fp);
  """

  V = loadSTDVector(fp,s_long_type)
  return V


def loadSTDVectorDouble(fp):
  """ V = loadSTDVectorDouble(fp)
  %
  % load a std::vector<double> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % V = loadSTDVectorDouble(fp);
  % fclose(fp);
  """
  V = loadSTDVector(fp,'d')
  return V
	 
def loadSTDVectorComplexDouble(fp):
  """ V = loadSTDVectorComplexDouble(fp)
  %
  % load a std::vector< std::complex<double> > from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % V = loadSTDVectorComplexDouble(fp);
  % fclose(fp);
  """
  V = loadSTDVector(fp,'dd')
  return V

def loadSTDVectorSTDVectorLong(fp):
  VV = loadSTDVectorSTDVector(fp,s_long_type)
  return VV

def loadSTDVectorSTDVectorDouble(fp):
  VV = loadSTDVectorSTDVector(fp,'d')
  return VV

def loadSTDVectorSTDVectorComplexDouble(fp):
  VV = loadSTDVectorSTDVector(fp,'dd')
  return VV
  
def loadSTDString(fp):
  """ s = loadSTDloadSTDString(fp)
  %
  % load a std::string from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % s = loadSTDString(fp);
  % fclose(fp);
  """
  s = ''
  # read <unsigned int> size:
  nSize = load_size_t(fp)
  if (0 < nSize):
    fmt = '%dc' % nSize
    s = fp.read(calcsize(fmt))			
  return s

		
def loadGMMDenseMatrixLong(fp):
  """ A = loadGMMDenseMatrixLong(fp)
  %
  % load a gmm::dense_matrix<long> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadGMMDenseMatrixLong(fp);
  % fclose(fp);
  """
  A = loadGMMDenseMatrix(fp,s_long_type)
  return A

	 
def loadGMMDenseMatrixDouble(fp):
  """ A = loadGMMDenseMatrixDouble(fp)
  %
  % load a gmm::dense_matrix<double> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadGMMDenseMatrixDouble(fp);
  % fclose(fp);
  """
  A = loadGMMDenseMatrix(fp,'d')
  return A


def loadGMMDenseMatrixComplexDouble(fp):
  """ A = loadGMMDenseMatrixComplexDouble(fp)
  %
  % load a gmm::dense_matrix< std::complex<double> > from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadGMMDenseMatrixComplexDouble(fp);
  % fclose(fp);
  """
  A = loadGMMDenseMatrix(fp,'dd')
  return A


def loadntuple2Long(fp):
  """ A = loadntuple2Long(fp)
  %
  % load a ntuple<long,2> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple2Long(fp);
  % fclose(fp);
  """
  A = loadntuple2(fp,s_long_type)
  return A	 


def loadntuple2_ulong(fp):
  """ A = loadntuple2_ulong(fp)
  %
  % load a ntuple<unsigned long,2> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple2_ulong(fp);
  % fclose(fp);
  """
  A = loadntuple2(fp,s_ulong_type)
  return A	 

def loadntuple2_size_t(fp):
  """ A = loadntuple2_size_t(fp)
  %
  % load a ntuple<size_t,2> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple2_size_t(fp);
  % fclose(fp);
  """
  A = loadntuple2(fp,s_size_t_type)
  return A	 
  
    
def loadntuple2Double(fp):
  """ A = loadntuple2Double(fp)
  %
  % load a ntuple<double,2> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple2Double(fp);
  % fclose(fp);
  """
  A = loadntuple2(fp,'d')
  return A	 
 

def loadntuple2ComplexDouble(fp):
  """ A = loadntuple2ComplexDouble(fp)
  %
  % load a ntuple< std::complex<double>,2> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple2ComplexDouble(fp);
  % fclose(fp);
  """
  A = loadntuple2(fp,'dd')
  return A	 

def loadntuple3Long(fp):
  """ A = loadntuple3Long(fp)
  %
  % load a ntuple<long,3> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple3Long(fp);
  % fclose(fp);
  """
  A = loadntuple3(fp,s_long_type)
  return A	 


def loadntuple3_ulong(fp):
  """ A = loadntuple3_ulong(fp)
  %
  % load a ntuple<unsigned long,3> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple3_ulong(fp);
  % fclose(fp);
  """
  A = loadntuple3(fp,s_ulong_type)
  return A	

def loadntuple3_size_t(fp):
  """ A = loadntuple3_size_t(fp)
  %
  % load a ntuple<size_t,3> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple3_size_t(fp);
  % fclose(fp);
  """
  A = loadntuple3(fp,s_size_t_type)
  return A	 
  
   
def loadntuple3Double(fp):
  """ A = loadntuple3Double(fp)
  %
  % load a ntuple<double,3> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple3Double(fp);
  % fclose(fp);
  """
  A = loadntuple3(fp,'d')
  return A	 
 

def loadntuple3ComplexDouble(fp):
  """ A = loadntuple3ComplexDouble(fp)
  %
  % load a ntuple< std::complex<double>,3> from file specified by file-pointer fp.
  %
  % example of usage:
  %
  % fp = fopen(sFileName,'rb');
  % if (-1 == fp)
  %   error(sprintf('unable to open %s for binary read',sFileName));
  % end
  % A  = loadntuple3ComplexDouble(fp);
  % fclose(fp);
  """
  A = loadntuple3(fp,'dd')
  return A	 

  
def load_size_t(fp):
  """ load a size_t from file fp 
  """
  s = fp.read(calcsize(s_size_t_type))
  n = unpack(s_size_t_type,s)[0]
  return n

    
def load_ulong(fp):
  """ load an unsigned long from file fp 
  """
  s = fp.read(calcsize(s_ulong_type))
  ul = unpack(s_ulong_type,s)[0]
  return ul

 
def load_long(fp):
  """ load a long from file fp 
  """
  s = fp.read(calcsize(s_long_type))
  l = unpack(s_long_type,s)[0]
  return l


def load_double(fp):
  """ load a double from file fp 
  """
  s = fp.read(calcsize('d'))
  d = unpack('d',s)[0]
  return d


def load_complex_double(fp):
  """ load a std::complex<double> from file fp 
  """
  s = fp.read(calcsize('dd'))
  c_ = unpack('dd',s)
  return complex(c_[0],c_[1])

def load_int32_t(fp):
  """ load an int32_t from file fp 
  """
  s = fp.read(calcsize('=l'))
  n = unpack('=l',s)[0]
  return n

def loadSTDVector(fp,fmt):
  """ 
     V = loadSTDVector(fp,fmt)
     load STD vector of \"struct\" fmt from file-pointer fp
  """
  if (fmt == 'd' or fmt =='L' or fmt == 'l' or fmt == 'Q' or fmt == 'q'):
    dataType = float64
  elif (fmt == 'dd'):
    dataType = complex128
  else:
    raise Error, 'unknown data type %s' % fmt

  nSize = load_size_t(fp)
  V = zeros((nSize,),dtype=dataType)
  if (0 < nSize):
    if (dataType != complex128):
      fmt_ = '%d%s' % (nSize, fmt)
      s = fp.read(calcsize(fmt_))
      V[:] = unpack(fmt_,s)
    else:
      fmt_ = '%d%s' % (nSize*2,'d')
      s = fp.read(calcsize(fmt_))
      V_ = unpack(fmt_,s)  # there has _got_ to be a better way...
      V[:].real = V_[0:size(V_):2]
      V[:].imag = V_[1:size(V_):2]

  return V
    
def loadSTDVectorSTDVector(fp,fmt):
  """ 
     VV = loadSTDVectorSTDVector(fp,fmt)
     load STD vector of STD vector of \"struct\" fmt from file-pointer fp
  """

  nSize = load_size_t(fp)
  VV = []
  if (0 < nSize):
    for n in range(nSize):
      VV.append(loadSTDVector(fp,fmt))
    
  return VV
      
def loadGMMDenseMatrix(fp,fmt):
  """ A = loadGMMDenseMatrix(fp,fmt)
   load GMM::dense_matrix of \"struct\" fmt from file-pointer fp
  """

  if (fmt == 'd' or fmt =='L' or fmt == 'l' or fmt == 'Q' or fmt == 'q'):
    dataType = float64
  elif (fmt == 'dd'):
    dataType = complex128
  else:
    raise Error, 'unknown data type %s' % fmt

  nRows = load_size_t(fp);
  nCols = load_size_t(fp);

  A = zeros((nRows, nCols),dtype=dataType)

  if (0 < nRows*nCols):
    fmt_ = '%d%s' % (nCols, fmt)
    # read row-major order:
    for row in range(nRows):
      s = fp.read(calcsize(fmt_))
      A[row,:] = unpack(fmt_,s)

  return A


def loadntuple2(fp,fmt):
  """ n2 = loadntuple2(fp,fmt)
     load ntuple<...,2> of \"struct\" fmt from file-pointer fp
  """
  if (fmt == 'd' or fmt =='L' or fmt == 'l' or fmt == 'Q' or fmt == 'q'):
    dataType = float64
  elif (fmt == 'dd'):
    dataType = complex128
  else:
    raise Error, 'unknown data type %s' % fmt

  nSize = 2
  n2 = zeros((nSize,),dtype=dataType)
  if (dataType != complex128):
    fmt_ = '%d%s' % (nSize, fmt)
    s = fp.read(calcsize(fmt_))
    n2[:] = unpack(fmt_,s)
  else:
    fmt_ = '%d%s' % (nSize*2,'d')
    s = fp.read(calcsize(fmt_))
    n3_ = unpack(fmt_,s)  # there has _got_ to be a better way...
    n2[:].real = n2_[0:size(V_):2]
    n2[:].imag = n2_[1:size(V_):2]

  return n2
  
def loadntuple3(fp,fmt):
  """ n3 = loadntuple(fp,fmt)
     load ntuple<...,3> of \"struct\" fmt from file-pointer fp
  """
  if (fmt == 'd' or fmt =='L' or fmt == 'l' or fmt == 'Q' or fmt == 'q'):
    dataType = float64
  elif (fmt == 'dd'):
    dataType = complex128
  else:
    raise Error, 'unknown data type %s' % fmt

  nSize = 3
  n3 = zeros((nSize,),dtype=dataType)
  if (dataType != complex128):
    fmt_ = '%d%s' % (nSize, fmt)
    s = fp.read(calcsize(fmt_))
    n3[:] = unpack(fmt_,s)
  else:
    fmt_ = '%d%s' % (nSize*2,'d')
    s = fp.read(calcsize(fmt_))
    n3_ = unpack(fmt_,s)  # there has _got_ to be a better way...
    n3[:].real = n3_[0:size(V_):2]
    n3[:].imag = n3_[1:size(V_):2]

  return n3
 
def loadStringList(fp):
  """ SL = loadStringList(fp)
     load list of strings from file-pointer fp
  """
  nSize = load_size_t(fp)
  SL = [ loadSTDString(fp) for n in range(nSize) ]
  return SL
        
def readHeaderDict(fp):
  """
  A work-around for the problem that the lambda-expression in "loadMappedBinary" confuses
  the execution of the "exec" clause
  """
  MAX_HEADERLEN = 256
  DICT = None

  header_str = fp.read(MAX_HEADERLEN)
  header_beg = header_str.find('$HEADER$')
  if (-1 != header_beg):
    # format of header: "$HEADER$ <any excluding '{'|'}'> <python dict syntax>  <any>"
    pos = header_str[header_beg:].find('{')
    if (-1 != pos):
      header_beg += pos
    else:
      header_beg = -1
    header_end = header_beg
    pos = header_str[header_beg:].find('}')     
    if (-1 != pos):
      header_end += pos + 1
    if (-1 != header_beg) and (-1 != header_end):
      exec 'DICT = ' + header_str[header_beg: header_end]
  
  return DICT
  

def loadMappedBinary(fp, header=False):
  """ DICT = loadMappedBinary(fp)
     load all mapped binary objects from file-pointer fp
  """
  DICT = None

  if header:
    DICT = readHeaderDict(fp)
  else:
    DICT = {}
    map = loadStringList(fp)
    for format in map:
      """
      # *** DEBUG ***
      print format
      """

      if (-1 != format.find('$HEADER$')):
        continue # skip the header field

      format_ = split(format) # break format string into tokens [<variable name>, <variable type: 'r','c','s'>, <container type: 's','v','a','n3'>]
      var_name = format_[0]       # name of the mapped variable
      var_type = format_[1]       # type: one of "r":real, "c":complex, "s":string
      container_type = format_[2] # container type: one of "s":scalar, "v":vector, "a":array, "n3":ntuple<...,3>
      var_value = {
        's': lambda: {
               'r': lambda: load_double(fp),
               'c': lambda: load_complex_double(fp),
               's': lambda: loadSTDString(fp)
             }[var_type](), 
        'v': lambda: {
               'r': lambda: loadSTDVectorDouble(fp),
               'c': lambda: loadSTDVectorComplexDouble(fp),
               's': lambda: error('load std::vector<std::string> not implemented')
             }[var_type](),       
        'a': lambda: {
               'r': lambda: loadGMMDenseMatrixDouble(fp),
               'c': lambda: loadGMMDenseMatrixComplexDouble(fp),
               's': lambda: error('load gmm::dense_matrix<std::string> not implemented')
             }[var_type](), 
        'n3': lambda: {
               'r': lambda: loadntuple3Double(fp),
               'c': lambda: loadntuple3ComplexDouble(fp),
               's': lambda: error('load ntuple<std::string,3> not implemented')
             }[var_type](), 
        'vv': lambda: {
               'r': lambda: loadSTDVectorSTDVectorDouble(fp),
               'c': lambda: loadSTDVectorSTDVectorComplexDouble(fp),
               's': lambda: error('load std::vector<std::vector<std::string> > not implemented')
             }[var_type]()            
      }[container_type]()
      # add the variable to the dictionary:
      DICT[var_name] = var_value

  return DICT  



def load_simple_object(fp):
  """
  load_simple_object(fp):
    method corresponding to "python_util::simple_object_base<std::complex<double>, double, long>::readBinaryVirtual" 
  """
  etype = load_ulong(fp)
  result = {
    0: lambda: error('load_simple_object: unknown data type enum: %d' % etype),
    1: lambda: load_complex_double(fp),
    2: lambda: load_double(fp),
    3: lambda: load_long(fp),
    4: lambda: bool(load_ulong(fp)),
    5: lambda: loadSTDString(fp),
    6: lambda: load_object_list(fp),
    7: lambda: load_object_map(fp)
  }[etype]()

  return result
  
def load_object_list(fp):
  """
  load_object_list(fp):
    load list of simple_object
  """  
  N = load_ulong(fp)
  result = [ load_simple_object(fp) for n in range(N) ] 

  return result
  
def load_object_map(fp):
  """
  load_object_map(fp):
    load string keyed map of simple_object
  """  
  N = load_ulong(fp)
  result = dict((loadSTDString(fp), load_simple_object(fp)) for n in range(N)) 

  return result

  
def error(msg):
  raise Error, msg
  return None
    
def interp_data(x, y, N, s_=0.0):
  """ [x_i, y_i] = interp_data(x, y, N)
      smoothly interpolate ordinate data "y" along abscissa "x" with "N" values
  """
  from scipy import linspace,interpolate
  x_i = linspace(min(x), max(x), N)
  sy_i = interpolate.splrep(x, y, s=s_) # spline control points  
  y_i = interpolate.splev(x_i, sy_i, der=0)
  return (x_i,y_i)
  
  
