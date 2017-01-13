"""
  Python classes and functions for dealing with data generated by Pol-MC.
  
  $Source: /usr/data0/leipzig_work/tmat_cvs/src/polarization_MC_port/polmc_data.py,v $
"""

import sys, os, pdb

from numpy import *
import cStringIO as STR

try:
  import cElementTree as ET
except ImportError:
  try:
    import elementtree.ElementTree as ET
  except ImportError:
    raise ImportError, 'Module \"polmc_data\" requires module \"cElementTree\"(preferable) or module \"elementtree\"'

"""
  # quasi-BNF for data to be extracted from XML:
  # global elements:
  <XML data> ::= <"energyBinNumber" <N-dimensional indexed array> element> # N = 3
                 <"fluenceBinNumber" <N-dimensional indexed array> element>
                 <"source/wavelength" <real number> element>
                 <object data list> 
     
  <object data list>     ::= <object data> [ <object data list> ]
    <object data>        ::= "<object" <"id" integer tag> <"name" string tag> ">"  <global object data> <interface data list>  "</object>
    <global object data>  ::= <"origin": 3-tuple element> <"absorbance": real number element> 
                                   <"scattererType" string element>
                                   <"index": real number element> <"mu_s": real number element> <"mu_a": real number element>
                                   <"indexScatterer": real number element> <"scattererRadius": real number element> 
                                   <"anisotropy": real number element>
                                   <"sampledAnisotropy": real number element>
                                   [ <scattering stats list> ]
    <scattering stats list>  ::= <scattering stats> [ <scattering stats list> ]
      <scattering stats>     ::=  <scattering stats 1D> | <scattering stats 2D>  # you see here why people should learn to "generalize"!  :(
        <scattering stats 1D> ::= "<scatteringStats" <"type"  "1D" tag> ">" <M-dimensional N-tuple array>  "</scatteringstats>"   # M = 1, N = 2
        <scattering stats 2D> ::= "<scatteringStats" <"type"  "2D" tag> ">" <N-dimensional row major array>  "</scatteringstats>"   # N = 2
                                       
  <interface data list>  ::= <interface data> [ <interface data list> ]
     <interface data>     ::= "<interface" <"id" integer tag> <"name" string tag> ">"   
                                  <discrete interface data>  <list interface data> "</interface>"          
     <discrete interface data> ::= <"normal": 3-tuple> <"origin": 3-tuple> 

     <list interface data> ::=  <base vector list> <detection basis list> <Stoke's data list> <scattering stats> 
       <base vector list>    ::=  <base vector data> [ <base vector list> ]
         <base vector data>    ::="<basevector" <"id" string tag> <3-tuple> ">" "</basevector>"
       <detection basis list>  ::=  <detection basis data> [ <base vector list> ] 
         <detection basis data> ::=  "<detectionbasis" <"id" string tag> <3-tuple> ">" "</detectionbasis>"
       <Stoke's data list>      ::=  <Stoke's data> [ <Stoke's data list> ]
         <Stoke's data>         ::=  "<StokesV" <Stoke's data header> <Stoke's data body> "</StokesV>"
           <Stoke's data header> ::= <Stoke's type tag> <"acceptanceCosineIndex" integer tag> ">"
                                         <"minAcceptanceCosine" real number element>
                                         <"maxAcceptanceCosine" real number element>
             <Stoke's type tag>  ::= <"transmittance" real number tag>
           <Stoke's data body>   ::= "<I>" <N-dimensional array> "</I> # N = 2  
                                     "<Q>" <N-dimensional array> "</Q>               
                                     "<U>" <N-dimensional array> "</U>                
                                     "<V>" <N-dimensional array> "</V>                
                                     "<betalin>" <N-dimensional array> "</betalin>                
                                     "<betacirc>" <N-dimensional array> "</betacirc>
      

  <N-dimensional indexed array>  ::= "<" ID ">" <ND indexed array element list> "< /" ID ">"
    <ND indexed array element list> ::= <ND indexed array element> [ "\n" <ND indexed array element list> ]
    <ND indexed array element>      ::= <index list> <real number>
    <index list>         ::= <integer> [ <index list>]     # note: here an N-dimensional array has N entries in the index list
  
  <N-dimensional row major array>  ::= "<" ID ">" <ND array element list> "< /" ID ">"
    <ND array element list> ::= <ND array element> [ "\n" <ND array element list> ]
    <ND array element>      ::= <index list> <real number>
    <index list>         ::= <integer> [ <index list>]     # note: here an N-dimensional array has N entries in the index list
    
  <M-dimensional N-tuple array>  ::= "<" ID ">" <MD-N array element list> "< /" ID ">"
    <MD-N  array element list> ::= <MD-N array element> [ "\n" <MD-N  array element list> ]
    <MD-N  array element>      ::= <index list> <N-tuple block>
    <index list>         ::= <integer> [ <index list>]     # note: here an M-dimensional array has M entries in the index list
  <N-tuple block>  ::= <real number> [ <N-tuple block> ]   #   this is just an "N-tuple", but without the enclosing parenthesis.

  # ==========================================================================================================
  # ==== other XML elements (e.g. "/<object>/<interface>/stats" should be supported on an as-needed basis ====
  # ==========================================================================================================

"""

# ------------------------ typed variables from XML "text": ------------------------------------------------
# ------------ all "histogram" type data is converted to "float" -------------------------------------------

def int_from_text(text):
  return int(text)

def real_from_text(text):
  return float(text)

def tuple_from_text(text):
  return eval(text)

def ND_indexed_array_from_text(text):
  istr = STR.StringIO(text)
  elements = []
  prev_tokens = None
  for line in istr:
    tokens = line.split()
    if len(tokens) > 0:
      elements.append(float(tokens[-1]))
      prev_tokens = tokens
  istr.close()
  
  rval = None
  if (len(elements) > 0):
    # shape is bin counts at last line containing tokens:
    shape = [int(prev_tokens[n])+1  for n in range(len(prev_tokens)-1)]    
    rval = array(elements).reshape(shape)  
  
  return rval

def ND_row_major_array_from_text(text, shape = None):
  # unspecified shape => 2D array, where each line in the text is an array row:
  nrows, ncols = 0, 0
  
  istr = STR.StringIO(text)  
  elements = []
  
  for line in istr:
    tokens = line.split()
    if (len(tokens) == 0):
      continue;
      
    if (ncols == 0):
      ncols = len(tokens)
    elif (ncols != len(tokens)):
      raise  RuntimeError, 'for ND row-major array, length of each line in \"text\" must match'   
    for t in tokens:
      elements.append(float(t))
    nrows += 1      
  istr.close()
  
  if not shape:
    shape = (nrows, ncols)
    
  return array(elements).reshape(shape)


def MD_Ntuple_array_from_text(text, N):
  # <bin number list>  <raw n-tuple>
  #   <bin number list> ::=  <integer> [ <bin number list> ] #  M entries
  #   <raw n-tuple>     ::=  <number list>
  #   <number list>     ::= <number>  [ <number list> ]      #  N entries
  
  istr = STR.StringIO(text)
  elements = []
  prev_tokens = None
  
  for line in istr:
    tokens = line.split()
    if len(tokens) == (N + 1):
      elements.extend([float(tokens[n]) for n in range(-N,0)])
      prev_tokens = tokens
  istr.close()

  # compute shape from last complete line:
  shape = [int(prev_tokens[n])+1 for n in range(len(prev_tokens)-N)]
  shape.append(N)  # last dimension is n-tuple dimension
  
  return array(elements).reshape(shape)

# ----------------------------------------------------------------------------------------------------------

def ignore_multiple_elements(elements, key):
  """
   This bypasses a Pol-MC output bug, in which objects output interfaces that they do not own (e.g. for interior objects)
     (this will allow us to use data XML files generated prior to fixing this bug).
     
   return an iterable which includes only the first elements with distinct values for the attribute
  """
  
  unique = {}
  for elt in elements:
    if not (elt.get(key) in unique):
      unique[elt.get(key)] = elt
  return unique.itervalues()

  
def optional_value(dict, key):
  # return the value if the dict exists and the key is present, otherwise return None
  rval = None
  if ((dict != None) and (key in dict)):
    rval = dict[key]
  return rval

def convert_from_element(key, element, convertor, accumulator = None, convertor_parms = None):
  """
  convert from ElementTree element to arbitrary type
  input parms:
    key: sub-element name
    element: source element to convert
    convertor: conversion method
    accumulator: previous value to accumulate to
    convertor_parms: additional parameters for convertor
  notes:
    specified accumulator makes presence of sub-element mandatory, otherwise "None" will be returned if sub-element is not present
    convertor signature: <arbitrary type> = convertor(element, accumulator = None, convertor_parms = None)
  """
  rval = None
  sub_element = element.find(key)

  if (accumulator != None): # test for presence, not value
    if (sub_element == None):
      raise RuntimeError, 'convert_from_element: accumulation: required sub-element \"%s\" is not present' % key
    rval = convertor(sub_element, accumulator, convertor_parms)
  else:
    if (sub_element != None):  # for some reason "if <element>:" doesn't work as expected with "ElementTree.element"
      rval = convertor(sub_element, None, convertor_parms)
      
  return rval  
  
def text_from_element(element, accumulator = None, convertor_parms = None):
  if (accumulator != None):  # test for presence, not value
    raise RuntimeError, 'text_from_element: \"accumulator\" should be None for text'
  return element.text
  
def int_from_element(element, accumulator = None, convertor_parms = None):
  rval = 0
  if (accumulator != None):  # test for presence, not value
    rval = accumulator
    rval += int(element.text)
  else:
    rval = int(element.text)
  
  return rval
  
def float_from_element(element, accumulator = None, convertor_parms = None):
  rval = 0
  if (accumulator != None):  # test for presence, not value
    rval = accumulator
    rval += float(element.text)
  else:
    rval = float(element.text)
  
  return rval

def tuple_from_element(element, accumulator = None, convertor_parms = None):
  rval = 0
  if (accumulator != None):  # test for presence, not value
    raise RuntimeError, 'tuple_from_element: \"accumulator\" should be None for tuple'
  else:
    rval = eval(element.text)
  
  return rval

def ND_indexed_array_from_element(element, accumulator = None, convertor_parms = None):
  istr = STR.StringIO(element.text)
  elements = []
  prev_tokens = None
  for line in istr:
    tokens = line.split()
    if len(tokens) > 0:
      elements.append(float(tokens[-1]))
      prev_tokens = tokens
  istr.close()
  
  rval = None
  if (len(elements) > 0):
    # shape is bin counts at last line containing tokens:
    shape = tuple([int(prev_tokens[n])+1  for n in range(len(prev_tokens)-1)]) # compare "tuple" shapes only 
    if (accumulator != None):  # test for presence, not value
      if (shape != accumulator.shape):
        raise RuntimeError, 'ND_indexed_array_from_element: \"accumulator\" shape does not match that of array from element'
      rval = accumulator
      rval +=  array(elements).reshape(shape)
    else: 
      rval = array(elements).reshape(shape)  
  
  return rval


def ND_row_major_array_from_element(element, accumulator = None, shape = None):
  # shape takes "convertor_parms" argument position
  # unspecified shape => 2D array, where each line in the text is an array row:
  nrows, ncols = 0, 0
  
  istr = STR.StringIO(element.text)  
  elements = []
  
  for line in istr:
    tokens = line.split()
    if (len(tokens) == 0):
      continue;
      
    if (ncols == 0):
      ncols = len(tokens)
    elif (ncols != len(tokens)):
      raise  RuntimeError, 'for ND row-major array, length of each line in \"text\" must match'   
    for t in tokens:
      elements.append(float(t))
    nrows += 1      
  istr.close()
  
  if not shape:
    shape = (nrows, ncols) # compare "tuple" shapes only
  
  rval = None
  if (len(elements) > 0):
    if (accumulator != None):  # test for presence, not value
      if (shape != accumulator.shape):
        # print 'shape %s, accu shape %s' % (shape, accumulator.shape)
        raise RuntimeError, 'ND_row_major_array_from_element: \"accumulator\" shape does not match that of array from element'
      rval = accumulator
      rval +=  array(elements).reshape(shape)
    else: 
      rval = array(elements).reshape(shape)  
  
  return rval


def MD_Ntuple_array_from_element(element, accumulator = None, N = None):
  # <bin number list>  <raw n-tuple>
  #   <bin number list> ::=  <integer> [ <bin number list> ] #  M entries
  #   <raw n-tuple>     ::=  <number list>
  #   <number list>     ::= <number>  [ <number list> ]      #  N entries
  # N takes "convertor_parms" argument position
  # unspecified N => 2 element tuples
  
  if not N:
    N = 2
    
  istr = STR.StringIO(element.text)
  elements = []
  prev_tokens = None
  
  for line in istr:
    tokens = line.split()
    if len(tokens) == (N + 1):
      elements.extend([float(tokens[n]) for n in range(-N,0)])
      prev_tokens = tokens
  istr.close()

  # compute shape from last complete line:
  shape = [int(prev_tokens[n])+1 for n in range(len(prev_tokens)-N)]
  shape.append(N)  # last dimension is n-tuple dimension
  shape = tuple(shape) # compare "tuple" shapes only
  
  rval = None
  if (len(elements) > 0):
    if (accumulator != None):  # test for presence, not value
      if (shape != accumulator.shape):
        raise RuntimeError, 'MD_Ntuple_array_from_element: \"accumulator\" shape does not match that of array from element'
      rval = accumulator
      rval +=  array(elements).reshape(shape)
    else: 
      rval = array(elements).reshape(shape)  
  
  return rval  

class scattering_stats(object):
  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator == None):  # test for presence, not value
      data = None
    else:
      data = accumulator  # return value
    
    if (element.get('type') == '1D'):
      data = MD_Ntuple_array_from_element(element, accumulator, 2)
    elif (element.get('type') == '2D'):
      data = ND_row_major_array_from_element(element, accumulator) # default is 2D array   
    else:
      raise RuntimeError, 'scattering_stats: convert_from_element: unknown value for \"type\" attribute: %s' % element.get('type')
    
    return data
    
    
class Stokes_data(object):
  @staticmethod
  def zeros(minAcceptanceCosine, maxAcceptanceCosine, shape = None):
    # return a "Stokes_data" object of specified shape (None => scalar) for accumulation:
    data = {'minAcceptanceCosine':minAcceptanceCosine, 'maxAcceptanceCosine':maxAcceptanceCosine}
    if (shape != None):
      data['I'], data['Q'], data['U'], data['V'], data['betalin'], data['betacirc'] = \
        zeros(shape,dtype=double), zeros(shape,dtype=double), zeros(shape,dtype=double), \
        zeros(shape,dtype=double), zeros(shape,dtype=double), zeros(shape,dtype=double)
    else:
      data['I'], data['Q'], data['U'], data['V'], data['betalin'], data['betacirc'] = \
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    return data

  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator == None):  # test for presence, not value
      data = {} # return value dict
    else:
      data = accumulator
    
    # accumulate from an attribute (as opposed to an element) value:
    if (accumulator == None):
      data['transmittance'] = element.get('transmittance')
    else:    
      if 'transmittance' in accumulator:
        data['transmittance'] += element.get('transmittance')
      else:
        raise RuntimeError, 'Stokes_data: accumulate: \"transmittance\" attribute not present in element'
    
    data['minAcceptanceCosine'] = convert_from_element('minAcceptanceCosine', element, float_from_element) # no accumulator
    data['maxAcceptanceCosine'] = convert_from_element('maxAcceptanceCosine', element, float_from_element)
    
    # if (optional_value(accumulator, 'I') != None):
    #   print '%s <-- ' % (accumulator['I'].shape,),
    
    data['I'] = convert_from_element('I', element, ND_row_major_array_from_element, optional_value(accumulator,'I'))
    
    # if (optional_value(data, 'I') != None):
    #   print '%s' % (data['I'].shape,),
    
    data['Q'] = convert_from_element('Q', element, ND_row_major_array_from_element, optional_value(accumulator,'Q'))
    data['U'] = convert_from_element('U', element, ND_row_major_array_from_element, optional_value(accumulator,'U'))
    data['V'] = convert_from_element('V', element, ND_row_major_array_from_element, optional_value(accumulator,'V'))
    
    # recalculate degree of linear and circular polarization:
    data['betalin']  = zeros_like(data['I'])
    data['betacirc'] = zeros_like(data['I'])

    for nd_x in ndindex(data['I'].shape):     
      if (data['I'][nd_x] != 0.0):
        data['betalin'][nd_x]  = sqrt(data['Q'][nd_x]**2 + data['U'][nd_x]**2) / data['I'][nd_x]
        data['betacirc'][nd_x] = abs(data['V'][nd_x])/data['I'][nd_x]

    return data
              
class scatterer(object):
  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator != None):  # test for presence, not value
      data = accumulator # return value dict
    else:
      data = {}

    if (accumulator == None):
      data['scatteringStats'] = {}
    for S in element.findall('scatteringStats'):
      data['scatteringStats'][S.get('type')] = scattering_stats.convert_from_element(S, optional_value(
                                                                                          optional_value(accumulator, 'scatteringStats'), S.get('type')))
    
    data['scattererType'] = convert_from_element('scattererType', element, text_from_element)
    data['anisotropy'] = convert_from_element('anisotropy', element, float_from_element)
    data['indexScatterer'] = convert_from_element('indexScatterer', element, float_from_element)
    data['scattererRadius'] = convert_from_element('scattererRadius', element, float_from_element)
    data['mu_a'] = convert_from_element('mu_a', element, float_from_element)
    data['mu_s'] = convert_from_element('mu_s', element, float_from_element)
    data['index'] = convert_from_element('index', element, float_from_element)
     
    if (accumulator == None):
      data['sampledAnisotropy'] = convert_from_element('sampledAnisotropy', element, float_from_element)
      if (isnan(data['sampledAnisotropy'])):
        data['sampledAnisotropy'] = None  # "pickle" can't deal with NAN ( play nice... )
    elif (('sampledAnisotropy' in accumulator) 
              and (accumulator['sampledAnisotropy'] != None)
              and ('1D' in data['scatteringStats'])):  # note: 'sampledAnisotropy' may be NAN if there's no scatteringStats['1D'] at all...
      # recalculate sampled anisotropy from the "1D" scattering stats:
      try:
        stats = data['scatteringStats']['1D']
        vt_W = linspace(0.0, pi, stats.shape[0])  # 1D stats is (N_bins, 2) shape
        vt_W *= cos(stats[:,0]) # cos(\theta) weighted by bin counts
        data['sampledAnisotropy'] = vt_W.sum() / float((stats[:,0]).sum())
      except:
        raise RuntimeError, 'scattering_stats: accumulate: unable to recalculate \"sampledAnisotropy\"'  
    
    return data 

class interface(object):
  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator == None):  # test for presence, not value
      data = {} # return value dict
    else:
      data = accumulator

    # print '   interface: %s' % element.get('name')
    
    data['id'] = int(element.get('id'))
    data['normal'] = convert_from_element('normal', element, tuple_from_element) # no accumulator
    data['origin'] = convert_from_element('origin', element, tuple_from_element)
    
    if (accumulator == None):
      data['basevector'] = {}
      for V in element.findall('basevector'):
        data['basevector'][V.get('id')] = tuple_from_text(V.text)

    if (accumulator == None):
      data['detectionbasis'] = {}
      for V in element.findall('detectionbasis'):
        data['detectionbasis'][V.get('id')] = tuple_from_text(V.text)
      
    if (accumulator == None):
      data['StokesV'] = []
    for n, SV in enumerate(element.findall('StokesV')):
    
      # print '      -- StokesV[%d] -- ' % n,
    
      if (accumulator == None):
        stokes_data = None
      else:
        stokes_data = accumulator['StokesV'][n]
                
      stokes_data = Stokes_data.convert_from_element(SV, stokes_data)
      
      # print ''
      
      if (accumulator == None):
        data['StokesV'].append(stokes_data)
      else:
        data['StokesV'][n] = stokes_data  # probably redundant (i.e. data has already been accumulated to correct destination)

    if (accumulator == None):
      data['scatteringStats'] = {}
    for S in element.findall('scatteringStats'):
      data['scatteringStats'][S.get('type')] \
        = scattering_stats.convert_from_element(S, optional_value(optional_value(accumulator, 'scatteringStats'), S.get('type')))    

    return data
      
    
class object_data(object):
  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator == None):  # test for presence, not value
      data = {} # return value dict
    else:
      data = accumulator

    # print '------ object: %s ----------' % element.get('name')
    
    if (accumulator == None):
      data['interfaces'] = {}
      
    for I in ignore_multiple_elements(element.findall('interface'), 'name'):
      data['interfaces'][I.get('name')] \
        = interface.convert_from_element(I, optional_value(optional_value(accumulator, 'interfaces'), I.get('name')))

    data['scattering'] = scatterer.convert_from_element(element, optional_value(accumulator, 'scattering')) 
    
    data['id'] = int(element.get('id'))
    data['origin'] = convert_from_element('origin', element, tuple_from_element)         # no accumulator
    data['absorbance'] = convert_from_element('absorbance', element, float_from_element) # no accumulator  
    
    return data     

class world(object):
  def __init__(self):
    self._data = {}
    self.initialized = False
  
  def __getitem__(self, key):
    return self._data[key]
      
  def read_from_XML(self, fp_or_filename, accumulate=False):
    self._data = world.convert_from_element(ET.parse(fp_or_filename), accumulate and self._data)
    self.initialized = True
    
  @staticmethod
  def convert_from_element(element, accumulator = None):
    if (accumulator != None):  # test for presence, not value
      data = accumulator # return value dict
    else:
      data = {}

    data['energyBinNumber'] = convert_from_element('/energyBinNumber', element, 
                                                    ND_indexed_array_from_element,                          # convertor
                                                    optional_value(accumulator, 'energyBinNumber'))         # accumulator

    data['fluenceBinNumber'] = convert_from_element('/fluenceBinNumber', element, 
                                                    ND_indexed_array_from_element,           
                                                    optional_value(accumulator, 'fluenceBinNumber')) 
    if (not 'objects' in data):
      data['objects'] = {}
    for O in element.findall('/object'):
      data['objects'][O.get('name')] = object_data.convert_from_element(O, optional_value(optional_value(accumulator, 'objects'), O.get('name')))
    
    return data
    
    
