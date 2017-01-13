"""
  Classes and utility methods associated with binary format photon path files as produced by Pol-MC.
  
  A photon_path_file, as output in "kBinaryFormat64" format from Pol-MC, is presented as a list of photon_path lists.
  Each photon_path is a list of photon_event.
"""

import os, sys, copy, pdb

from numpy import *

from binaryIO import *
set_IO_sizeBits(32)

# global constant:
PHOTON_SAFETY_DISTANCE = 5.0e-8
OBJECT_SAFETY_DISTANCE = PHOTON_SAFETY_DISTANCE/2.0



def euler_rotation_matrix(alpha, beta, gamma):
  """
    generate Euler rotation matrix
    to rotate cartesian coordinates following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
  """
  ca = cos(alpha); sa = sin(alpha)
  cb = cos(beta); sb = sin(beta)
  cg = cos(gamma); sg = sin(gamma)
  aRot = matrix([[cg*cb*ca - sg*sa, cg*cb*sa + sg*ca, -cg*sb],
                 [-sg*cb*ca - cg*sa, -sg*cb*sa + cg*ca, -sg*sb],
                 [sb*ca, sb*sa, cb]])
  return aRot
  

def inverse_euler_rotation_matrix(alpha, beta, gamma):
  """
  generate inverse Euler rotation matrix
    to produce the inverse of the cartesian coordinate rotation following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
    (for the _forward_ transformation)
  """
  ca = cos(alpha); sa = sin(alpha)
  cb = cos(beta); sb = sin(beta)
  cg = cos(gamma); sg = sin(gamma)
  aInvRot = matrix([[cg*cb*ca - sg*sa, -sg*cb*ca - cg*sa, sb*ca],
                 [cg*cb*sa + sg*ca, -sg*cb*sa + cg*ca, sb*sa],
                 [-cg*sb, -sg*sb, cb]])
  return aInvRot


def normalized_dot_product(v1, v2):
  if not isinstance(v1, ndarray):
    v1 = array(v1)
  if not isinstance(v2, ndarray):  
    v2 = array(v2)
  n1 = linalg.norm(v1)
  n2 = linalg.norm(v2)  
  v = v1.copy()
  v *= v2
  rval = v.sum()
  rval /= n1 * n2
  return rval

class StokesV(object):
  """
    Stoke's vector
    attributes:
      I, Q, U, V:    Stoke's parameters
      E_p, E_s, E_0: Electric field vectors in parallel, perpendicular, and propagation direction, respectively,
      phase:         Total propagation phase
  """
  __slots__ = ['I','Q','U','V','E_p','E_s','E_0','phase']
  
  def __init__(self, I = 1.0, Q = 0.0, U = 0.0, V = 0.0, E_p = (1.0, 0.0, 0.0), E_s = (0.0, 1.0, 0.0), E_0 = (0.0, 0.0, 1.0), phase = 0.0):
    self.I = I
    self.Q = Q
    self.U = U
    self.V = V
    self.E_p = E_p
    self.E_s = E_s
    self.E_0 = E_0
    self.phase = phase
  
  def read_binary(self, fp):
    self.I, self.Q, self.U, self.V = load_double(fp), load_double(fp), load_double(fp), load_double(fp)
    self.E_p, self.E_s, self.E_0 = loadntuple3Double(fp), loadntuple3Double(fp), loadntuple3Double(fp)
    self.phase = load_double(fp)
  
  def normalized(self):
    # factory method: create a normalized Stoke's vector from the current object.
    rval = None
    if (self.I > finfo(double).eps):
      rval = StokesV(1.0, self.Q/self.I, self.U/self.I, self.V/self.I, self.E_p, self.E_s, self.E_0, self.phase)
    else:
      print 'photon_path: StokesV: normalized: *** WARNING: Stoke\'s vector is not normalizable: (I, Q, U, V): (%f, %f, %f, %f) ****'\
        % (self.I, self.Q, self.U, self.V)
        
    return rval
    
  def inLabFrame(self, e_p, e_s, e_0):
    """
    modify the reference frame of the Stoke's vector to be consistent with the specified surface e_p, e_s, and normal vectors.
      input parms:
        e_p: reference direction (at surface) for parallel polarization
        e_s: reference direction (at surface) for perpendicular polarization
        e_0: surface normal
        
      NOTE: If you rotate the surface, you also need to rotate the basis vectors that you give to this method.
    """
    
    """
    Implementation note:
    
      1) A complex 2D electric field vector in the photon propagation frame
        $ E = (E_1 \exp(i \delta_1), E_2 \exp(i \delta_2) $\ 
        may be represented in terms of its Stoke's parameters:
        \begin{align*}
          I =&  |E|^2 \\
          Q =&  |E_1|^2 - |E_2|^2 \\
          U =&  \Re(E_2 \cdot E_1^*) \\
            =&  \frac{1}{2} ( E_2 \cdot E_1^* + E_2^* \cdot E_1 )
          V =&  =&  \Im(E_2 \cdot E_1^*) \\
            =&  \frac{1}{2} ( E_2 \cdot E_1^* - E_2^* \cdot E_1 )
        \end{align*}
        
      2) This vector is then rotated using the appropriate Euler rotation matrix into reference frame of the surface;
      
      3) The projection in the surface is taken;
      
      4) The projection is reconverted into the corresponding Stoke's parameters. 
      
      Note that the operation  of the rotation of a vector is the inverse to that of the rotation of its co-ordinate system, 
        and for this reason, if we compute the required angles to rotate the photon reference frame into that of the surface, 
        and then build the inverse Euler rotation matrix. 
    """
    # compute complex field vector corresponding to Stoke's parameters:
    delta = arctan2(self.V, self.U)  # $ \delta_y - \delta_x $\ W.L.O.G. take $\delta_x \equiv 0$
    E_1 = sqrt((self.I + self.Q)/2.0)
    E_2 = sqrt((self.I - self.Q)/2.0)
    E = (E_1, E_2 * exp(1j * delta), 0.0)
    
    # compute Euler angles for rotating the photon reference frame into the lab reference frame (using the Z, X', Z'' convention):
    
    # components of  photon propagation vector in lab frame:
    n1, n2, n3 = normalized_dot_product(self.E_0, e_p), normalized_dot_product(self.E_0, e_s), normalized_dot_product(self.E_0, normal)
    
    alpha = -arctan2(n1,n2) #  Z rotation (rotate s.t. the normal has no x component)
    beta  = -arccos(n3)   #  X' rotation (rotate s.t. the normal has no y component)
    gamma = 0.0
    
    aInvRot = inverse_euler_rotation_matrix(alpha, beta, gamma)
    # rotate into the lab frame:
    E_lab = aInvRot * matrix(E).T
    # take the projection onto the xy-plane:
    E_lab[2] = 0.0
    E_lab2 = multiply(E_lab, conj(E_lab))  # (abs(E_1)**2, abs(E_2)**2)
     
    # recalculate the Stoke's parameters:
    self.I = dot(E_lab.T, conj(E_lab))
    self.Q = E_lab2[0] - E_lab2[1]
    self.U = real(E_lab[1]*conj(E_lab[0]))
    self.V = imag(E_lab[1]*conj(E_lab[0]))
    self.E_p, self.E_s, self.E_0 = e_p, e_s, normal  # copy from the lab frame
     
    
class photon_event(object):
  """ 
    Information associated with a Monte Carlo "photon" scattering event.
    attributes:
      pos: position tuple (3D)
      dir: direction unit vector      
  """
  __slots__ = ['pos','dir','Stokes','W','t','event_type']
  
  def __init__(self, pos = (0.0, 0.0, 0.0), dir = (0.0, 0.0, 1.0), Stokes = StokesV(), W = 0.0, t = 0.0, event_type = 0):
    self.pos = pos
    self.dir = dir
    self.Stokes = Stokes
    self.W = W
    self.t = t
    self.event_type = event_type
    
  def read_binary(self, fp):
    self.pos, self.dir = loadntuple3Double(fp), loadntuple3Double(fp)
    self.Stokes.read_binary(fp)
    self.W, self.t = load_double(fp), load_double(fp)
    self.event_type = load_int32_t(fp)  

def photon_path_read_binary(fp):
  photon_path = None
  
  hdr = fp.read(len('PHOTON'))
    
  if not hdr:
    return # OK: probably end-of-file
  elif (hdr != 'PHOTON'):
    raise IOError, 'photon_path_read_binary: header not found: file format not \"kBinaryFormat64\"'
      
  N_events = load_size_t(fp) # note: "set_IO_sizeBits(32)" at module load
  try:
    if (N_events > 0):
      photon_path = []
    
    for nevent in range(N_events):
      # create a distinct event object for each event:
      event = photon_event()    
      event.read_binary(fp)
      photon_path.append(event) 
  except:
    raise IOError, 'photon_path_read_binary: file format error, or incomplete photon-path record'
  
  return photon_path
  
  
def load_photon_paths(input_filename):
  photon_paths = []
  fp = open(input_filename, 'rb')
  
  while True:
    path = photon_path_read_binary(fp)
    if path:
      photon_paths.append(path)
    else:
      break 
      
  fp.close()
  return photon_paths


class photon_paths_iter(object):
  def __init__(self, fp_or_filename):
    if isinstance(fp_or_filename, file):
      self.fp = fp_or_filename
      self.close_fp = False
    else:
      self.fp = open(fp_or_filename, 'rb')
      self.close_fp = True
      
  def __del__(self):
    if self.close_fp:
      self.fp.close()
    self.close_fp = False
              
  def __iter__(self):
    return self
    
  def next(self):
    path = photon_path_read_binary(self.fp)
    if not path:
      raise StopIteration
    return path
       
  def reset(self):
    self.fp.seek(0)

"""
# ================================ obsolete methods: ==================================================
# ============= Implementation note: "interface_bins" method requires knowledge of non-scattering "event_type";
# =============   this indicated that actually removing non-scattering events from a photon-path wasn't
# =============   a good idea.  Methods are now implemented to _ignore_ events not matching the specified type.

def select_event_by_event_type(photon_event, event_type):
  # photon_event if event is of specified type:
  rval = None
  if photon_event.event_type == event_type:
    rval = photon
  return rval
  
def restrict_path_by_event_type(photon_path, event_type):
  # photon path including only events of specified type:
  rval = []
  for event in photon_path:
    if event.event_type == event_type:
      rval.append(event)
  return rval  

def restrict_paths_by_event_type(photon_paths, event_type):
  # photon-paths including only events of specified type:
  rval = []
  for photon in photon_paths:
    restrict_photon = restrict_path_by_event_type(photon, event_type)
    if (len(restrict_photon) > 0):
      rval.append(restrict_photon)
  
  return rval
# ================================== end: obsolete methods =============================================
"""


def photon_depth_classes(photon_paths, depth_ranges, event_type = 1):
  """
    separate photon paths by maximum vertical depth
    (note: implementation of the general case with rotation is easy).
    input parms:
      paths: photon-paths iterable
      depth_ranges: list of (<z min>, <z max>) tuples
      event_type: 
  """
  rval = [[] for n in range(len(depth_ranges))]
  
  for photon in photon_paths:
    # get depth-range for the path:
    z_min, z_max = finfo(double).max, finfo(double).min
    any_events = False
    for event in photon:
      if (event.event_type != event_type):  # ignore photon_events not of specified event_type
        continue
      any_events = True
      z_min = min(z_min, event.pos[2])
      z_max = max(z_max, event.pos[2]) # here we would project along a normal for the general case
    
    # classify path by its maximum depth:
    if any_events:
      for nr, rtuple in enumerate(depth_ranges):
        if ((z_max > rtuple[0] + OBJECT_SAFETY_DISTANCE) and (z_max < rtuple[1] - OBJECT_SAFETY_DISTANCE)):
          rval[nr].append(photon)
        
  return rval
    
def volume_bins(photon_paths, bin_range, shape, 
                bin_type = 'fluence', smooth = False, mu_a = 1.0, mu_s = 0.0, 
                event_type = 1, detector_normal = (0.0, 0.0, 1.0), minAcceptanceCosine = 0.0, maxAcceptanceCosine = 1.0):
  """
  accumulate volume bins from photon paths
  note: paths should _already_ include only photon events of desired type (use "restrict" methods above).
  
  input parms:
    
    photon_paths: iterable
    
    bin_range: list containing range tuples for each dimension: [(min_x, max_x), (min_y, max_y), (min_z, max_z),...]
    
    shape: shape of the bins array (number of dimensions corresponds to tuple entries in "bin_range")
    
    bin_type: value to accumulate to the bins: 'fluence' or 'energy'
    
    smooth: also accumulate at inter-event positions using interpolation with specified mu_a
    
    mu_a: mu_a to use for scaling energy, or for interpolation
    
    mu_s: mu_s to use for scaling energy
   
    event_type: type of photon_event to add to fluence or energy accumulation
    
    detector_normal: normal vector to the detector surface
    
    minAcceptanceCosine, maxAcceptanceCosine: where antiparallel to detector normal is 1.0
    
      NOTES: 
        -- energy deposition is scaled by mu_a / mu_t, so default values return non-rescaled energy (i.e. which is the same as fluence);
        -- fluence deposition is as unscaled photon weight, which from sampling considerations, 
             should give fluence in units of "mu_t" (times whatever the units of photon weight are); 
             that is , fluence should be rescaled by "1.0/mu_t" to get physical units.
        .
  returns:
    array of bins  
  """
  
  def bin_index(pos, bin_range, shape):
    # calculate ND index mapping position pos into an ND-array with specified coordinate ranges and shape.
    #   notes: 
    #     -- bin_range tuples shall be ordered
    #     -- out-of-range positions return "None"
    
    nx = [int(floor((pos[n] - bin_range[n][0])/(bin_range[n][1] - bin_range[n][0]) * float(shape[n]))) for n in range(len(shape))]

    if (any([n < 0 for n in nx]) or any([nx[n] >= shape[n] for n in range(len(nx))])):
      nx = None
    else:
      nx = tuple(nx) # support index "by-tuple" (i.e. for many containers: "by-list" indexing doesn't work)
      
    return nx
  
  if (bin_type == 'fluence'):
    bFluence = True
  elif (bin_type == 'energy'):
    bFluence = False
  else:
    raise RuntimeError, 'photon_path: volume_bins: unknown accumulation type %s requested' % bin_type
  
  # ordered acceptance-cosine tuple:
  acceptanceCosine = (min(minAcceptanceCosine, maxAcceptanceCosine), max(minAcceptanceCosine, maxAcceptanceCosine))
  
  # extinction coefficient:
  mu_t = mu_a + mu_s 
    
  bins = zeros(shape,dtype=float) # return value
  
  # pdb.set_trace()
  
  for photon in photon_paths:
    # ignore photons outside of the acceptance-angle:
    #   cosign at the detector surface (s.t. antiparallel is 1.0):
    cosine = -normalized_dot_product(photon[-1].dir, detector_normal)
    if ((acceptanceCosine[0] > cosine) or (cosine > acceptanceCosine[1])):
      continue;      
    for event in photon:
      if (event.event_type != event_type):  # ignore photon events not of the specified type
        continue
      if not smooth:
        nx = bin_index(event.pos, bin_range, shape)
        if (nx != None):
          # note: index will be None if pos is out-of-range, 
          #   which is not an error here (bins are accumulated only for the region of interest specified via bin_range).
          if bFluence:
            bins[nx] += event.W
          else:
            bins[nx] += event.W * (mu_a / mu_t)
      else:
        pass
  
  return bins

def detector_fluence(photon_paths, e_p, e_s, e_0, minAcceptanceCosine, maxAcceptanceCosine, N_cosine_bins):
  """
    Total fluence at the detector from photon paths: same format as in polmc_data (w.r.t. accumulated cosine-bins, except only one spatial bin).
      Last two photon_events in the path shall correspond to the interface transition, photon weight _prior_ to this transition
        is used as the value for accumulation.

    input parms:
      photon_paths: iterable
      
      e_p, e_s, e_0: surface basis vectors (i.e. "el", "er", and "normal", using Pol-MC's notation)
        NOTE: if you rotate the detector surface, these basis vectors should be rotated as well (using the same Euler rotation matrix);
          WARNING: Pol-MC uses the X, Y', Z'' convention (which is very rare and strange for Euler rotations);
            this method uses the Z, X', Z'' convention.  
            For example Pol-MC (0.0, pi/4.0, 0.0) is (-pi/2.0, pi/4.0, 0.0) using the Z, X', Z'' convention.
             
      minAcceptanceCosine, maxAcceptanceCosine:       
        range of photon-event direction-cosines (with respect to surface normal) corresponding to accumulated events;
          note: sign of dot-product is reversed (i.e. anti-parallel photon_event direction and surface normal correspond to cosine = 1.0).
          
      N_cosine_bins: number of accumulation bins.
    
    return value:
      list of cosine bins (in Pol-MC, "total increasing angle" format) of Stokes vector data in "Stokes_data" format (i.e. same as "polmc_data" object)
      
    implementation notes:
      -- this is a detector-specific implementation, it could also be generalized to calculate interface fluence;
      -- as the value for detection-type "mDetect" is unknown here, "minAcceptanceCosine" and "maxAcceptanceCosine" must be set correctly 
           with respect to the surface normal (note: _sign_ of dot-product is inverted for Pol-MC): for example: 
             to accumulate all entry-events: (minAcceptanceCosine, maxAcceptanceCosine) = (0.0,  1.0);
             to accumulate all exit-events:  (minAcceptanceCosine, maxAcceptanceCosine) = (0.0, -1.0);
           there is no requirement that (minAcceptanceCosine <= maxAcceptanceCosine).
      -- removing photon-events that are _not_ scattering events will screw-up INTERFACE FLUENCE values
         (see note above, with respect to "obsolete methods").
      .
  """
  cosine_bin_edges = linspace(minAcceptanceCosine, maxAcceptanceCosine, N_cosine_bins)
  bins = [polmc_data.Stokes_data.zeros(cosine_bin_edges[n-1], cosine_bin_edges[n], shape = None)  for n in range(1, N_cosine_bins)]
  
  def cosine_bin_index(cosine, bin_range, N_bins):
    n = int(floor((cosine - bin_range[0])/(bin_range[1] - bin_range[0]) * float(N_bins)))
    if (n >= N_cosine_bins):
      n = None
    return n
    
  for photon in photon_paths:
    event = photon[-1]
    if (event.event_type != 3):  # kInterfaceTransmitEvent: note: photon statistics are recorded _prior_ to killing the photon.
      raise RuntimeError, 'photon_path: detector_fluence: last event in photon_path is not \"kInterfaceTransmitEvent\"'
    
    cosine = -normalized_dot_product(event.dir, normal)  # "cosine < 0" => exit-events, "cosine > 0" => entry-events
    localStokes = event.Stokes.normalized()
    localStokes.toLabFrame(e_p, e_s, e_0)
    
    n = cosine_bin_index(cosine, (minAcceptanceCosine, maxAcceptanceCosine), N_cosine_bins)
    if (n != None):
      # note: index will be None if cosine is out-of-range,
      #   which is not an error here (bins are accumulated only for the region of interest specified via (minAcceptanceCosine, maxAcceptanceCosine)).
      bins[n]['I'] += localStokes.I * event.W 
      bins[n]['Q'] += localStokes.Q * event.W
      bins[n]['U'] += localStokes.U * event.W
      bins[n]['V'] += localStokes.V * event.W
  
  # accumulate cosine bins to Pol-MC "expanding angle" format:
  for n in range(N_cosine_bins-1,0,-1):
    bins[n-1]['I'] += bins[n]['I']
    bins[n-1]['Q'] += bins[n]['Q']
    bins[n-1]['U'] += bins[n]['U']
    bins[n-1]['V'] += bins[n]['V']
    bins[n-1]['maxAcceptanceCosine'] = 1.0 
     
  return bins
  
      
