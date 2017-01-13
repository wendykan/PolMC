"""
  Classes and utility methods associated with binary format photon path files as produced by Pol-MC.
  
  A photon_path_file, as output in "kBinaryFormat64" format from Pol-MC, is presented as a list of photon_path lists.
  Each photon_path is a list of photon_event.
"""

import os, sys, pdb
import copy as CP

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

# -----------------------------------------------------------------------------------------------------------------------------------
# ------------------------------- clones of Pol-MC functions of the same name: ------------------------------------------------------
# ----------- KAT 31.07.2010: these functions should match the implementation of the Pol-MC functions -------------------------------
# -----------   NO CORRECTIONS are imposed if errors were found in Pol-MC (e.g. "MeasureStokesVectorInLabFrame") --------------------
# -----------------------------------------------------------------------------------------------------------------------------------

def normalized_dot_product(v1, v2):
  if not isinstance(v1, ndarray):
    v1 = array(v1)
  if not isinstance(v2, ndarray):  
    v2 = array(v2)
  rval = dot(v1, v2)
  rval /= sqrt(dot(v1,v1) * dot(v2,v2))
  return rval


def normalized_cross_product(v1, v2):
  if not isinstance(v1, ndarray):
    v1 = array(v1)
  if not isinstance(v2, ndarray):  
    v2 = array(v2)
  rval = cross(v1, v2)
  rval /= sqrt(dot(v1,v1) * dot(v2,v2))
  return rval


def oriented_angle_between(u, v, w):
  if not isinstance(u, ndarray):
    u = array(u)
  if not isinstance(v, ndarray):
    v = array(v)
  if not isinstance(w, ndarray):
    w = array(w)
  
  vSin = cross(u, v)
  phi = arctan2(linalg.norm(vSin), dot(u, v))

  """  
  sinPhi = normalized_cross_product(u, v)
    
  phi = arcsin(linalg.norm(sinPhi))
  
  if dot(u, v) <= 0.0:
    phi = pi - phi
  """
    
  if dot(vSin, w) <= 0.0:
    phi *= -1.0
    
  return phi

def intensity_through_linear_polarizer(S, e_s_lab, e_p_lab, e_0_lab):
  """
  I_p, I_s = intensity_through_linear_polarizer(<StokesV>, e_s_lab, e_p_lab, e_0_lab)
    warning: input vectors in CYCLIC ORDER (not in same order as Pol-MC method)!
  """
  e_s, e_p, e_0 = S.E_s, S.E_p, S.E_0
  vectorInParaNormalplane = normalized_cross_product(e_0, e_s_lab);
  vectorInPerpNormalplane = normalized_cross_product(e_p_lab, e_0);

  phi_para = oriented_angle_between(e_p, vectorInParaNormalplane, e_0);
  phi_perp = oriented_angle_between(e_p, vectorInPerpNormalplane, e_0);

  s = CP.deepcopy(S)
  s.rotate_ref_frame_around_e_0(phi_para)

  I_p = (s.I+s.Q)/2.0

  s = CP.deepcopy(S)
  s.rotate_ref_frame_around_e_0(phi_perp)

  I_s = (s.I+s.Q)/2.0  

  return I_s, I_p
  

def intensity_through_circular_polarizer(S, e_s_lab, e_p_lab, e_0_lab):
  """
  I_right, I_left = intensity_through_circular_polarizer(<StokesV>, e_s_lab, e_p_lab, e_0_lab)
    warning: input vectors in CYCLIC ORDER (not in same order as Pol-MC method)!
  """
  e_s, e_p, e_0 = S.E_s, S.E_p, S.E_0
  e_minus45 = (e_s_lab - e_p_lab)/sqrt(2.0)
  e_plus45 = (e_s_lab + e_p_lab)/sqrt(2.0)

  vectorInParaNormalplane = normalized_cross_product(e_0, e_s_lab);
  # vectorInPerpNormalplane = normalized_cross_product(e_p_lab, e_0);
  phi_para = oriented_angle_between(e_p, vectorInParaNormalplane, e_0);

  # In para-perp reference frame, apply right quarter waveplate
  s = CP.deepcopy(S)
  s.rotate_ref_frame_around_e_0(phi_para)
  s.U = S.V
  s.V = -S.U
  s.rotate_ref_frame_around_e_0(-phi_para)
  dummy, I_right = intensity_through_linear_polarizer(s, e_minus45, e_plus45, e_0_lab)
    
  # In para-perp reference frame, apply left quarter waveplate
  s = CP.deepcopy(S)
  s.rotate_ref_frame_around_e_0(phi_para)
  s.U = -S.V
  s.V = S.U
  s.rotate_ref_frame_around_e_0(-phi_para)
  dummy, I_left = intensity_through_linear_polarizer(s, e_minus45, e_plus45, e_0_lab)

  return I_left, I_right

# -----------------------------------------------------------------------------------------------------------------------------------

class StokesV(object):
  """
    Stoke's vector
    attributes:
      I, Q, U, V:    Stoke's parameters
      E_s, E_p, E_0: Electric field vectors in parallel, perpendicular, and propagation direction, respectively,
      phase:         Total propagation phase
  """
  __slots__ = ['I','Q','U','V','E_s','E_p','E_0','phase']
  
  def __init__(self, I = 1.0, Q = 0.0, U = 0.0, V = 0.0, E_s = (1.0, 0.0, 0.0), E_p = (0.0, 1.0, 0.0), E_0 = (0.0, 0.0, 1.0), phase = 0.0):
    self.I = I
    self.Q = Q
    self.U = U
    self.V = V
    self.E_s = E_s
    self.E_p = E_p
    self.E_0 = E_0
    self.phase = phase
    
  def __copy__(self): 
    return StokesV(self.I, self.Q, self.U, self.V, 
                   CP.copy(self.E_s), CP.copy(self.E_p), CP.copy(self.E_0),
                   self.phase)
  
  def __deepcopy__(self, memo):
    return StokesV(self.I, self.Q, self.U, self.V, 
                   CP.deepcopy(self.E_s, memo), CP.deepcopy(self.E_p, memo), CP.deepcopy(self.E_0, memo), 
                   self.phase)
    
  def read_binary(self, fp):
    self.I, self.Q, self.U, self.V = load_double(fp), load_double(fp), load_double(fp), load_double(fp)

    # warning (note order change): E_s, E_p, E_0 NOT saved in CYCLIC ORDER:
    self.E_p, self.E_s, self.E_0 = loadntuple3Double(fp), loadntuple3Double(fp), loadntuple3Double(fp)

    self.phase = load_double(fp)
  
  def normalized(self):
    # factory method: create a normalized Stoke's vector from the current object.
    rval = None
    if (self.I > finfo(double).eps):
      rval = self.__copy__()
      rval.I, rval.Q, rval.U, rval.V = 1.0, self.Q/self.I, self.U/self.I, self.V/self.I
    else:
      print 'photon_path: StokesV: normalized: *** WARNING: Stoke\'s vector is not normalizable: (I, Q, U, V): (%f, %f, %f, %f) ****'\
        % (self.I, self.Q, self.U, self.V)
        
    return rval

  def rotate_ref_frame_around_e_0(self, phi):
    """
    clone of Pol-MC "StokesV::RotateReferenceFrameAroundPropagationDirection" method.
    """
    e_s, e_p, e_0 = CP.copy(self.E_s), CP.copy(self.E_p), CP.copy(self.E_0)
    cos_phi, sin_phi = cos(phi), sin(phi)    
    
    self.E_s[0] = e_s[0] * cos_phi + e_p[0] * sin_phi
    self.E_s[1] = e_s[1] * cos_phi + e_p[1] * sin_phi
    self.E_s[2] = e_s[2] * cos_phi + e_p[2] * sin_phi
    
    self.E_p[0] = - e_s[0] * sin_phi + e_p[0] * cos_phi
    self.E_p[1] = - e_s[1] * sin_phi + e_p[1] * cos_phi
    self.E_p[2] = - e_s[2] * sin_phi + e_p[2] * cos_phi

    self.rotate_polarization_state(-phi)

  def rotate_polarization_state(self, phi):
    # clone of Pol-MC "StokesV::RotatePolarizationStateBy" method
    cos_2phi = cos(2.0 * phi);
    sin_2phi = sin(2.0 * phi);
    
    Q = self.Q
    self.Q = cos_2phi * self.Q - sin_2phi * self.U
    self.U = sin_2phi * Q  + cos_2phi * self.U
              
  def transformSelfToLabFrame(self, e_s_lab, e_p_lab, e_0_lab):
    """
    modify the reference frame of the Stoke's vector to be consistent with the specified surface e_p, e_s, and normal vectors.
      input parms:
        e_s: reference direction (at surface) for perpendicular polarization
        e_p: reference direction (at surface) for parallel polarization
        e_0: surface normal
        
      warning: input vectors in CYCLIC ORDER (not in same order as Pol-MC method)!
        
      NOTE: If you rotate the surface, you also need to rotate the basis vectors that you give to this method.

      WARNING: this is a clone of Pol-MC: "PhotonCote::MeasureStokesVectorInLabFrame".
      
      KAT: 31.07.2010: "PhotonCote::MeasureStokesVectorInLabFrame": does not respect $ I^2 = Q^2 + U^2 + V^2 $; 
        for this reason, I do not believe that it may not be fully correct (the data's not in yet).
        Stokes parameters are dealt with in two distinct usage scenarios: 
          the first scenario assumes that the parameters are a representation of a coherent optical field vector (e.g. the electric field);
          the second scenario assumes that the parameters are a statistical representation of the optical polarization properties 
          (that is, as might be measured from incoherent fields).
        The first usage requires that the identity: $ I^2 = Q^2 + U^2 + V^2 $\ always be satisfied; 
          whereas the second usage allows the inequality $ I^2 \ge Q^2 + U^2 + V^2 $.
        In polarized Monte-Carlo, it appears that the whole confusion about what to do with Stoke's vectors at an interface relates to
          a disagreement between different scientists about which usage is appropriate for the manipulation of the Stoke's parameters
          during Monte Carlo simulation.  Even more unfortunately, this should not be a very difficult question to answer from the transport equation;
          the correct answer requires only classical electrodynamics, and an understanding of how statistics (or more precisely stochastic integration)
          relates to the application of the Monte Carlo technique to solve the transport equation.
        I personally do not immediately see the answer, and I don't have time to derive it now.  I note that the primary reason that I am 
          uncomfortable with the Pol-MC treatment of the Stoke's parameters is that it is inconsistent; sometimes the first usage scenario is
          applied, and at other times the second.
        For this reason, here I simply clone the Pol-MC implementations of appropriate methods, in order to calculate the fluence at the detector
          from the photon paths.
    """
    if not isinstance(e_s_lab, ndarray):
      e_s_lab = array(e_s_lab)
    if not isinstance(e_p_lab, ndarray):
      e_p_lab = array(e_p_lab)
    if not isinstance(e_0_lab, ndarray):
      e_0_lab = array(e_0_lab)
 
    
    # *** DEBUG *** NAN
    save = self.__copy__()
 
            
    e_s, e_p, e_0 = self.E_s, self.E_p, self.E_0
    e_minus45 = (e_s_lab - e_p_lab)/sqrt(2.0)
    e_plus45 = (e_s_lab + e_p_lab)/sqrt(2.0)
    
    I_s, I_p = intensity_through_linear_polarizer(self, e_s_lab, e_p_lab, e_0_lab)
    I_minus45, I_plus45 = intensity_through_linear_polarizer(self, e_minus45, e_plus45, e_0_lab)
    I_circ_minus45, I_circ_plus45 = intensity_through_circular_polarizer(self, e_s_lab, e_p_lab, e_0_lab)

    # KAT: 31.07.2010: here is where I believe the inconsistency is in this calculation:
    #   the outgoing Stoke's I is retained from the incoming value, but the Q, U, V no longer correspond to this
    #     total intensity value.
    # However, as stated above, this is a _clone_ of Pol-MC: "PhotonCote::MeasureStokesVectorInLabFrame"
    # === possibly "self.I = I_p + I_s" (or something like that?) should be added here. ===
    
    self.Q = (I_p - I_s)
    self.U = (I_plus45 - I_minus45)
    self.V = (I_circ_plus45 - I_circ_minus45)
    
    if (isnan(self.I) or isnan(self.Q) or isnan(self.U) or isnan(self.V)):
      raise RuntimeError, 'NAN Stokes parameters'
    
class photon_event(object):
  """ 
    Information associated with a Monte Carlo "photon" scattering event.
    attributes:
      pos: position tuple (3D)
      dir: direction unit vector      
  """
  __slots__ = ['pos','dir','Stokes','W','t','event_type']
  
  def __init__(self, pos = (0.0, 0.0, 0.0), dir = (0.0, 0.0, 1.0), Stokes = None, W = 0.0, t = 0.0, event_type = 0):
    self.pos = pos
    self.dir = dir
    if (Stokes == None):
      self.Stokes = StokesV()  # create _distinct_ StokesV object
    else:
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
        continue                            # (important: these events are _still_ retained in the path,
                                            #    they just do not participate in the testing)
      any_events = True
      z_min = min(z_min, event.pos[2])
      z_max = max(z_max, event.pos[2]) # here we would project along a normal for the general case
    
    # classify path by its maximum depth:
    if any_events:
      for nr, rtuple in enumerate(depth_ranges):
        if ((z_max > rtuple[0] + OBJECT_SAFETY_DISTANCE) and (z_max < rtuple[1] - OBJECT_SAFETY_DISTANCE)):
          rval[nr].append(photon)
        
  return rval


class photon_depth_classes_iter(object):

  def __init__(self, fp_or_filename, depth_ranges, class_index, event_type = 1):
    self.photon_paths = photon_paths_iter(fp_or_filename)
    self.depth_ranges = depth_ranges
    self.class_index = class_index
    self.event_type = event_type

  def __iter__(self):
    return self
    
  def next(self):
  
    photon_found = False
    photon = None
    
    for photon in self.photon_paths:

      # get depth-range for the path:
      z_min, z_max = finfo(double).max, finfo(double).min
      any_events = False
      for event in photon:
        if (event.event_type != self.event_type):  # ignore photon_events not of specified event_type
          continue                            # (important: these events are _still_ retained in the path,
                                              #    they just do not participate in the testing)
        any_events = True
        z_min = min(z_min, event.pos[2])
        z_max = max(z_max, event.pos[2]) # here we would project along a normal for the general case

      # classify path by its maximum depth:
      if any_events:
        if ((z_max > self.depth_ranges[self.class_index][0] + OBJECT_SAFETY_DISTANCE) 
            and (z_max < self.depth_ranges[self.class_index][1] - OBJECT_SAFETY_DISTANCE)):
          photon_found = True
          break
      
    if not photon_found:
      raise StopIteration
      
    return photon
             
  def reset(self):
    self.photon_paths.reset()


# global constant (maximum number of scattering (or other) events for by-event-number fluence):
MAX_NUMBER_EVENTS = 100

    
def volume_bins(photon_paths, bin_range, shape, 
                bin_type = 'fluence', smooth = False, mu_a = 1.0, mu_s = 0.0, 
                event_type = 1, 
                detector_perp = (1.0, 0.0, 0.0), detector_para = (0.0, 1.0, 0.0), detector_normal = (0.0, 0.0, 1.0), 
                minAcceptanceCosine = 0.0, maxAcceptanceCosine = 1.0,
                count_event_number = False):
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
   
    event_type: type of photon_event to add to fluence or energy accumulation (1:scattering_event)
    
    detector_perp, detector_para, detector_normal: detector surface basis vectors
    
    minAcceptanceCosine, maxAcceptanceCosine: where antiparallel to detector normal is 1.0
    
    count_event_number: tabulate 2D fluence by number of events of specified type in the photon path
    
      NOTES: 
        -- energy deposition is scaled by mu_a / mu_t, so default values return non-rescaled energy (i.e. which is the same as fluence);
        -- fluence deposition is as unscaled photon weight, which from sampling considerations, 
             should give fluence in units of "mu_t" (times whatever the units of photon weight are); 
             that is , fluence should be rescaled by "1.0/mu_t" to get physical units.
        .
  returns:
    <3D array of volume bins>, <total fluence at the detector in "polmc_data" "Stokes_data" dict format>
    
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
  
  # order the acceptance-cosine tuple:
  acceptanceCosine = (min(minAcceptanceCosine, maxAcceptanceCosine), max(minAcceptanceCosine, maxAcceptanceCosine))
  
  # extinction coefficient:
  mu_t = mu_a + mu_s 
  
  # ----------- return values: ----------------------------------------------  
  bins = zeros(shape,dtype=float) 
  
   # detector_fluence is in "polmc_data" "Stokes_data" format (but containing scalar entries instead of arrays):
  if not count_event_number:
    detector_fluence = {'I':0.0, 'Q':0.0, 'U':0.0, 'V':0.0, 
                        'minAcceptanceCosine':minAcceptanceCosine, 
                        'maxAcceptanceCosine':maxAcceptanceCosine} 
  else:
    detector_fluence = 
      [ {'I':0.0, 'Q':0.0, 'U':0.0, 'V':0.0, 
         'minAcceptanceCosine':minAcceptanceCosine, 
         'maxAcceptanceCosine':maxAcceptanceCosine} for N_events in range(MAX_NUMBER_EVENTS) ]   
  # -------------------------------------------------------------------------                    
  
  for photon in photon_paths:
    # ignore empty photon_paths
    if not photon:
      continue

    transmit_event = photon[-1]
    if (transmit_event.event_type != 3):
      print 'photon_path: volume_bins: WARNING: last event:%d in the photon path is not \"kInterfaceTransmitEvent\":3'\
        % transmit_event.event_type

    # ignore photons outside of the acceptance-angle:
    #   cosign at the detector surface (s.t. antiparallel is 1.0):
    cosine = -normalized_dot_product(photon[-1].dir, detector_normal)
    if ((acceptanceCosine[0] > cosine) or (cosine > acceptanceCosine[1])):
      continue
     
    N_events = 0            
    for event in photon:
      if (event.event_type != event_type):  # ignore photon events not of the specified type
        continue
        
      N_events += 1  
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
      
    # accumulate the detector fluence: (important: see notes at "StokesV: transformToLabFrame" above):
    localStokes = transmit_event.Stokes.normalized()
    localStokes.transformSelfToLabFrame(detector_perp, detector_para, detector_normal)
    
    if not count_event_number:    
      detector_fluence['I'] += localStokes.I * transmit_event.W 
      detector_fluence['Q'] += localStokes.Q * transmit_event.W
      detector_fluence['U'] += localStokes.U * transmit_event.W
      detector_fluence['V'] += localStokes.V * transmit_event.W  
    else:
      if N_events >= MAX_NUMBER_EVENTS:
        raise RuntimeError: 'photon_path: volume_bins: too many events of specified type in the path'
      detector_fluence[N_events]['I'] += localStokes.I * transmit_event.W 
      detector_fluence[N_events]['Q'] += localStokes.Q * transmit_event.W
      detector_fluence[N_events]['U'] += localStokes.U * transmit_event.W
      detector_fluence[N_events]['V'] += localStokes.V * transmit_event.W  
   
  return bins, detector_fluence


def detector_fluence(photon_paths, e_s, e_p, e_0, minAcceptanceCosine, maxAcceptanceCosine, 
                     N_cosine_bins,
                     event_type = 1,
                     count_event_number = False):
  """
    Total fluence at the detector from photon paths: same format as in polmc_data (w.r.t. accumulated cosine-bins, except only one spatial bin).
      Last photon_events in the path shall correspond to the interface transition.

    input parms:
      photon_paths: iterable
      
      e_s, e_p, e_0: surface basis vectors (i.e. "er", "el", and "normal", using Pol-MC's notation)
        NOTE: if you rotate the detector surface, these basis vectors should be rotated as well (using the same Euler rotation matrix);
          WARNING: Pol-MC uses the X, Y', Z'' convention (which is very rare and strange for Euler rotations);
            this method uses the Z, X', Z'' convention.  
            For example Pol-MC (0.0, pi/4.0, 0.0) is (-pi/2.0, pi/4.0, 0.0) using the Z, X', Z'' convention.
             
      minAcceptanceCosine, maxAcceptanceCosine:       
        range of photon-event direction-cosines (with respect to surface normal) corresponding to accumulated events;
          note: sign of dot-product is reversed (i.e. anti-parallel photon_event direction and surface normal correspond to cosine = 1.0).
          
      N_cosine_bins: number of accumulation bins.
      
      event_type: type of event for which to tabulate 2D fluence (1:scattering_event)
      
      count_event_number: tabulate 2D fluence by number of events of specified type
    
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
  
  if not count_event_number:
    bins = [polmc_data.Stokes_data.zeros(cosine_bin_edges[n-1], cosine_bin_edges[n], shape = None)  for n in range(1, N_cosine_bins)]
  else:
    bins = [ [polmc_data.Stokes_data.zeros(cosine_bin_edges[n-1], cosine_bin_edges[n], shape = None)  for n in range(1, N_cosine_bins)]
             for N_events in range(MAX_NUMBER_EVENTS) ]
  
  def cosine_bin_index(cosine, bin_range, N_bins):
    n = int(floor((cosine - bin_range[0])/(bin_range[1] - bin_range[0]) * float(N_bins)))
    if (n >= N_cosine_bins):
      n = None
    return n
    
  for photon in photon_paths:
    transmit_event = photon[-1]
    if (transmit_event.event_type != 3):  # kInterfaceTransmitEvent: note: photon statistics are recorded _prior_ to killing the photon.
      raise RuntimeError, 'photon_path: detector_fluence: last event in photon_path is not \"kInterfaceTransmitEvent\"'
    
    cosine = -normalized_dot_product(event.dir, e_0)  # "cosine < 0" => exit-events, "cosine > 0" => entry-events

    # only count number of events if requested:
    N_events = 0
    if count_event_number:
      for event in photon:
        if event.event_type == event_type:
          N_events += 1

    localStokes = transmit_event.Stokes.normalized()
    localStokes.transformSelfToLabFrame(e_s, e_p, e_0)
    
    n = cosine_bin_index(cosine, (minAcceptanceCosine, maxAcceptanceCosine), N_cosine_bins)
    if (n != None):
      # note: index will be None if cosine is out-of-range,
      #   which is not an error here (bins are accumulated only for the region of interest specified via (minAcceptanceCosine, maxAcceptanceCosine)).
      if not count_event_number:
        bins[n]['I'] += localStokes.I * transmit_event.W 
        bins[n]['Q'] += localStokes.Q * transmit_event.W
        bins[n]['U'] += localStokes.U * transmit_event.W
        bins[n]['V'] += localStokes.V * transmit_event.W
      else:
        bins[N_events][n]['I'] += localStokes.I * transmit_event.W 
        bins[N_events][n]['Q'] += localStokes.Q * transmit_event.W
        bins[N_events][n]['U'] += localStokes.U * transmit_event.W
        bins[N_events][n]['V'] += localStokes.V * transmit_event.W
        
  
  # accumulate cosine bins to Pol-MC "expanding angle" format:
  for n in range(N_cosine_bins-1,0,-1):
    if not count_event_number:
      bins[n-1]['I'] += bins[n]['I']
      bins[n-1]['Q'] += bins[n]['Q']
      bins[n-1]['U'] += bins[n]['U']
      bins[n-1]['V'] += bins[n]['V']
      bins[n-1]['maxAcceptanceCosine'] = 1.0 
    else:
      for N in range(MAX_NUMBER_EVENTS):
        bins[N][n-1]['I'] += bins[n]['I']
        bins[N][n-1]['Q'] += bins[n]['Q']
        bins[N][n-1]['U'] += bins[n]['U']
        bins[N][n-1]['V'] += bins[n]['V']
        bins[N][n-1]['maxAcceptanceCosine'] = 1.0        
     
  return bins
  
      
