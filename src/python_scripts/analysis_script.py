import photon_path as PP
from numpy import *

test_filename = 'photons.dat'
photon_paths = PP.load_photon_paths(test_filename)  



N_events = [0, 0, 0, 0, 0, 0]

for photon in photon_paths:
  for event in photon:
    N_events[event.event_type] += 1

print N_events

# fluence for different layers



#### note for extracting polmc_data pickle file
world._data.key()
world['fluenceBinNumber'].shape
world['objects'].keys()
world['objects']['angle_detector']['interfaces']['forward']['StokesV'][0]['I'].sum()  #--->this gives us intensity at detector

