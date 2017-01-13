# PolMC
Polarized Monte Carlo simulation used in my doctoral dissertation
- C++ code simulating photonics traveling in biological tissues in custom geometry, with geometries such as objects, layers, interfaces, etc in the user defined space.
- C++ code was parallelizable and ran with MPI on clusters at TACC (Texas Advanced Computing Center)
- All physical objects/geometries are defined with their own photonic properties saved in .plist (XML) files
- Results are saved in large binary files and analyzed (aggregated) in python

Some published results at
- http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=719199
- https://repositories.lib.utexas.edu/bitstream/handle/2152/ETD-UT-2010-12-2047/KAN-DISSERTATION.pdf?sequence=1
