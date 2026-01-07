This is a software package developed in ASEN 5050 Spaceflight Dynamics. It is designed to house several 
functions and objects to make orbit analysis easier. It uses 5 main objects and has a tutorial file, as 
well as a few utility functions. 

  -  ephemeris.py: This file houses the Ephemeris object, which is responsible for storing satellite 
  position, time, orbital elements, frame and epoch. This is built to work with the functions. 
  -  graphing_utils.py: This file houses the Graph object, which is designed to use Plotly to graph 
  whatever is needed from Ephemeris objects. 
  -  math_tools.py: This file has a few functions I wrote including datetime interpolaters, cubic spline 
  interpolators and sin and cosine functions to import in complex functions for efficiency. 
  -  numerical.py: This file contains two objects- the Integration object and the ForcingFunction object. 
  These are used for numerical propagation and the ForcingFunction object allows for adding different 
  forcing functions together and studying how perturbations affect satellite orbit. 
  -  orbit_package: This contains the Orbit object, which is used for all of the analysis and analytical 
  math on orbits. It houses functions to create ephemeris from orbit elements, to go between reference
  frames, to solve Lambert's Problem and many other things. 
  -  threebp_objects.py: This file houses the Cr3bp object which is used in the circular restricted 3 body 
  problem. 
  -  tutorial.ipynb: This is a Jupyter Notebook for demonstrating the capabilities of the objects. This is 
  still in development. 