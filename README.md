# ISMaRTS24
A collection of models and interactive MATLAB scripts for variational calculus problems related to the classical brachistochrone problem. Some of these files (.m) may be compatible with GNU Octave, but this is not guaranteed. Some are not because they use proprietary MATLAB toolboxes. I am going to attempt to make it compatible so that these results and tools are available to anyone.  

These are not complete models and they need some work to be polished up. Further, I have additional ideas for expanding these models into more complex scenarios (e.g., a brachistochrone under drag and variable gravity in spherical coordinates where the user specifies the endpoints). More importantly, I need to create comprehensive documentation explaining the code for others (and even myself at times!). 

----- **INDIVIDUAL FILE DESCRIPTIONS** -----

--- Brachistochrone_Plotter.m ---
   
   This will plot N brachistochrones, lines, and parabolas (specified with the b input variable). These are spaced at an interval of 1 so that every whole number between 1 and b (inclusive) will have a set of three curves terminating there. The y value is fixed at -1, so we are running a loop through (n,-1) from 1 to b. For this reason, it is best to keep the b value less than, say, 25.
   Further, a step size needs to be specified for the x axis partition and the partition of the parameter theta (which is used to plot the brachistochrone according to a set of parametric equations). I wouldn't go below 0.001 as you are not really getting any further accuracy of importance and are only slowing down the computations. I usually run it between 0.1 and 0.005.
   The purpose of this script is to make a nice visualization for comparing how curves less optimal than the brachistochrone match up with this "path of least time". Other scripts dive deeper into this, comparing descent times and arc lengths of non-optimal curves to that of the brachistochrone.

--- bfun.m ---
