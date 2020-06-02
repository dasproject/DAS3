# DAS3: 3D shoulder and elbow model

Developer team: Ton van den Bogert (a.vandenbogert@csuohio.edu), Dimitra Blana (dimitra.blana@abdn.ac.uk), Ed Chadwick (edward.chadwick@abdn.ac.uk)

DAS3 is the third generation of the Dynamic Arm Simulator, a musculoskeletal model of the shoulder and arm. The arm simulator is developed as part of the NIH contract “Brain-Controlled Hybrid Functional Electrical Stimulation”, led by Dr. Robert F. Kirsch at Case Western Reserve University.

The model was originally developed for real-time simulation.  When used together with a suitable visualisation environment, the DAS3 model can provide a virtual reality tool, with a human in the control loop, for testing brain interfaces for functional electrical stimulation (FES) [1].

However, the model comes with an API that allows other types of applications.  Offline simulations (faster or slower than real time) can be performed to simulate experiments or do optimization of control systems.  Trajectory optimizations can be done using collocation methods to do state estimation, predictive simulations, or design optimization [2].

This repository has work in progress to port the project from Opensim 3.3 to Opensim 4.0.  The Opensim 3.3 version is still at https://simtk.org/projects/das.  Binaries and documentation are available, and source code is accessible for team members.

When the port is completed, the most recent version can be downloaded here.  We hope that users will contribute to the project.

To use DAS3, the following environment is needed:
- Matlab 2018b or later
- Opensim 4.0 [3] 
- Opensim Matlab API [4].  Run the Opensim GUI once before installing the Matlab API to ensure that the user folders are created properly.

### Contents of documentation

- File repository
- How to install and run the model
- Model reference
- Model testing

### References
[1] Chadwick EK, Blana D, Kirsch RF, van den Bogert AJ (2014) Real-time simulation of three-dimensional shoulder girdle and arm dynamics. *IEEE Trans Biomed Eng* **61**(7): 1947-1956.

[2] van den Bogert AJ, Blana D, Heinrich D (2011) Implicit methods for efficient musculoskeletal simulation and optimal control. *Procedia IUTAM* **2**: 297-316. 

[3] https://opensim.stanford.edu/

[4] https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab 



