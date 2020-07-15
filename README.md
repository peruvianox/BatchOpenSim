# BatchOpenSim

## Overview

**What:** Automate processing of walking biomechanical data (markers and ground reaction forces) through OpenSim musculosleletal modelling using the Matlab to OpenSim API

**Why:** automated processing of multiple trials speeds up processing time and reduces the likelihood of errors

**Who:** human biomechanics researchers, especially those studying walking while on a dual belt treadmill with embedded force platforms

**How:** as described in ___


### Setup Requirements
- downladed Opensim (https://simtk.org/projects/opensim) and Matlab (https://www.mathworks.com/products/matlab.html)
- setup OpenSim API, instructions at: https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

### Input Requrements
- matching .trc and .anc files for each trial recorded (filenames matched by spelling)
- OpenSim model and set up files for each task (see OpenSimProcessingFiles folder) 

 
## Instructions

### Setup and Code

First you'll need to download OpenSim 4.0 and then follow the instructions below to set up the Matlab API. https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab
 
Processing Codes can be found at:  https://github.com/peruvianox/BatchOpenSim
 
#### Data Formatting
You’ll need to arrange your participant ANC and TRC files within a single folder. Then place all your ANC & TRC files for each trial within the first level of each participant folder. The batch script will create an OpenSim folder and FORCES files.  
 
#### Calibration and Coordinate System 
Within the ABL_OpenSim_Setup_Batch matlab script, you'll need your own force plate calibration files (lines 65-74), and set your own coordinate system conversions (lines 245-256).
 
#### Selecting Time Window to Process
If you are going to run computed muscle control (CMC) or static optimization (SO), it is helpful to select a small window (<2 seconds) for which to run the simulation. Longer windows may lead to inaccurate results and require high computing power. 


### OpenSim Data processing
#### Setup Files
In the OpenSimProcessingFiles folder, you'll find the generic setup files for each state of opensim processing (scale, inverse kinematics, residual reduction, and computed muscle control), which are copied into each subject's folder and rewritten with their own parameters and filenames. When you have your own data processed, you can load and run these manually to see how they work. 
#### Scaling
Scale each model to match each subjects size, using marker locations to scale each body segment individually. 

#### Inverse Kinematics
Fit the scaled model to recorded marker trajectories for each trial. 

#### Residual Reduction Algorithm
Align ground reaction forces with body dynamics during recorded motions in an iterative fashion. 

#### Computed Muscle Control
Calculate the muscle-tendon unit forces required to generate the recorded movement. 

#### Muscle Analysis


#### Static Optimization
Calculate the 

#### Inverse Dynamics


## Overview of Scripts
- ABL_Batch_OpenSim.m - parent code for batch processing multiple subjects and trails in OpenSim (including scaling, inverse kinematics, residual reduction algorithm, computed muscle control, static optimization, muscle analysis, and inverse dynamics)
- ABL_OpenSim_Setup_Batch.m - convert TRC and ANC files from original data formats and coordinate systems to formats usable by OpenSim. Options to identify crossover steps, parse windows from full trial, and add torso markers to lower-body only walking data. 
- ABL_Scale.m - Scale each subject using anatomical marker definitions to specify segment dimensions along their primary axes. See https://www.youtube.com/watch?v=ZG7wzvQC6eU&t=3083s at 13:45 for specifics

## OpenSim Resources
OpenSim Documentation Site
https://simtk-confluence.stanford.edu:8443/display/OpenSim/Documentation
  
Informational Videos
- OpenSim Webinar: Tips and Tricks for Data Collection, Scaling and Inverse Kinematics in OpenSim https://www.youtube.com/watch?v=ZG7wzvQC6eU&t=3083s 
- The Scale Tool: Evaluating the Results https://www.youtube.com/watch?v=7STRMJefzpI

## License


