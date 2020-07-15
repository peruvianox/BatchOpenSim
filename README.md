# BatchOpenSim

OVERVIEW
What: Automate processing of walking biomechanical data (markers and ground reaction forces) through OpenSim musculosleletal modelling using the Matlab to OpenSim API

Why: automated processing of multiple trials speeds up processing time and reduces the likelihood of errors

Who: human biomechanics researchers, especially those studying walking while on a dual belt treadmill with embedded force platforms

How: as described in ___


SETUP CONFIGURATION
- downladed Opensim (https://simtk.org/projects/opensim) and Matlab (https://www.mathworks.com/products/matlab.html)
- setup OpenSim API, instructions at: https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

INPUT REQUIREMENTS
- matching .trc and .anc files for each trial recorded (filenames matched by spelling)
- OpenSim model and set up files for each task (see OpenSimProcessingFiles folder) 


 
INSTRUCTIONS
Setup and Code
First you'll need to download OpenSim 4.0 and then follow the instructions below to set up the Matlab API. https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab
 
Processing Codes can be found at:  https://github.com/peruvianox/BatchOpenSim
 
Data Formatting
You’ll need to arrange your participant ANC and TRC files within a single folder: 
 
and then have all your files (ANC & TRC)  for each trial within the first level of each participant folder: 

 
The batch script will create the OpenSim folder and FORCES files shown above. 
 
 
 
Calibration and Coordinate System 
Within the ABL_OpenSim_Setup_Batch matlab script, you'll need your own force plate calibration files (lines 65-74), and set your own coordinate system conversions (lines 245-256).
 
 
Selecting Time Window to Process



OpenSim Data processing
Setup Files
In the OpenSimProcessingFiles folder, you'll find the generic setup files for each state of opensim processing (scale, inverse kinematics, residual reduction, and computed muscle control), which are copied into each subject's folder and rewritten with their own parameters and filenames. When you have your own data processed, you can load and run these manually to see how they work. 
Scaling


Inverse Kinematics

Residual Reduction Algorithm

Computed Muscle Control


Muscle Analysis

Static Optimization

Inverse Dynamics


OVERVIEW OF SCRIPTS
- ABL_Batch_OpenSim.m - parent code for batch processing multiple subjects and trails in OpenSim (including scaling, inverse kinematics, residual reduction algorithm, computed muscle control, static optimization, muscle analysis, and inverse dynamics)
- ABL_Setup_Batch.m

OpenSim RESOURCES
OpenSim Documentation Site
https://simtk-confluence.stanford.edu:8443/display/OpenSim/Documentation
 
INFORMATIONAL VIDEOS
https://www.youtube.com/watch?v=7STRMJefzpI
https://www.youtube.com/watch?v=ZG7wzvQC6eU&t=3083s 

LICENSE

