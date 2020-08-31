# BatchOpenSim (BOS)

## Overview

**What:** automate processing of walking biomechanical data (markers and ground reaction forces) through OpenSim musculosleletal modelling

**Why:** automated processing saves time, reduces errors, and improves consistency

**Who:** human biomechanics researchers, especially those studying walking while on a dual belt treadmill with embedded force platforms

**How:** as described here: https://www.biorxiv.org/content/10.1101/2020.07.31.230698v1

## Instructions

### Setup
- download Opensim (https://simtk.org/projects/opensim) and Matlab (https://www.mathworks.com/products/matlab.html)
- setup OpenSim API for Matlab, instructions at: https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab
- Download  BOS code at:  https://github.com/peruvianox/BatchOpenSim

### Inputs
- matching .trc and .anc files for each trial recorded (filenames matched by spelling)
- OpenSim model and set up files for each task (see OpenSimProcessingFiles folder) 
 
### Format Input Data
First, build a BOS database by creating a study folder containing folders for each subject, like this: 
![Study Folder](https://github.com/peruvianox/BatchOpenSim/blob/master/doc/Study_Folder.PNG)

Then place all your ANC & TRC files for each trial within the first level of each participant folder. 
![Input Files](https://github.com/peruvianox/BatchOpenSim/blob/master/doc/Input_Files.PNG)

BatchOpenSim will create the *OpenSim folder*, convert all original files to .mot format, and populate the new OpenSim folder with that subject's data. All future modelling results will be saved to each subject's *OpenSim folder*. 

### Calibration and Coordinate System 
Within the ABL_OpenSim_Setup_Batch matlab script, you'll need your own force plate calibration files (lines 65-74), and set your own coordinate system conversions (lines 245-256).

### Data Structure
BatchOpenSim will save a heirarchical structure containing data for each subject (called *Subjects*), containing its file path, directory information, demographics, and trials completed. 
![Subject View](https://github.com/peruvianox/BatchOpenSim/blob/master/doc/Subject_view.png)

Going into each subject's "Trials" data, you can then see its folder, type, associated files, ground reaction forces, temporal spatial data, and other info)
![Trial View](https://github.com/peruvianox/BatchOpenSim/blob/master/doc/Trial_view.png)
 
### Time Window
If you are going to calculate muscle-tendon unit control (computed muscle control or static optimization), it is helpful to select a small window of time (<2 seconds) over which to run the simulation. Long windows of processing require substantial computing power and may results may be questionable. CMC and SO require a 0.03 s buffer to determine inital parameters, so we add 0.05 s to the start and end of the selected window. Processing window times are added in the *ABL_OpenSim_SetupBatch.m* file, search for "BestTimes" variable. 

### OpenSim Data processing
#### Setup Files
In the OpenSimProcessingFiles folder, you'll find the generic setup files which are copied into each subject's folder and rewritten with their own parameters and filenames. After you have run BatchOpenSim, you can manually load and run the setup files to check results. 

#### Scaling
Scale each model to match each subjects size, using marker locations to scale each body segment individually. 

#### Inverse Kinematics
Fit the scaled model to recorded marker trajectories for each trial. 

#### Residual Reduction Algorithm
Align ground reaction forces with body dynamics during recorded motions in an iterative fashion. 

#### Computed Muscle Control
Calculate the muscle-tendon unit forces required to generate the recorded movement. 

#### Muscle Analysis
more info coming soon

#### Static Optimization
more info coming soon

#### Inverse Dynamics
more info coming soon

## Major BOS Functions
- ABL_Batch_OpenSim.m - parent code for batch processing multiple subjects and trails in OpenSim (including scaling, inverse kinematics, residual reduction algorithm, computed muscle control, static optimization, muscle analysis, and inverse dynamics)
- ABL_OpenSim_Setup_Batch.m - convert TRC and ANC files from original data formats and coordinate systems to formats usable by OpenSim. Options to identify crossover steps, parse windows from full trial, and add torso markers to lower-body only walking data. 
- ABL_Scale.m - Scale each subject using anatomical marker definitions to specify segment dimensions along their primary axes. See https://www.youtube.com/watch?v=ZG7wzvQC6eU&t=3083s at 13:45 for specifics

## OpenSim Resources
OpenSim Documentation Site
https://simtk-confluence.stanford.edu:8443/display/OpenSim/Documentation
  
Informational Videos
- OpenSim Webinar: Tips and Tricks for Data Collection, Scaling and Inverse Kinematics in OpenSim https://www.youtube.com/watch?v=ZG7wzvQC6eU&t=3083s 
- The Scale Tool: Evaluating the Results https://www.youtube.com/watch?v=7STRMJefzpI

## License & Redistribution
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

