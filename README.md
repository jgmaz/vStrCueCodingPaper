# vStrCueCodingPaper

This repository hosts the code and preprocessed data files used for Gmaz,
Carmichael & van der Meer, "Persistent coding of outcome-predictive cue features in the rat nucleus accumbens" (2018) ([preprint](https://www.biorxiv.org/content/early/2018/08/27/300251)).

The code makes use of the [vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab), and was tested on Windows 7 using MATLAB R2015a.

Please see [here](https://github.com/jgmaz/vStrCueCodingPaper/blob/master/LICENSE.md) for licensing information.

To reproduce the results in the paper, start by cloning this repository. Then, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared')); % lab codebase
addpath(genpath('\GitHub\vandermeerlab\code-matlab\toolboxes\MClust-3.5')); % used for utility functions
addpath(genpath('\GitHub\vStrCueCodingPaper\Code')); % paper repo
```
The workflow for the analysis and figure generation is found in **MASTER_workflow.m**. This script is divided into cells that run the various analyses reported in the paper, saving the output at each step to be used subsequently for higher-level analysis and figure generation. 

The top of this master workflow file contains two constants whose values you will need to set. They are:
  * `GITHUB_PATH`: location of the top-level folder containing both this repo and the vandermeerlab codebase (example: `C:\Users\mvdm\Documents\GitHub\`)
  * `DATA_ROOT`: location of the top-level folder containing the data (example: `C:\data\vStrCueCoding\`, which should have subfolders `R053`, `R060`, etc.). Note that the analysis will create some folders in this location to store intermediate results files.

After that, you should be all set to run the analysis. Below is a brief summary of the function of each cell.

### Behavior summary performance (needed for Figure 2)

This cell takes the metadata for each recording session (see [here](https://github.com/jgmaz/vStrCueCodingPaper/blob/eeb72f6ecce44c0bcd4c6d7a3a3cf97342ffbbe1/METADATA.md) for description) and generates the summary performance variable used for **Figure 2C**.

### Plot behavior example and summary

This cell takes the behavioral performance summary variable generated above, and a previously generated learning curve example for R060, and generates **Figure 2B,C**.

### Preprocess units (needed for Figures 3,4,5,6,7) 

This cell performs the preprocessing step of the spike .t files, extracting trial-by-trial firing rates, presence or absence of cue-modulation, PETHs, and rasters for each unit. 

The output is a .mat file that is saved in a common *spike data* folder for all units, defined at the top of the main workflow script. 

### Cue-onset GLM (needed for Figure 4)

This cell performs the sliding window GLM for time bins surrounding cue-onset.

The **genGLM** and **genGLMshuff** functions take the processed spike data, and output the result of the GLMs in the *analysis files* folder, where the **genGLMstats** function generates the summary variable used for visualization. 

### Plot cue-onset GLM

This cell plots the output variable of the previous cell (*GLM_window*) to generate **Figure 4A** and **Figure 4 supplement 1**.

### Cue-onset LDA (needed for Figure 4)

This cell performs the classifier used in **Figure 4B**.

The **genLDA** and **genLDAshuff** functions take the processed spike data, and output classifier performance in *analysis files*, where *genLDAstats* creates a summary variable used for plotting.

### Plot cue-onset LDA

This cell plots the output variable of the previous cell (*Class_accuracy*) to generate **Figure 4B**.

### Nosepoke & outcome GLM (needed for Figure 7)

This cell performs a similar function as the **_Cue-onset GLM_** cell, performing the sliding window GLM for time bins surrounding nosepoke and outcome receipt.

### Plot nosepoke & outcome GLM

This cell plots the output variable of the previous cell (*GLM_window*) to generate **Figure 7A,B** and **Figure 7 supplement 1A-D**.

### Recoded coefficients (needed for Figure 4,7)

This cell performs the analyses that generate the correlation matrices for recoded coefficient arrays.

The **genRecode** function first gathers information about which units show modulation by the cue, accessing the spiking data from *spike data*, then generates recoded coefficient arrays and correlation matrices using the output of the various sliding window GLMs in *analysis files*, and outputs a summary variable (*GLM_coeff*) in *analysis files*.

**genRecodeShuff** takes the recoded arrays from **genRecode** and shuffles their across-unit ordering to compare *GLM_coeff* to chance levels.

### Plot recoded correlation matrices

This cell plots the output variable of the previous cell (*GLM_coeff*) to generate **Figure 4D**, **Figure 7D-F**, and **Figure 7 supplement 1E-F**.

### Plot cue-onset examples

This cell takes the processed spiking data for the examples highlighted in the text, and visualizes the PETHs and rasterplots seen in **Figure 3** and **Figure 3 supplements 1,2**.
d
### Plot nosepoke examples

This cell takes the processed spiking data for the examples highlighted in the text, and visualizes the PETHs and rasterplots seen in **Figure 6** and **Figure 6 supplements 1,2**.

### Cue-onset distributed firing (needed for Figure 5)

This cell uses the PETHs generated from genPROCESS in a cell above, and orders them according to their maximum and minumum firing for identity, location, and outcome at time of cue-onset. The output is used for **Figure 5** and **Figure 5 supplement 1**.

### Plot cue-onset distributed firing

This cell uses the variables created in the previous cell to generate **Figure 5** and **Figure 5 supplement 1**.

### Nosepoke distributed firing (needed for Figure 7 supplement 2)

This cell uses the PETHs generated from genPROCESS in a cell above, and orders them according to their maximum and minimum firing for identity, location, and outcome at time of nosepoke. The output is used for **Figure 7 supplement 2**.

### Plot nosepoke distributed firing

This cell uses the variables created in the previous cell to generate **Figure 7 supplement 2**.

### Plot scatterplot

This cell uses the output of the sliding window GLM centered on cue-onset to generate **Figure 4 supplement 2**, a scatterplot showing the firing rate for the light and sound block for each of the cue-modulated units.

### Attributions

This codebase makes use of [errorbarxy](https://www.mathworks.com/matlabcentral/fileexchange/4065-errorbarxy), [heatmap](https://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps), [shadedErrorBar](https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar), and [subtightplot](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot).
