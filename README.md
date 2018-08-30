# vStrCueCodingPaper


Code used for Gmaz,
Carmichael & van der Meer, "Persistent coding of outcome-predictive cue features in the rat nucleus accumbens" (2018) ([preprint](https://www.biorxiv.org/content/early/2018/08/27/300251)).

This repo makes use of the [vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);

Once you have downloaded the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
cd('\promoted\analysis files\'); % replace paths with yours
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared')); % lab codebase
addpath(genpath('\GitHub\vandermeerlab\code-matlab\toolboxes\MClust-3.5')); % MClust-3.5
addpath(genpath('\GitHub\vStrCueCodingPaper\Code')); % paper repo
```
The workflow for the analysis and figure generation is found in **MASTER_workflow**. It is divided into cells that run the various analyses reported in the paper, saving the output at each step to be used subsequently for higher-level analysis and figure generation. 

Below is a brief summary of the function of each cell.

### Behavior example learning curves (NOT COMPLETED)

This cell takes the metadata from each recording session for R060 and generates the learning curves seen in **Figure 2B**.

### Behavior summary performance (NOT COMPLETED)

This cell takes the metadata for each recording session and generates the summary performance figure and statistics seen in **Figure 2C**.

### Preprocess units (needed for Figures 3,4,5,6,7)

This cell performs the preprocessing step of the .t files, loading the spikes, trial-by-trial firing rates, presence or absence of cue-modulation, PETHs, and rasters for each unit. 

The output is a .mat file that is saved in a common *destination* for all units. 

The *directory* is the location of the data (e.g. if data for a session is in ‘E:\vStr-cuecoding\promoted\R057\R057-2015-02-25’, the directory would be ‘E:\vStr-cuecoding\promoted\’).

### Cue-onset GLM (needed for Figure 4)

This cell performs the sliding window GLM for time bins surrounding cue-onset.

The **genGLM** and **genGLMshuff** functions take the processed spike data in *directory*, and output the result of the GLMs in the *destination*, where the **genGLMstats** function generates the summary variable used for visualization. 

### Plot cue-onset GLM

This cell plots the output variable of the previous cell (*GLM_window*) to generate **Figure 4A** and **Figure 4 supplement 1**.

### Cue-onset LDA (needed for Figure 4)

This cell performs the classifier used in **Figure 4B**.

The **genLDA** and **genLDAshuff** functions take the processed spike data in *directory*, and output classifier performance in *destination*, where *genLDAstats* creates a summary variable used for plotting.

### Plot cue-onset LDA

This cell plots the output variable of the previous cell (*Class_accuracy*) to generate **Figure 4B**.

### Nosepoke & outcome GLM (needed for Figure 7)

This cell performs a similar function as the **_Cue-onset GLM_** cell, performing the sliding window GLM for time bins surrounding nosepoke and outcome receipt.

### Plot nosepoke & outcome GLM

This cell plots the output variable of the previous cell (*GLM_window*) to generate **Figure 7A,B** and **Figure 7 supplement 1A-D**.

### Recoded coefficients (needed for Figure 4,7)

This cell performs the analyses that generate the correlation matrices for recoded coefficient arrays.

The **genRecode** function first gathers information about which units show modulation by the cue, accessing the spiking data from *spike_directory*, then generates recoded coefficient arrays and correlation matrices using the output of the various sliding window GLMs in *directory*, and outputs a summary variable (*GLM_coeff*) in *destination*.

**genRecodeShuff** takes the recoded arrays from **genRecode** and shuffles their across-unit ordering to compare *GLM_coeff* to chance levels.

### Plot recoded correlation matrices

This cell plots the output variable of the previous cell (*GLM_coeff*) to generate **Figure 4D**, **Figure 7D-F**, and **Figure 7 supplement 1E-F**.

### Plot cue-onset examples

This cell takes the processed spiking data for the examples highlighted in the text, and visualizes the PETHs and rasterplots seen in **Figure 3** and **Figure 3 supplements 1,2**.
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % where spiking data resides

### Plot nosepoke examples

This cell takes the processed spiking data for the examples highlighted in the text, and visualizes the PETHs and rasterplots seen in **Figure 6** and **Figure 6 supplements 1,2**.
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % where spiking data resides

### Cue-onset distributed firing (needed for Figure 5)

In-progress. Generates the data for **Figure 5**.

### Nosepoke distributed firing (needed for Figure 7 supplement 2)

### Plot nosepoke distributed firing

In-progress. Generates the data for **Figure 7 supplement 2**.

### Plot scatterplot

In-progress. Generates **Figure 4 supplement 2**.

### Attributions

This codebase makes use of (heatmap)[https://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps], (shadedErrorBar)[https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar], and (subtightplot)[https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot].