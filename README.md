# vStrCueCodingPaper


Figure 4:

To get Figure 4 from the paper, first run MASTER_behavior_example to extract the variables needed for the example plots in A and B. Then run MASTER_behavior to extract the summary performance measures used in C and D. Then run MASTER_behavior_figure to generate Figure 4.

**Note: Access to behavioral data in MASTER_behavior_example is currently hardcoded in the script. Can save all the .mat files it uses and put them in a folder for the paper.**

**Note: In the current form of MASTER_behavior, each behavioral session was manually run in the first half of the script, before generated the summary values in the second half of the script. Will include a call to the .mat files containing each session in the final version.**


Recording pre-analysis:

All the neural analyses use variables generated from MASTER_main_workflow, whose purpose is to take the .t files for each cell, extract the spikes, sync with the behavioral data, and generate the firing measures used for subsequent analyses.

**Note: There are multiple functions that sync the behavioral data with spiking data as the organization of inputs/outputs to the TTL box changed when we switched running rooms.** 


Figure 5:

To get Figure 5 from the paper, run MASTER_GLM_example_cells.


Figure 6:

To get Figure 6 from the paper, run MASTER_GLM.

**Note: Running MASTER_GLM in its current form will also generate several figures that are not used in Figure 6 (leftover from when developing the pipeline).**

**Note: There is an option specifying which rat to use, as my structuring of variables slightly changed midway through analysis (e.g. storing each block of trials in a separate cell within the same struct in metadata versus two separate structs in MASTER_main_workflow. Will rerun MASTER_main_workflow to standardize this and remove this option**


Figure 7:

To get Figure 7 from the paper, run MASTER_population_average.


Figure 8:

To get Figure 8 from the paper, run MASTER_task_tiling.

**Note: There are currently mainly scripts for this figure that need to be consolidated. Running MASTER_task_tiling on its own will not generate the figure, at the moment need to run all of the below scripts to get the figure (but will consolidate into one script now that analysis is done):**

MASTER_task_tiling
MASTER_task_tiling_MIN
MASTER_task_tiling_control
MASTER_task_tiling_control_MIN
MASTER_task_tiling_location
MASTER_task_tiling_location_MIN
MASTER_task_tiling_location_control
MASTER_task_tiling_location_control_MIN
MASTER_task_tiling_corrCoeff
MASTER_task_tiling_figure


Figure 9:

To get Figure 9 from the paper, run MASTER_GLM_NP.

**Note: Like MASTER_GLM, running MASTER_GLM_NP in its current form will also generate several figures that are not used in Figure 9 (leftover from when developing the pipeline).**


Figure 10:

To get Figure 10 from the paper, run MASTER_population_average selecting the nosepoke-aligned option.

**Note: MASTER_population_average is currently only set up to generate Figure 7, but will add a switch command that detects which time-aligned data to pull out.** 


Figure 11:

To get Figure 11 from the paper, run MASTER_task_tiling selecting the nosepoke-aligned option.

**Note: Like with MASTER_population_average, MASTER_task_tiling is currently only set up to generate Figure 8, but will had the switch command, or alternatively a new script.**
