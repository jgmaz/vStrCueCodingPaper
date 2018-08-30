# vStrCueCodingPaper


Figure 2:

To get Figure 2 from the paper, first run MASTER_behavior_example to extract the variables needed for the example plots in B. Then run MASTER_behavior to extract the summary performance measures used in C. Then run MASTER_behavior_figure to generate Figure 2.

**Note: Access to behavioral data in MASTER_behavior_example is currently hardcoded in the script. Can save all the .mat files it uses and put them in a folder for the paper.**

**Note: In the current form of MASTER_behavior, each behavioral session was manually run in the first half of the script, before generated the summary values in the second half of the script. Will include a call to the .mat files containing each session in the final version.**


Recording pre-analysis:

All the neural analyses use variables generated from MASTER_main_workflow, whose purpose is to take the .t files for each cell, extract the spikes, sync with the behavioral data, and generate the firing measures used for subsequent analyses.


Figure 3:

To get Figure 3 from the paper, run MASTER_GLM_example_cells.



Figure 3 supplements 1 and 2:

To get Figure 3 supplements 1 and 2 form the paper, run MASTER_GLM_example_cells_expanded_for_SUPP.


Figure 4:

To get Figure 4 from the paper, run MASTER_GLM to perform the sliding window GLM, and MASTER_GLM_Shuff to generate the shuffled data. Run LDA_SCRIPT to perform the LDA, and LDA_shuff to generate the shuffled data. To get Figure 4A,B, run FIGURE script. To get Figure 4D run the RECODE script.

Figure 4 supplements 1 and 2:

To get Figure 4 supplement 1 from the paper, run FIGURE script after running MASTER_GLM as above.

To get Figure 4 supplement 2 from the paper, run SCATTER script.


Figure 5:

To get Figure 5 from the paper, run MASTER_task_tiling.

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



Figure 5 supplement 1:

To get Figure 5 supplement 1 from the paper, run FIGURE script.


Figure 6:

To get Figure 6 from the paper, run MASTER_GLM_NP_example_cells.



Figure 6 supplements 1 and 2:

To get Figure 6 supplements 1 and 2 form the paper, run MASTER_GLM_NP_example_cells_expanded_for_SUPP.


Figure 7:

To get Figure 7 from the paper, run MASTER_GLM_NP to perform the sliding window GLM, and MASTER_GLM_NP_Shuff to generate the shuffled data. To get Figure 4A,B, run FIGURE script. To get Figure 4D-F run the RECODE script.



Figure 7 supplements 1 and 2:


To get Figure 7 supplement 1 A-D from the paper, run FIGURE script after running MASTER_GLM_NP as above. To get Figure 7 supplement 1 E,F, run FIGURE script after running RECODE script as above.

To get Figure 7 supplement 2 from the paper, run MASTER_task_tiling_NP.
