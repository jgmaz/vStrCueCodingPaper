%% set up paths
restoredefaultpath; % start with clean slate
GITHUB_PATH = 'C:\Users\mvdm\Documents\GitHub\';
addpath(genpath([GITHUB_PATH 'vandermeerlab\code-matlab\shared'])); % lab codebase
addpath(genpath([GITHUB_PATH 'vandermeerlab\code-matlab\toolboxes\MClust-3.5'])); % MClust-3.5
addpath(genpath([GITHUB_PATH 'vStrCueCodingPaper\Code'])); % paper repo

%% locate the data
% top-level folder that contains subdirs with data from each subject (R053, R056, R057, R060)
% replace this with wherever your data is located
DATA_ROOT = 'C:\data\vStrCueCoding\'; 

%% Behavior example learning curves [NOT COMPLETED]

%% Behavior summary performance [NOT COMPLETED]

%% Preprocess units (needed for Figures 3,4,5,6,7)
directory = DATA_ROOT;
destination = [DATA_ROOT 'spike data\']; % where to save .mat files
PETH_generation = 1; % generating PETHs used for Figures 3,5,6,7-supplement is time consuming. Switch to 0 to bypass this step.

genProcess(directory,destination,PETH_generation);

%% Cue-onset GLM (needed for Figure 4)
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % working directory
destination = 'E:\vStr-cuecoding\promoted\analysis files\'; % where to save .mat files
num_Shuffs = 2;
select_Epoch = 1; % 1 = cue-onset only, 2 = nosepoke & outcome, 3 = all time points

genGLM(directory,destination);
genGLMshuff(directory,destination,num_Shuffs);
GLM_window = genGLMstats(destination,num_Shuffs,select_Epoch);

%% Plot cue-onset GLM
% Figure 4A
which_plot = 1; % 1 = cue-onset GLM in Figure 4A; 2 = cue-onset GLM in Figure 4 supplement 1; 3 = NP & outcome GLM in Figure 7A,B; 4 = NP & outcome GLM in Figure 7 supplement 1A-D

plotGLM(GLM_window,which_plot);

% Figure 4 supplement 1
which_plot = 2; % 1 = cue-onset GLM in Figure 4A; 2 = cue-onset GLM in Figure 4 supplement 1; 3 = NP & outcome GLM in Figure 7A,B; 4 = NP & outcome GLM in Figure 7 supplement 1A-D

plotGLM(GLM_window,which_plot);

%% Cue-onset LDA (needed for Figure 4)
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % working directory
destination = 'E:\vStr-cuecoding\promoted\analysis files\'; % where to save .mat files
num_Iterations = 5; % number of times to perform 10x cross-validation
num_Shuffs = 5;

genLDA(directory,destination,num_Iterations);
genLDAshuff(directory,destination,num_Shuffs,num_Iterations);
Class_accuracy = genLDAstats(destination);

%% Plot cue-onset LDA
% Figure 4B
plotLDA(Class_accuracy);

%% Nosepoke & outcome GLM (needed for Figure 7)
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % working directory
destination = 'E:\vStr-cuecoding\promoted\analysis files\'; % where to save .mat files
num_Shuffs = 2;
select_Epoch = 2; % 1 = cue-onset only, 2 = nosepoke & outcome, 3 = all time points

genGLMNP(directory,destination);
genGLMNPshuff(directory,destination,num_Shuffs);
GLM_window_NP = genGLMstats(destination,num_Shuffs,select_Epoch);

%% Plot nosepoke & outcome GLM
% Figure 7A,B
which_plot = 3; % 1 = cue-onset GLM in Figure 4A; 2 = cue-onset GLM in Figure 4 supplement 1; 3 = NP & outcome GLM in Figure 7A,B; 4 = NP & outcome GLM in Figure 7 supplement 1A-D

plotGLM(GLM_window_NP,which_plot);

% Figure 7 supplement 1A-D
which_plot = 4; % 1 = cue-onset GLM in Figure 4A; 2 = cue-onset GLM in Figure 4 supplement 1; 3 = NP & outcome GLM in Figure 7A,B; 4 = NP & outcome GLM in Figure 7 supplement 1A-D

plotGLM(GLM_window_NP,which_plot);

%% Recoded coefficients (needed for Figure 4,7)
spike_directory = 'E:\vStr-cuecoding\promoted\spike data\'; % where spiking data resides
directory = 'E:\vStr-cuecoding\promoted\analysis files\'; % working directory
destination = 'E:\vStr-cuecoding\promoted\analysis files\'; % where to save .mat files
num_Shuffs = 10;

GLM_coeff = genRecode(spike_directory,directory,destination);
genRecodeShuff(directory,destination,num_Shuffs);

%% Plot recoded correlation matrices
% Figure 4D
which_plot = 1; % 1 = cue-onset correlation matrix in Figure 4D; 2 = cue feature correlation matrices in Figure 7D-F; 3 = NP & outcome correlation matrices in Figure 7 supplement 1E-F

plotRecode(GLM_coeff,which_plot);

% Figure 7D-F
which_plot = 2; % 1 = cue-onset correlation matrix in Figure 4D; 2 = cue feature correlation matrices in Figure 7D-F; 3 = NP & outcome correlation matrices in Figure 7 supplement 1E-F

plotRecode(GLM_coeff,which_plot);

% Figure 7 supplement 1E-F
which_plot = 3; % 1 = cue-onset correlation matrix in Figure 4D; 2 = cue feature correlation matrices in Figure 7D-F; 3 = NP & outcome correlation matrices in Figure 7 supplement 1E-F

plotRecode(GLM_coeff,which_plot);

%% Plot cue-onset examples
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % where spiking data resides

% Figure 3
plotExamples(directory)

% Figure 3 supplements 1,2
plotExamplesSUPP(directory)

%% Plot nosepoke examples
directory = 'E:\vStr-cuecoding\promoted\spike data\'; % where spiking data resides

% Figure 6
plotExamplesNP(directory)

% Figure 6 supplements 1,2
plotExamplesNPSUPP(directory)

%% Cue-onset distributed firing (needed for Figure 5)

%% Plot cue-onset distributed firing
% Figure 5

% Figure 5 supplement 1

%% Nosepoke distributed firing (needed for Figure 7 supplement 2)

%% Plot nosepoke distributed firing
% Figure 7 supplement 1

%% Plot scatterplot
% Figure 4 supplement 2