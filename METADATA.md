## Metadata

Each recording session contains a metadata file. Contained within this .mat file is the cell array *TrialInfo* that contains all the trial-related information for he first (TrialInfo{1,1}) and second block (TrialInfo{1,2}) of a session. Below is a brief summary of each subfield:

cue_type: whether it is a *light* or *sound* block.

cue_ID: which of the two cues within a block was presented during that trial.

rewarded: whether the current trial was a reward-available (‘1’) or reward-unavailable (‘0’) trial.
trialT: time of start of trial.

photosensorID: which arm of the track the trial occurred on.

nosepokeT: time of start of nosepoke (‘0’ if no nosepoke that trial).

nosepokeID: which receptacle of the track the trial occurred on.

trial_to_nosepokeT: length of time from trial initiation to nosepoke.

unnosepokeT: time of end of nosepoke.

nosepoke_length: length of nosepoke.

unnosepoke_to_trialT: time from end of nosepoke to initiation of following trial.

approached: whether the rat approached (‘1’) or skipped (‘0’) that trial.

quick_approach: trials in which the rat nosepoked for less than the 1 s delay from nosepoke start to outcome receipt.

trial_length: time from trial start to trial start.

backwards: if the rat ran backwards in between any trials.

trial_length_analysis: time from trial start to nosepoke (approach trials) or trial start (skip trials).

offsetT: used with *dataPoint* structure (see below).

Additionally, the *dataPoint* structure contains timestamps that the various scripts use to align behavioral and neural data for each trial at cue-onset (metadata.dataPoint.Trials) and subsequent nosepoke (metadata.dataPoint.Nosepokes).
