# vStrCueCoding data

Preprocessed data used for Gmaz,
Carmichael & van der Meer, "Persistent coding of outcome-predictive cue features in the rat nucleus accumbens" (2018) ([preprint](https://www.biorxiv.org/content/early/2018/08/27/300251)).

This data set contains data from four subjects (R053, R056, R057,
R060), organized into daily recording sessions (one session per
subfolder).

Each recording session contains the following:
  * `ExpKeys.m` file, one per session, describing some basic properties of that session.
  * `*.t` files, one per putative single neuron, containing spike times. Load with `LoadSpikes()`.
  * `metadata.mat` file, one per session, containing task events such as cue presentation times, nosepokes, and reward delivery times.

### Description of ExpKeys fields

### Description of metadata fields


