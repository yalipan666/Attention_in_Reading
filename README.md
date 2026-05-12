# Attention_in_Reading

This repository contains the code used for a co-registered eye tracking and MEG study, investigating attention in reading with rapid invisible frequency tagging (RIFT), together with the MATLAB, R, and Psychtoolbox scripts used for stimulus presentation, data processing pipeline (preprocessing, epoching, coherence analysis, behavioural analysis), and figure generation.

The pre-print manuscript is: https://www.biorxiv.org/content/10.1101/2025.05.27.656336v1.full.pdf

The raw datasets associated with this project are: https://figshare.com/projects/Attention_in_Reading/250316

For a brief introduction on how to display a visual stimulus in RIFT and how to analyse tagging response signals, please see the "RIFT_tutorial.pdf" in this repo. 

## Repository structure

### `PTB_codes/`
Psychtoolbox experiment code (PTB-3) and stimulus materials. This folder contains the main experiment script, EyeLink support files, and several `.mat` files with experimental sentences and comprehension questions.

### `Analyse_codes/`
Analysis code for MEG, eye-tracking, behavioural metrics, coherence analysis, group-level aggregation, correlation analyses, and figure generation.

### Top-level files
- `OccipSens_full.mat`: a supporting MATLAB file used by the analysis pipeline.
- `.gitignore`: ignores shared or generated folders.
- `README.md`: short repository summary on GitHub.
- `RIFT tutorial.pdf`: short description of the RIFT stimuli presentation and tagging response analysis.

## Processing and analysis workflow

A typical workflow in this codebase is:

1. **Prepare experiment metadata** with `S1_Get_ExpInfo.m`.
2. **Preprocess one participant** with `S2_PreProcessing.m`.
3. **Reject ICA artifacts and epoch the data** with `S3_Get_all_epoches.m`.
4. **Compute single-subject tagging coherence** with `Coh_RIFT.m` or `Coh_RIFT_parallel.m`.
5. **Aggregate subjects at group level** with `Group_Coh.m, `TaggingResponseCurve.m`, or `Group_Coh_parallel.m`.
6. **Eye movement data analysis and statistical tests** with `Ana_behavioral.m`.
7. **Relate tagging coherence to reading behavior** with `Corr_Tag_EM.m`.
8. **Single-trial tagging response** with `Coh_RIFT_parallel_ControlAnalysis.m`

## `Analyse_codes/` script guide

### `S1_Get_ExpInfo.m`
Creates the central experiment-information structure (`ExpInfo`). It stores subject IDs, PTB file names, EyeLink file names, participant metadata, bad sensors, sentence/frequency version assignments, eye-tracking side, MRI codes, trigger values, condition labels, and the event-header definition. It saves this structure to `Analyse_data/ExpInfo.mat`.

### `Get_Paths.m`
Returns the project root directory for either server or local execution, adds the analysis folder to the MATLAB path, and loads a diverging colourmap used in plotting.

### `S2_PreProcessing.m`
Runs the subject-level preprocessing step. It loads raw MEG FIF files, creates the MEG header/data/trigger objects, applies basic FieldTrip preprocessing for ICA, runs ICA on the downsampled data, parses EyeLink `.asc` files, builds `EyeData`, and constructs the trial-level `Event` table that aligns eye movements with MEG triggers.

### `ICA4rawdata.m`
Downsamples the data to 200 Hz and runs `ft_componentanalysis` with the `runica` method. This is used to identify artifact components for later rejection.

### `Get_MEGData.m`
Loads raw MEG data and trigger information from the FIF files and returns the header, continuous data matrix, and MEG trigger table. This is the MEG-data import helper used by preprocessing.

### `Get_EyeData.m`
Parses EyeLink ASCII output, extracts fixations, saccades, blinks, and sentence triggers, maps fixations to words using word-location coordinates, and stores all extracted eye-movement events in a structured format.

### `Get_Event.m`
Combines EyeLink timing with MEG trigger timing to build the main per-fixation event table. It computes word location, location relative to the target word, saccade duration before fixation, fixation onset in MEG time, first-pass status, next/previous word order, condition, and pupil size.

### `Get_Epoch.m`
Epoching helper used throughout the analysis pipeline. It is called to extract MEG segments around event triggers and is used for baseline epochs as well as reading-epoch extraction.

### `S3_Get_all_epoches.m`
Performs ICA component rejection and then builds the analysis-ready epochs. It creates baseline epochs, extracts reading epochs around fixation onsets or fixation offsets, applies trial selection rules, and saves the resulting FieldTrip epoch structures.

### `Ana_behavioral.m`
Computes behavioural summaries from the event files. It derives gaze-duration measures, total gaze duration, regression probability, and condition-wise behavioural metrics, then saves them to the behavioural results folder.

### `Coh_RIFT.m`
Single-subject RIFT coherence analysis. It loads the reading epochs, removes trials without valid photodiode signals, selects the tagging channels, computes time-resolved coherence in the tagging-frequency range, and stores the result in a subject-level coherence structure.

### `Coh_RIFT_parallel.m`
Parallelised version of the single-subject coherence analysis. It performs the same core RIFT coherence computation as `Coh_RIFT.m`, but is set up for batch or cluster execution.

### `Coh_RIFT_parallel_ControlAnalysis.m`
Parallel control-analysis version of the coherence pipeline. This version is intended for alternative analysis settings or control comparisons derived from the same tagging framework. Please note that, in this analysis, tagging responses were quantified using power rather than coherence. This is because coherence is a cross-trial measure, whereas power analysis can be performed at the single-trial level.

### `Group_Coh.m`
Group-level aggregation script for coherence results. It loads the per-subject `TagCoh` outputs, aligns the two frequency versions used across participants, collects significant tagging sensors, and saves the combined group structure to `TagCoh.mat`.

### `Group_Coh_parallel.m`
Parallel/group-scale version of the coherence aggregation workflow. It supports running the subject aggregation and group setup in a batch-friendly way.

### `TaggingResponseCurve.m`
Builds and plots the tagging-response curve across word positions. For each subject and tag frequency, it computes coherence around each word location, averages across the selected occipital tagging sensors and time window, and then summarises the curve at the group level.

### `Corr_Tag_EM.m`
Relates tagging-coherence measures to eye-movement behaviour. It calculates correlations between tagging coherence differences and reading speed, first-fixation duration, and other behavioural summaries, and stores the resulting statistics in a correlation structure.

### `Correlations_plot.py`
Python plotting helper for correlation results. It is used to visualise the relationships produced by the MATLAB correlation analyses.

### `SingleTrl_PowDiff.R`
R analysis script for single-trial power-difference data. It runs linear mixed-effects models and creates summary plots for power differences between tagging-related conditions and noise conditions.

## `PTB_codes/` contents

### `Read_FullAttention_2wrd.m`
Main Psychtoolbox experiment script. It sets the subject/session parameters, screen geometry, EyeLink options, tag-frequency assignment, stimulus timing, and display settings, then launches the reading experiment and writes the session output files.

### `Eyelink/`
Support files for EyeLink calibration, validation, recording, and file transfer during the experiment.

### Stimulus and condition files
- `Question.mat`: question material used in the task.
- `SentMat_1.mat`, `SentMat_2.mat`: sentence and word-layout materials for the two sentence versions.
- `TarWodLocFreq_1.mat`, `TarWodLocFreq_2.mat`: target-word location and frequency information for the two experimental versions.

### `io64.mexw64`
Windows I/O library used for low-level hardware communication, typically for trigger sending or parallel-port control.

## Outputs produced by the pipeline

The analysis code writes intermediate and final files into `Analyse_data/` and `Results/`, including:
- `ExpInfo.mat`
- `data.mat`, `hdr.mat`, `Trigger_MEG.mat`
- `EyeData.mat`, `Event.mat`
- `ica.mat`, `data_icaclean.mat`
- `epoch_BL_Cross.mat`, `epoch_WrdOn.mat`
- `TagCoh.mat`, `TagCurve.mat`
- `BehaData.mat`, `Corr.mat`

## Notes

- The codebase is built around **MEG + eye-tracking reading experiments** with **RIFT** tagging.
- The analysis pipeline depends heavily on **MATLAB/FieldTrip**, with **R** used for the single-trial mixed-effects analysis and **Python** used for some plots.
- The analysis scripts assume a fixed folder structure under `Full_Attention/` and use separate server/local root paths in `Get_Paths.m`.
