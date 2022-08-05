# EEG data processing codes listed below 
Codes associated with the analysis of auditory tones during smartphone usage


Auditory tones were created using: AudioTone.txt

Dependencies for EEG analysis: EEGlab, LIMO EEG, FOOOF

Prior to using the codes for ERP and CWT, the data was pre-procssed in the following manner:See methods section of the accompanying paper, inclduing the use og gettechnicallycleanedEEG.m

The following code was used for creating LIMO files for ERP (Level 1) : gatherEEGLIMOAUDIO_DataPool_Indoor_Out_ALLSubjects.m


The following code was used for creating LIMO files based on CWT (Level 1) (note, can take time and high (>180 GB) RAM workstation recommended): 
CWT alone:  gatherEEGLIMOAUDIO_DataPool_Indoor_Out_cwt_network
FOOOF with CWT:  gatherEEGLIMOAUDIO_DataPool_Out_cwt_fooof.m (was re-run for august revisions using: gatherEEGLIMOAUDIO_DataPool_Out_cwt_fooof_rerun.m)
Welch's method: gatherEEGLIMOAUDIO_DataPool_Outdoor_Welsh.m
FOOOF with Welch's method: gatherEEGLIMOAUDIO_DataPool_Out_welch_fooof_network_uploadCopy.m (run for august revisions)

The following code was used for performng second level one-sample-t-tests (Level 2)  : gatherfollowupOnesample_cwt_outdoor

The following code was used for creating the result figures:
CWT/FOOOF results: create_figure_spectralprops
AEP results: create_figure_AEP
Remaining individual figures: create_figure_other


#August 2022, Added the pre-processing steps: 

1. Data filtering (0.1 to 70 Hz) and removing of bad channels according to impedence measurements (using getdataorganised.m, a multipurpose data organiser used within CODELAB)
2. ICA was run using RUN_ICA.m
3. The missing channels were imputed using, re-referenced to the average, and blinks removed using (gettechnicallycleanedEEG.m)
4. For ERP analysis the data was further filtered between 0.1 and 45 Hz.




August 2022, Added the additional codes used for welch method: Indicate above in brackets as 'august revisions'. 
