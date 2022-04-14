# EEG data processing codes listed below 
Codes associated with the analysis of auditory tones during smartphone usage


Auditory tones were created using: AudioTone.txt

Dependencies for EEG analysis: EEGlab, LIMO EEG, FOOOF

Prior to using the codes for ERP and CWT, the data was pre-procssed in the following manner:See methods section of the accompanying paper, inclduing the use og gettechnicallycleanedEEG.m

The following code was used for creating LIMO files for ERP (Level 1) : gatherEEGLIMOAUDIO_DataPool_Indoor_Out_ALLSubjects.m


The following code was used for creating LIMO files based on CWT (Level 1) (note, can take time and high (>180 GB) RAM workstation recommended): 
CWT alone:  gatherEEGLIMOAUDIO_DataPool_Indoor_Out_cwt_network
FOOOF:  gatherEEGLIMOAUDIO_DataPool_Out_cwt_fooof


The following code was used for performng second level one-sample-t-tests (Level 2)  : gatherfollowupOnesample_cwt_outdoor
