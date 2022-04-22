# wm_flex
Spatial working memory task with random and predictable responses. 

## ExptScripts
Code to generate stimulus sequences and run our fMRI experiments. Includes main task as well as all localizers.
Behavioral data from the experiments can be found in our OSF repository.
## Preprocessing
### PreprocScripts
All preprocessing code, including unwarping of raw data, motion correction/registration/detrending, ROI identification. We provide the fully pre-processed data in our OSF repository, so this code is not needed to reproduce our analyses, but it is provided here for convenience.
### MakeSamples
Code to create "SampleFile" and "TimingFile" for each subject - these files are included in our OSF repository. <br>
Code to create "avgSignal_XX" files for each subject and each task - these files need to be re-created by running the code here.
## Analysis
Code to reproduce all main analyses of behavior as well as fMRI data. Also includes additional analyses not included in our paper.

## To access data
Preprocessed data files can be downloaded from our Open Science Framework repository, at https://osf.io/5pk2z/

## Other notes
- In all code here, subjects 2-7 correspond to the 6 subjects included in the final analyses. Subject 1 participated in an early pilot version of the experiment and did not complete the final version of the experiment.
- This code also does preprocessing and some initial analysis for ROIs in the frontoparietal Multiple-Demand network, which we did not analyze or include in our paper. 
- Also included here is code for a button pressing task with a delayed response (referred to as DWM Loc), which is not analyzed or included in our paper.
