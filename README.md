## Flexible utilization of spatial- and motor-based codes for the storage of visuo-spatial information.
Code related to our 2022 paper: <br>
https://doi.org/10.7554/eLife.75688

### To access data:

Preprocessed fMRI data files can be downloaded from our Open Science Framework repository, at: https://osf.io/te5g2/ <br><br>
That repository includes files "SampleFile_S0X.mat" and "TimingFile_S0X.mat" for each subject, as well as behavioral data (DataBehavioral.zip). It also includes the saved results of our decoding analyses (i.e. accuracy and shuffled label accuracy) in files "spat_decoding_results.zip", "resp_decoding_results.zip" and "bound_decoding_results.zip". <br> <br>
To re-create the files with names like "MainTaskSignalByTrial_S0X.mat", which are directly used by our decoding analyses, you should run the script at:
Preprocessing/MakeSamples/get_trial_avg_signal_all.m

### Contents of this repository:

- <b>ExptScripts</b>
  - Code to generate stimulus sequences and run our fMRI experiments. Includes main task as well as all localizers.
      Behavioral data from the experiments can be found in our OSF repository. </li>
      
- <b> Preprocessng </b> 
  - <b> PreprocScripts</b>
    - All preprocessing code, including unwarping of raw data, motion correction/registration/detrending, ROI identification. 
    - We provide the fully pre-processed data in our OSF repository, so this code is not needed to reproduce our analyses, but it is provided here for convenience.    
  - <b> MakeSamples </b>
    - Code to create "SampleFile" and "TimingFile" for each subject - these files are included in our OSF repository. 
    - Code to create "XX_SignalByTrial_XX" files for each subject and each task - these files need to be re-created by running the code here (get_trial_avg_signal_all.m)

- <b> Analysis </b>
  - <b> AnalyzeBehavior </b>
    - Code to analyze subject behavior on our main task as well as all localizer tasks. To run these, you'll need the behavioral data (located in our OSF repository as DataBehavior.zip). 
    - This should be unzipped into a folder called DataBehavior, on same level as the Analysis directory.
    
  - <b> Univariate </b>
    - Code to perform deconvolution and make plots of deconvolved BOLD signal for each ROI.
  - <b> Decode_space </b>
    - Code to perform multivariate decoding of remembered spatial position. 
    - The files starting with "Trn..." in this folder perform the decoding analyses, and files starting with "plot..." make the plots and perform statistical testing. 
    - You can re-run the decoding analyes from scratch using this code, or load the results of decoding from our OSF repository (spat_decoding_results.zip) and then start from the plotting stage. Results of decoding should be placed in a folder called "Decoding_results" inside the "Decode_space" directory.   
  - <b> Decode_response </b>
    - Code to perform multivariate decoding of upcoming motor actions (button press). 
    - The files starting with "Classify..." in this folder perform the decoding analyses, and files starting with "plot..." make the plots and perform statistical testing. 
    - You can re-run the decoding analyes from scratch using this code, or load the results of decoding from our OSF repository (resp_decoding_results.zip) and then start from the plotting stage. Results of decoding should be placed in a folder called "Decoding_results" inside the "Decode_response" directory. 
  - <b> Decode_boundary </b>
    - Code to perform multivariate decoding of preview disk boundary. 
    - Classify_boundary.m does the decoding analyses, and plotClassResults_Boundary.m makes the plots and performs statistical testing. 
    - You can re-run the decoding analyes from scratch using this code, or load the results of decoding from our OSF repository (bound_decoding_results.zip) and then start from the plotting stage. Results of decoding should be placed in a folder called "Decoding_results" inside the "Decode_boundary" directory. 
  - <b> stats_code </b>
    - General use code for statistical testing and for training/testing decoders.
  - <b> plotting_utils </b>
    - General use code for making plots.
    
### Other notes:

  - In all code here, participants 2-7 correspond to the 6 subjects included in the final analyses. Subject 1 participated in an early pilot version of the experiment and did not complete the final version of the experiment.
  - This code also does preprocessing and some initial analysis for ROIs in the frontoparietal Multiple-Demand network, which we did not analyze or include in our paper. 
  - Also included here is code for a button pressing task with a delayed response (referred to as DWM Loc), which is not analyzed or included in our paper.
  - Some files in this repo will have the string "OriSpin2" - this was the original codename for this experiment, but doesn't mean anything otherwise.
  - Sometimes in the code the conditions are referred to as "Predictable" and "Random" - these are older names for the conditions that we called "Informative" and "Uninformative"
  - In the behavioral data directories (contained in DataBehavior.zip), there are many sessions for each participant, even though only 3 sessions were performed as part of this experiment. The remaining sessions contain data from localizer task runs that were performed during separate sessions (usually, earlier in time as part of other experiments).


Please contact mmhender@cmu.edu with any questions or concerns. Thanks!
