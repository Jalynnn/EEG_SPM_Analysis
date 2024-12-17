# EEG_SPM_Analysis

# To-Do

* Check reference and ground - Maybe verify with EEGLab
* String warning with for loop for "conditionlabels" or "condition_labels" variables
* Create a baseline with the 15 second rest windows

# Updates

12/13/2024
* Prepare.m file is better separated by task and filled with explanation comments
* 1-second epoch code is working as expected with no dead code
* Reviewed channel selections (accel, ref, and gnd) and warnings to create a refined to-do list
* Removed participant data from GH and uploaded to the OneDrive
* Updated J_prepare.m files are found everywhere

12/9/2024
* Repository created with the most recent scripts and participant data (incorrect event_data.csv's)
* GIT LFS initialized with .gitattributes specifically with *.psd file test and *.dat files
* Almost all event_data.csv files have been updated except the following:

  * Problematic with Importing (N = 3)
    * P17
    * P24
    * P30

  * Problematic with No Timestamps (N = 1)
    * P26

  * Problematic with Dropped files (N = 9)
    * P08
    * P25
    * P27
    * P28
    * P29
    * P30
    * P31
    * P32
    * P36
