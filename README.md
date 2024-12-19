# EEG_SPM_Analysis

# To-Do

* Check reference and ground
* * Default was created by LiveAmp Connector, so none of our recorded data will explicitly show which channels are which. Searching for a wrap-around - in the lab later today, I'm going to look at the software we used to pick up signal quality.
* Create an automatic baseline instead of manual input.

# Updates

12/19/2024
* The order of steps was rearranged due to upcoming errors with renaming 'ACC' channels.
* Several attempts at adding a baseline correction allowed me to find cases of variables being overwritten, so this was also fixed, and several D and S variables exist.
* Baseline was fixed as it was already being successfully added in with .Dbaseln, but essentially, eeg_spm_epochs had its type of baseline that wasn't being added in, so it was quieted and set equal to 0.
* Using a baseline of 15 seconds is currently done manually, and this aims to be fixed by 12/20. However, when working on this, I had to start using P11 instead of P10 since P10 had no recorded rest periods.
* Looking at the string warning was a rabbit hole that circled back to our first error of unmatched trial sizes. This was due to the trl file having one line as anticipated but still using seven labels (Condition 1 x6 and Condition 1_TEST).
* After all these changes, I started looking at the other two files and found that no files were being rejected!

12/13/2024 - 12/17/2024
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

# Resources

https://jsheunis.github.io/2018-06-28-spm12-matlab-scripting-tutorial-1/

https://jsheunis.github.io/2018-06-28-spm12-matlab-scripting-tutorial-4/

https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/MEEG/mmn/

https://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf

https://github.com/neurodebian/spm12/blob/b86f0712da80daa4fae7939eeebebe0b2629cbc0/toolbox/MEEGtools/spm_eeg_tms_correct.m#L68

https://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf

https://github.com/brain-products/LSL-LiveAmp
