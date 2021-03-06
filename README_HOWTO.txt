1st step: import data & adopt a particular architecture
=======================================================
With the script by Christophe Pallier (to be run on Neurospin workstations)
1) create manually a parameter file, saved as 'list_subjects.txt' (example provided in EXAMPLEpreproc_paramfile.m)
NB: the session to analyze must be labelled with the prefix epi_ . The reason is that functional files will have the form XXXepi_YYY.nii (with XXX the preprocessing prefixes and YYY the specific file name, whatever the length) and realignment parameters will have the form rp_aepi_YYY.txt; functional files will be matched with realignment files based on the YYY which is unique to a given session. Not using epi_ as a prefix in the import step will mess up the identification of this unique YYY tag in later steps.
2) run the command IN THE FOLDER WHERE DATA SHOULD BE IMPORTED:
# import_subjects_from_trio -a anat -f fMRI < list_subjects.txt
3) check that FieldMap for magnitude data has 'B0_1' in the file name; and that the phase data has 'B0_2' in the file name.

2nd step: preprocess data
=========================
With the help of the fmspm12batch
1) specify the parameter file (example provided in EXAMPLEpreproc_paramfile.m)
2) write and run SPM batches (and also Topup), given this parameter file, with fmspm12batch_preproc('mypreprofile')
3) check the movement parameters with fmspm12batch_replotRealign
4) execute the fmspm12batch_preproc_newArchi to copy files into a new architecture.
5) to extract surface (currently not in a batch): 
	i) Normalize the extracted white & grey matter: Batch > Normalize Write
		> deformation: y
		> images: c1 & C2
		> resolution [1.5 1.5 1.5]
	ii) extract surface: Render > Extract Surface

NB: the routine retreive most acquisition parameters (TR, slice timing ...) automatically, from the DICOM header. 

3rd step: first level analysis
==============================
1) specify the parameter file, modelname_paramfile.m
2) create the 'regressor script', modelname_multicond.m
3) create the 'contrast file', modelname_contrast.m (optional if you don't want to run contrasts)
4) run the analysis with fmspm12batch_1stlevel.m

4th step: 2nd level analysis
============================
1) specify the parameter file, modelname_paramfile.m
2) run the analysis with fmspm12batch_2ndlevel.m

BEWARE: CHECK THAT THESE ROUTINES ARE STILL VALID WITH SUBJECT NUMBER >9 (pb of 01 1 10 order...)

List of functions
=================
To run the preprocessing:
fmspm12batch_preproc.m
fmspm12batch_preproc_sf_make1job1sub.m
fmspm12batch_AddTopupCorrection_job1sub.m
fmspm12batch_preproc_GetSliceTiming_NS.m
fmspm12batch_preproc_newArchi.m

To run the 1st level analysis:
fmspm12batch_1stlevel_contrast.m
fmspm12batch_1stlevel_specify_sf_make1job1sub.m

To run the 2nd level analysis:
fmspm12batch_2ndLevelCorr.m
fmspm12batch_2ndlevel.m
fmspm12batch_2ndlevel_WriteBatch_Classicalttest.m

To run jobs:
fmspm12batch_run1job.m
fmspm12batch_runParalleljobs.m

Other independent tools:
> fmspm12batch_CheckContrastSpecif.m
> fmspm12batch_CheckROI.m
> fmspm12batch_CompareRes.m
> fmspm12batch_GetClusterCoord.m
review results from multiple subjects
> fmspm12batch_multipleIndiv.m
> fmspm12batch_multipleMIP.m
> fmspm12batch_multipleRender.m
> fmspm12batch_multipleSlices.m
PhyPI
> fmspm12batch_PhyPI_script.m
> fmspm12batch_PhyPI_v2_script.m
fmspm12batch_replotRealign.m
fmspm12batch_reviewDesign.m

SHOULD BE IMPROVED
==================
-> make the phase encoding direction an input parameter: matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.ewrap (currently it is necessarily Y).


