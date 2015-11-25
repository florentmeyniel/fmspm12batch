1st step: import data & adopt a particular architecture
=======================================================
With the script by Christophe Pallier (to be run on Neurospin workstations)
1) create manually a parameter file, saved as 'list_subjects.txt' (example provided in EXAMPLEpreproc_paramfile.m)
2) run the command IN THE FOLDER WHERE DATA SHOULD BE IMPORTED:
# import_subjects_from_trio -a anat -f fMRI < list_subjects.txt
3) check that FieldMap for magnitude data has 'B0_1' in the file name; and that the phase data has 'B0_2' in the file name.

2nd step: preprocess data
=========================
With the help of the fmspm12batch
1) specify the parameter file (example provided in EXAMPLEpreproc_paramfile.m)
2) write and run SPM batches, given this parameter file, with fmspm12batch_preproc('mypreprofile')
2-bis) if interested, run the top-up correction in addition, with fmspm12batch_AddTopupCorrection('mypreprofile')
3) check the movement parameters with fmspm12batch_replotRealign
4) execute the fmspm12batch_preproc_newArchi to copy files into a new architecture.
5) to extract surface (currently not in a batch): 
	i) Normalize the extracted white & grey matter: Batch > Normalize Write
		> deformation: y
		> images: c1 & C2
		> revolution [1.5 1.5 1.5]
	ii) extract surface: Render > Extract Surface

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

dans fmspm12batch_preproc_sf_make1job1sub:
% TO IMPROVE: l. 166 : faire une boucle sur les session pour récupérer la première image automatiquement!!!
