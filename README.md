# ICA motor

This pipeline is meant to extract the dimension of resting state fMRI data representative of a functional network, given a template representing said network through independent component analysis (ICA). The method has 3 steps:

1) wIMM Component Generation - Generate ICA components at varying degrees of dimension reduction.
2) wIMM Creation - In case the network of interest has been split between components, combine all representative components into final map.
3) wIMM Task Validation - Calculating hit rate of top component and first level activation maps of task fMRI data, meant to elicit the same functional network.


# Installation
Download and unzip repository. The pipeline can be run in parts, or can be run all together from the ICA_template_pipeline script. Default values have been provided, and can be changed to fit the environment in which it is being run. 

SPM, AFNI, FSL, and antsApplyTransforms are required before running this pipeline and can be installed for free at the links below.
https://www.fil.ion.ucl.ac.uk/spm/software/download/
https://afni.nimh.nih.gov/download
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
https://sourceforge.net/projects/advants/files/


Before running, preprocess anatomical and resting state fMRI data through fMRIPrep. If planning on utilizing validation aspect of pipeline, generate first-level activation T-map of task-fMRI data in SPM.
