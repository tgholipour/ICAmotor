function ICA_template_pipeline(maindir, fslpath, afnipath, antspath, fmriprepdir, firstlevdir, templatepath, comp_sel_thresholds, numcomps, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel, meanshiftclusterpath)
if(nargin==0)
    maindir = '/home/mkrishnamu/Downloads/ICAmotor-main/';
    spmpath = '/nas/data/app/spm12/';
    fslpath = '/usr/local/software/FSL-6.0.3/fsl';
    afnipath= '/usr/local/software/afni-2021_01';
    antspath = '/hct/fmriprep/install/bin/antsApplyTransforms';
    icasmoothkernel = 2.49; %Smoothing kernel to smooth fmriprep preprocessed rs-fMRI data in fsl
    fmriprepdir='Data/'; %directory with fMRIprep outputs
    comp_sel_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]; %list of thresholds used to select components that comprise the final map
    numcomps = {'auto', '20', '30', '40', '50', '60'}; %list of Total Number of Components used to select components that comprise the final map
    templatepath = 'Templates/HMAT_pre_post_central.nii';
    templateMNIspace = 'MNI152NLin6Asym';
    firstlevdir = 'FirstLevel/';
    comp_val_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];%List of ICA map thresholds used to compare to task activation
    activ_val_thresholds = [ 0.05 0.01 0.005 0.001];%List of task activation p value thresholds (assuming task activation is a t map)
    taskname = 'motor';
    meanshiftclusterfolderpath = '/hct/fmriprep/ICA/';
end
cd(maindir);
addpath(spmpath);
addpath(meanshiftclusterfolderpath);
wIMM_component_generation(fslpath, fmriprepdir, numcomps, icasmoothkernel);
wIMM_creation(antspath, afnipath, fmriprepdir, templatepath, comp_sel_thresholds, numcomps, templateMNIspace);
wIMM_task_validation(fmriprepdir, comp_val_thresholds, activ_val_thresholds, taskname, firstlevdir, afnipath);
end
