function ICA_template_pipeline(wMCM_flag, fslpath, afnipath, antspath, fmriprepdir, firstlevdir, templatepath, comp_sel_thresholds, numcomps, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel, mMCM_maskpath)

if(nargin==0)
    wMCM_flag = 1;
    spmpath = '/nas/data/app/spm12/';
    fslpath = '/usr/local/software/FSL-6.0.3/fsl';
    afnipath= '/usr/local/software/afni-2021_01';
    antspath = '/hct/fmriprep/install/bin/antsApplyTransforms';
    icasmoothkernel = 2.49;
    fmriprepdir='Data/';
    comp_sel_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 1.96];
    numcomps = {'auto', '20', '30', '40', '50', '60'};
    templatepath = 'Templates/HMAT_pre_post_central.nii';
    firstlevdir = 'FirstLevel/';
    comp_val_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];
    activ_val_thresholds = [ 0.05 0.01 0.005 0.001];
    taskname = 'motor';
end

addpath(spmpath);
 if wMCM_flag
    wMCM_component_generation(fslpath, fmriprepdir, numcomps, icasmoothkernel);
    wMCM_creation(antspath, afnipath, fmriprepdir, templatepath, comp_sel_thresholds, numcomps);
    wMCM_task_validation(fmriprepdir, comp_val_thresholds, activ_val_thresholds, taskname, firstlevdir);
 end
 end
