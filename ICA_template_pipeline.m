%This pipeline requires spm12 which can be downloaded at https://www.fil.ion.ucl.ac.uk/spm/software/download/
function ICA_template_pipeline(wMCM_flag, mMCM_flag, fslpath, afnipath, antspath, fmriprepdir, firstlevdir, templatepath, comp_sel_thresholds, numcomps, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel, mMCM_maskpath)

if(nargin==0)
    spmpath = '';
    fslpath = '';
    afnipath='';
    antspath = '';
    smoothkernel = 2.49;
    mMCM_maskpath = '';
    fmriprepdir='Data/';
    comp_sel_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 1.96];
    numcomps = {'auto', '20', '30', '40', '50', '60'};
    templatepath = 'Templates/sample_template.nii';
    firstlevdir = 'FirstLevel/';
    comp_val_thresholds = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];
    activ_val_thresholds = [ 0.05 0.01 0.005 0.001];
    taskname = 'motor';
end
addpath(spmpath);
 if wMCM_flag
    wMCM_component_generation(fslpath, fmriprepdir, numcomps, smoothkernel);
    wMCM_creation(antspath, afnipath, fmriprepdir, templatepath, comp_sel_thresholds, numcomps);
    wMCM_task_validation(fmriprepdir, comp_val_thresholds, activ_val_thresholds, templatepath, taskname, firstlevdir);
 end
 if mMCM_flag
     mMCM_component_generation(fslpath, afnipath, antspath, fmriprepdir, mMCM_maskpath, smoothkernel);
     mMCM_creation(fmriprepdir, templatepath);
     mMCM_task_validation(fmriprepdir, templatepath, comp_val_thresholds, activ_val_thresholds, taskname, firstlevdir);
 end
 end
