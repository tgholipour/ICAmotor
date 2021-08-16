%This pipeline requires spm12 which can be downloaded at https://www.fil.ion.ucl.ac.uk/spm/software/download/
function ICA_template_pipeline(fslpath, fmriprepdir, firstlevdir, templatepath, comp_sel_thresholds, numcomps, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel)

if(nargin==0)
    fslpath = '';
    fmriprepdir='Data/';
    comp_sel_thresholds = [1 1.5 2 2.5 3 1.96];
    numcomps = {'auto', '20', '30', '40', '50', '60'};
    templatepath = 'Templates/sample_template.nii';
    firstlevdir = 'FirstLevel/';
    comp_val_thresholds = [1 1.5 2 2.5 3];
    activ_val_thresholds = [ 0.05 0.01 0.005 0.001];
    taskname = 'motor';
    icasmoothkernel = 6;
 end
    component_generation(fslpath, fmriprepdir, numcomps);
    component_selection(fmriprepdir, templatepath, comp_sel_thresholds, numcomps);
    component_validation(fmriprepdir, templatepath, firstlevdir, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel);
end
