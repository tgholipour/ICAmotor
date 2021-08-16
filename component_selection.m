%%This script selects the ICA component most likely to represent the
%%functional network of the template.
function component_selection(fmriprepdir, templatepath, comp_sel_thresholds, numcomps)
%%User input
clear all;
if nargin==0
    fmriprepdir='Data/';
    comp_sel_thresholds = [1 1.5 2 2.5 3 1.96];
    numcomps = {'auto', '20', '30', '40', '50', '60'};
    templatepath = 'Templates/sample_template.nii';
end
%end
%% Setup - get names of subjects and set up table
addpath('/nas/data/app/spm12');
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};
icasummaryvarnames={'SubjectID', 'template', 'numcomp', 'topdici', 'topcomp', 'topcomphr', 'topcompfar'};
icatable=table({'subjectID'},{'template'},{'numcomp'},0,0,0,0, 'VariableNames', icasummaryvarnames);
tablei = 1;
cd(fmriprepdir);
[~, nm, ~] = spm_fileparts(templatepath);
%%
for subi = 1:length(subs)
    subid = char(subs(subi));
    disp(subid);
    tic;
    inputfilename = dir([[subid '/rest/'], '*_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    if ~isempty(inputfilename)
        if ~isdir([subid '/Template/'])
            mkdir([subid '/Template/']);
        end
        nativetemplatepath = [subid '/Template/template_native.nii'];
        if isfile([subid '/anat/' subid '_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5'])
            nativize_template_str = ['/hct/fmriprep/antsApplyTransforms -i ' templatepath ' -t ' subid '/anat/' subid '_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5 '  ' -r ' subid '/anat/' subid '_desc-preproc_T1w.nii.gz -o ' nativetemplatepath];
            if ~exist(nativetemplatepath, 'file')
                system(nativize_template_str);
            end
        elseif isfile([subid '/anat/' subid '_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5'])
            nativize_template_str = ['/hct/fmriprep/antsApplyTransforms -i ' templatepath ' -t ' subid '/anat/' subid '_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 '  ' -r ' subid '/anat/' subid '_desc-preproc_T1w.nii.gz -o ' nativetemplatepath];
            if ~exist(nativetemplatepath, 'file')
                system(nativize_template_str);
            end
        end
        templatedeob = [subid '/Template/template_native_deoblique.nii'];
        if ~isfile(templatedeob)
            system(['module load AFNI/2021_01; 3dWarp -deoblique -prefix ' templatedeob ' -gridset ' subid '/restmelodic_auto/mean.nii.gz ' nativetemplatepath]);
        end
        template = spm_read_vols(spm_vol(templatedeob));
        
        
        topdicis = zeros(length(numcomps),length(comp_sel_thresholds));
        topdici_ind = 1;
        topdici_thresh = zeros(length(numcomps),length(comp_sel_thresholds));
        for numi = numcomps
            outdir = [subid '/restmelodic_' char(numi)];
            if ~isfile([outdir '/melodic_IC.nii'])
                gunzip([outdir '/melodic_IC.nii.gz']);
            end
            nummelodic=spm_read_vols(spm_vol([outdir '/melodic_IC.nii']));
            
            dici_template=zeros(size(nummelodic,4),1);
            hr_template=zeros(size(nummelodic,4),1);
            far_template=zeros(size(nummelodic,4),1);
            noise_template = zeros(size(nummelodic,4),1);
            hr_far_diff_template = zeros(size(nummelodic,4),1);
            for compi = 1:size(nummelodic,4)
                zind = 1;
                for zi = comp_sel_thresholds
                    component=squeeze(nummelodic(:,:,:,compi));
                    component = component>zi;
                    intersect_template = component & template;
                    hitrate_template = sum(intersect_template(:))/sum(template(:));
                    falsealarm_template = component & ~template;
                    falsealarmrate_template = sum(falsealarm_template(:))/sum(~template(:));
                    hr_far_template = norminv(hitrate_template, 0, 1)-norminv(falsealarmrate_template, 0, 1);
                    dici_template(compi, zind) = hr_far_template;
                    hr_template(compi, zind) = hitrate_template;
                    far_template(compi, zind) = falsealarmrate_template;
                    noise_template(compi, zind) = sum(falsealarm_template(:))/sum(component(:));
                    hr_far_diff_template(compi, zind) = (hitrate_template - falsealarmrate_template);
                    
                    zind = zind+1;
                end
            end
            [b,i] = sort(dici_template, 'descend');
            [comp, freq,vals] = mode(i(:,1:length(comp_sel_thresholds)-1),2);
            
            if(length(vals(1))>1)
                maxdici = -Inf;
                maxcomp = 0;
                for vali = 1:length(vals{1})
                    valcomp = vals{1}(vali);
                    dici196 = dici_template(valcomp, length(comp_sel_thresholds));
                    if dici196>maxdici
                        maxdici = dici196;
                        maxcomp = valcomp;
                    end
                end
                topdicis(topdici_ind, 1) = maxcomp;
                topdicis(topdici_ind, 2) = maxdici;
                topdicis(topdici_ind, 3) = hr_template(maxcomp, length(comp_sel_thresholds));
                topdicis(topdici_ind, 4) = far_template(maxcomp, length(comp_sel_thresholds));
                topdicis(topdici_ind, 5) = noise_template(maxcomp, length(comp_sel_thresholds));
                topdicis(topdici_ind, 6) = hr_far_diff_template(maxcomp, length(comp_sel_thresholds));
                
                
            else
                topdicis(topdici_ind, 1) = comp(1);
                topdicis(topdici_ind, 2) = dici_template(comp(1), length(comp_sel_thresholds));
                topdicis(topdici_ind, 3) = hr_template(comp(1), length(comp_sel_thresholds));
                topdicis(topdici_ind, 4) = far_template(comp(1), length(comp_sel_thresholds));
                topdicis(topdici_ind, 5) = noise_template(comp(1), length(comp_sel_thresholds));
                topdicis(topdici_ind, 6) = hr_far_diff_template(comp(1), length(comp_sel_thresholds));
                
            end
            topdici_thresh(topdici_ind, :) = dici_template(topdicis(topdici_ind, 1), :);
            
            topdici_ind = topdici_ind+1;
        end
        [M,I] = max(topdicis(:,2));
        outdir = ['/hct/fmriprep/ICA/Data/' subid '/restmelodic_' char(numcomps(I))];

        if ~isfile([outdir '/mean.nii'])
            gunzip([outdir '/mean.nii.gz']);
        end
        outdir = ['/hct/fmriprep/ICA/Data/' subid '/restmelodic_' char(numcomps(I))];
        nummelodic = spm_read_vols(spm_vol([outdir '/melodic_IC.nii']));
        compheader = spm_vol([outdir '/mean.nii']);
        compheader.fname = [outdir '/' subid '_topcomponent_' nm '.nii'];
        compheader.private.dat.fname = compheader.fname;
        spm_write_vol(compheader, squeeze(nummelodic(:,:,:,topdicis(I, 1))));
        
        for thresholdsi = 1:length(comp_sel_thresholds)
            icatable.SubjectID{tablei,1} = subid;
            icatable.thresholds(tablei)=comp_sel_thresholds(thresholdsi);
            icatable.dici(tablei) = topdici_thresh(I, thresholdsi);
            icatable.topmeandici(tablei) = M;
            icatable.topcomp(tablei) = topdicis(I, 1);
            icatable.topmeanhr(tablei) = topdicis(I, 3);
            icatable.topmeanfr(tablei) = topdicis(I, 4);
            icatable.topmeannoise(tablei) = topdicis(I, 5);
            icatable.topmeanhrfardiff(tablei) = topdicis(I, 6);
            icatable.numcomp{tablei,1}= char(numcomps(I));
            tablei = tablei+1;
        end
    end
    toc;
end
writetable(icatable,[fmriprepdir 'ICA_summary' nm '.csv'],'Delimiter',',');
end