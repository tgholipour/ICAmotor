function component_validation(fmriprepdir, templatepath, firstlevdir, comp_val_thresholds, activ_val_thresholds, taskname, icasmoothkernel)
clear all;
%% User Input
if nargin==0
    fmriprepdir = 'Data/';
    templatepath = 'Templates/sample_template.nii';
    firstlevdir = 'FirstLevel/';
    comp_val_thresholds = [1 1.5 2 2.5 3];
    activ_val_thresholds = [ 0.05 0.01 0.005 0.001];
    taskname = 'motor';
    icasmoothkernel = 6;
end
%% Setup - read in ica summary
addpath('/nas/data/app/spm12');
[~, nm, ~] = spm_fileparts(templatepath);
tic;
ica_topdici = readtable([fmriprepdir 'ICA_summary' nm '.csv']);
ica_topdici = ica_topdici(:,[1 3 5]);
ica_topdici = unique(ica_topdici);
subids = ica_topdici.SubjectID;
numcomps = ica_topdici.numcomp;
topcomps = ica_topdici.topcomp;
icaactsumvarnames={'SubjectID','actfile', 'componentfile','topcomp' 'network', 'z', 'task','activ_val_thresholds', 'activ_cutoff', 'mutual_info', 'dici', 'dice_coef', 'percent_act_in_networkmask'};
icaacttable=table({'subjectID'},{'actfile'},{'componentfile'},0,{'network'},0,{'task'}, 0, 0, 0, 0, 0, 0,'VariableNames', icaactsumvarnames);
tablei = 1;
cd(fmriprepdir);
%% Begin Validation
for i = 1:height(ica_topdici)
    subid = char(subids(i));
    disp(subid)
    
    numcomp = numcomps(i);
    if isnan(numcomp)
        numcomp = 'auto';
    end
    topcomp = topcomps(i);
    
    templatemaskfile = [subid '/Template/template_native_deoblique.nii'];
    templatemask = spm_read_vols(spm_vol(templatemaskfile));
    templatemask = templatemask>0;
    totaltemplatemask = sum(templatemask(:));
    
    
    component_filename = [subid '/restmelodic_' num2str(numcomp) '/' subid '_topcomponent_' nm '.nii'];
    output_filename = [subid '/restmelodic_' num2str(numcomp) '/' subid '_topcomponent_' nm '_smooth.nii'];
    spm_smooth(component_filename,output_filename , [icasmoothkernel,icasmoothkernel,icasmoothkernel]);
    
    for actthresh = activ_val_thresholds
        for zi = comp_val_thresholds
            disp(taskname)
            activationfile = dir([firstlevdir '/' subid '/' taskname '/*/spmT_0001.nii']);
            if ~isempty(activationfile)
                actfolder = activationfile.folder;
                activationfile = [activationfile.folder '/' activationfile.name];
                if isfile(activationfile)
                    activation = spm_read_vols(spm_vol(activationfile));
                    smoothcomp = spm_read_vols(spm_vol(output_filename));
                    smoothcompz = smoothcomp>zi;
                    
                    taskact = activation;
                    taskact(~templatemask) = 0;
                    [~, idx] = max(taskact(:));
                    [r,c,p] = ind2sub(size(taskact), idx);
                    xs=[r-1:r+1];
                    ys=[c-1:c+1];
                    zs=[p-1:p+1];
                    count=0;
                    for x=xs
                        for y=ys
                            for z=zs
                                if(smoothcompz(x,y,z)~=0)
                                    count=count+1;
                                    
                                end
                                
                            end
                        end
                    end
                    icaacttable.motorsensitivity(tablei) = count;
                    
                    actpos = activation(activation>0);
                    actcutoff = prctile(actpos,100-actthresh);
                    activationposthrsh = activation>actcutoff;
                    
                    maskcomp = smoothcomp;
                    maskcomp(~templatemask) = 0;
                    [maxv, idx] = max(maskcomp(:));
                    [r,c,p] = ind2sub(size(maskcomp), idx);
                    if(activationposthrsh(r,c,p)~=0)
                        icaacttable.motorspecificity(tablei) = 1;
                    else
                        icaacttable.motorspecificity(tablei) = 0;
                    end
                    
                    
                    sigact = activation>1.65;
                    
                    
                    activation(~templatemask) = 0;
                    smoothcomp(~templatemask) = 0;
                    sigact(~templatemask) = 0;
                    
                    totalsigactmask = sum(sigact(:));
                    
                    if actthresh==0.05
                        actcutoff = 1.65;
                    elseif actthresh==0.01
                        actcutoff = 2.35;
                    elseif actthresh==0.005
                        actcutoff = 2.612;
                    elseif actthresh==0.001
                        actcutoff = 3.15;
                    else
                        actpos = activation(activation>0);
                        actcutoff = prctile(actpos,100-actthresh);
                    end
                    activation = activation>actcutoff;
                    
                    spm=load([actfolder '/SPM.mat']);
                    spm_df=spm.SPM.xX.erdf;
                    p=1-cdf('t',actcutoff,spm_df);
                    icaacttable.activ_val_thresholds_p(tablei) = p;
                    
                    intersect_activationz = smoothcompz & activation;
                    hitrate_activationz = sum(intersect_activationz(:))/sum(activation(:));
                    falsealarm_activationz = smoothcompz & ~activation;
                    
                    notactmask = ~activation&templatemask;
                    
                    falsealarmrate_activationz = sum(falsealarm_activationz(:))/sum(notactmask(:));
                    hr_far_activationz = norminv(hitrate_activationz, 0, 1)-norminv(falsealarmrate_activationz, 0, 1);
                    icaacttable.hitrate(tablei) = hitrate_activationz;
                    icaacttable.falsealarmrate(tablei) = falsealarmrate_activationz;
                    
                    icaacttable.actpropofmask(tablei) = sum(activation(:))/sum(templatemask(:));
                    icaacttable.comppropofmask(tablei) = sum(smoothcompz(:))/sum(templatemask(:));
                    
                    
                    icaacttable.SubjectID{tablei,1} = subid;
                    icaacttable.task{tablei,1} = [actfolder];
                    icaacttable.actfile{tablei,1} = activationfile;
                    icaacttable.componentfile{tablei,1} = output_filename;
                    icaacttable.topcomp(tablei) = topcomp;
                    icaacttable.activ_val_thresholds(tablei) = actthresh;
                    icaacttable.activ_cutoff(tablei) = actcutoff;
                    icaacttable.z(tablei) = zi;
                    icaacttable.dici(tablei) = hr_far_activationz;
                    icaacttable.percent_act_in_networkmask(tablei) = totalsigactmask/totaltemplatemask*100;
                    tablei = tablei+1;
                end
            end
            
        end
    end
    
end
rmpath(genpath('/nas/data/app/spm12'));
cd('/hct/fmriprep/ICA/');
writetable(icaacttable,[fmriprepdir 'ICA_validation' nm '.csv'],'Delimiter',',');
toc;
end