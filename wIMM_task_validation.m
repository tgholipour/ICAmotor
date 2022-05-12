function wIMM_task_validation(fmriprepdir, thresholds, activation_thresh, taskname, firstlevdir, afnipath)
subids = dir(fmriprepdir);
subids = {subids.name};
subids = subids(3:length(subids));
icaactsumvarnames={'SubjectID','actfile', 'componentfile', 'z', 'task','activation_thresh', 'activ_cutoff', 'percent_act_in_networkmask'};
icaacttable=table({'subjectID'},{'actfile'},{'componentfile'},0,{'task'}, 0, 0, 0,'VariableNames', icaactsumvarnames);
tablei = 1;
for i = 1:length(subids)
    subid = char(subids(i));
    disp(subid)
    component_filename = dir([fmriprepdir subid '/' subid '_final_ICA_map.nii']);
    if ~isempty(component_filename)
        component_filename = [component_filename.folder '/' component_filename.name];
        templatepath = [fmriprepdir subid '/template/template_native_resample.nii'];
        templatemask = spm_read_vols(spm_vol(templatepath));
        templatemask = templatemask>0;
        totaltemplatemask = sum(templatemask(:));
        
        for actthresh = activation_thresh
            for zi = thresholds
                disp(taskname);
                activationfile = dir([firstlevdir subid '/*' taskname '*/*/spmT_0001.nii']);
                if ~isempty(activationfile)
                    activationfolder = activationfile.folder;
                    activationfile = [activationfolder '/' activationfile.name];
                    if isfile(activationfile)
                        disp(activationfile)
                        disp(actthresh)
                        disp(zi)
                        activation = spm_read_vols(spm_vol(activationfile));
                        smoothcomp = spm_read_vols(spm_vol(component_filename));
                        if sum(size(activation)==size(smoothcomp))<3
                            resampleactivationfile = [activationfolder '/resample_spmT_0001.nii'];
                            resamplestr = [afnipath '/3dWarp -deoblique -gridset ' component_filename ' -overwrite -prefix ' resampleactivationfile ' ' activationfile];
                            system(resamplestr);
                            activation = spm_read_vols(spm_vol(resampleactivationfile));
                        end
                        smoothcompz = smoothcomp>zi;
                        
                        templateact = activation;
                        templateact(~templatemask) = 0;
                        [~, idx] = max(templateact(:));
                        [r,c,p] = ind2sub(size(templateact), idx);
                        xs=[max(r-1,1):min(r+1,size(templateact,1))];
                        ys=[max(c-1,1):min(c+1,size(templateact,1))];
                        zs=[max(p-1,1):min(p+1,size(templateact,1))];
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
                        icaacttable.sensitivity(tablei) = count;
                        
                        sigact = activation>1.65;
                        
                        activation(~templatemask) = 0;
                        smoothcompz(~templatemask) = 0;
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
                        end
                        activation = activation>actcutoff;
                        
                        intersect_activationz = smoothcompz & activation;
                        hitrate_activationz = sum(intersect_activationz(:))/sum(activation(:));
                        icaacttable.hitrate(tablei) = hitrate_activationz;
                       
                        icaacttable.SubjectID{tablei,1} = subid;
                        icaacttable.task{tablei,1} =taskname;
                        icaacttable.actfile{tablei,1} = activationfile;
                        icaacttable.componentfile{tablei,1} = component_filename;
                        icaacttable.activation_thresh(tablei) = actthresh;
                        icaacttable.activ_cutoff(tablei) = actcutoff;
                        icaacttable.z(tablei) = zi;
                        icaacttable.percent_act_in_networkmask(tablei) = totalsigactmask/totaltemplatemask*100;
                        tablei = tablei+1;
                    end
                end
                
            end
        end
    end
    
end
rmpath(genpath('/nas/data/app/spm12'));
writetable(icaacttable,'wIMM_validation.csv','Delimiter',',');
end
