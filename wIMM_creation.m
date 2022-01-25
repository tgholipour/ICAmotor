function wIMM_creation(antspath, afnipath, fmriprepdir, templatepath, thresholds, numcomps)
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};
subs = subs(3:length(subs));
thresholds = [thresholds 1.96];
for subi = 1:length(subs)
    disp(char(subs(subi)))
    subid = char(subs(subi));
    tic;
    if ~isfolder([fmriprepdir subid '/template/'])
        mkdir([fmriprepdir subid '/template/'])
    end
    inputfilename = dir([[fmriprepdir char(subs(subi)) '/func/'], '*_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    
    if ~isempty(inputfilename)
        nativetemplatepath = [fmriprepdir char(subs(subi)) '/template/template_native.nii'];
        if isfile([fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5'])
            nativize_template_str = [antspath ' -i ' templatepath ' -t ' fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5 '  ' -r ' fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_desc-preproc_T1w.nii.gz -o ' nativetemplatepath];
            if ~exist(nativetemplatepath, 'file')
                system(nativize_template_str);
            end
        elseif isfile([fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5'])
            nativize_template_str = [antspath ' -i ' templatepath ' -t ' fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 '  ' -r ' fmriprepdir char(subs(subi)) '/anat/' char(subs(subi)) '_desc-preproc_T1w.nii.gz -o ' nativetemplatepath];
            if ~exist(nativetemplatepath, 'file')
                system(nativize_template_str);
            end
        end
        templatedeob = [fmriprepdir char(subs(subi)) '/template/template_native_deoblique.nii'];
        templateresample = [fmriprepdir char(subs(subi)) '/template/template_native_resample.nii'];
        if ~isfile(templatedeob)
            system([afnipath '/3dWarp -deoblique -overwrite -prefix ' templatedeob ' ' nativetemplatepath]);
        end
        if ~isfile(templateresample)
            system([afnipath '/3dresample -overwrite -master ' fmriprepdir char(subs(subi)) '/restmelodicsmooth_auto/mean.nii.gz -prefix ' templateresample ' -input ' templatedeob]);
        end
        template = spm_read_vols(spm_vol(templateresample));
        template = template>0;
        
        
        topdicis = zeros(length(numcomps),length(thresholds));
        topdici_ind = 1;
        topdici_thresh = zeros(length(numcomps),length(thresholds));
        for numi = numcomps
            outdir = [fmriprepdir char(subs(subi)) '/restmelodicsmooth_' char(numi)];
            if ~isfile([outdir '/melodic_IC.nii'])
                gunzip([outdir '/melodic_IC.nii.gz']);
            end
            nummelodic=spm_read_vols(spm_vol([outdir '/melodic_IC.nii']));
            brainmaskfile = [outdir '/mask.nii'];
            if ~isfile (brainmaskfile)
                gunzip([brainmaskfile '.gz']);
            end
            
            brainmask = spm_read_vols(spm_vol(brainmaskfile));
            dici_template=zeros(size(nummelodic,4),1);
            hr_template=zeros(size(nummelodic,4),1);
            far_template=zeros(size(nummelodic,4),1);
            noise_template = zeros(size(nummelodic,4),1);
            hr_far_diff_template = zeros(size(nummelodic,4),1);
            for compi = 1:size(nummelodic,4)
                zind = 1;
                for zi = thresholds
                    component=squeeze(nummelodic(:,:,:,compi));
                    component = component>zi;
                    intersect_template = component & template;
                    hitrate_template = sum(intersect_template(:))/sum(template(:));
                    if hitrate_template==1
                        hitrate_template = (sum(template(:))-1)/sum(template(:));
                    elseif hitrate_template==0
                        hitrate_template = 1/sum(template(:));
                    end
                    falsealarm_template = component & ~template & brainmask;
                    outsidetemplate_insidebrainmask = ~template & brainmask;
                    falsealarmrate_template = (sum(falsealarm_template(:))/sum(outsidetemplate_insidebrainmask(:)));
                    if falsealarmrate_template==1
                        falsealarmrate_template = (sum(outsidetemplate_insidebrainmask(:))-1)/sum(outsidetemplate_insidebrainmask(:));
                    elseif falsealarmrate_template==0
                        falsealarmrate_template = 1/sum(outsidetemplate_insidebrainmask(:));
                    end
                    hr_far_template = norminv(hitrate_template, 0, 1)-norminv(falsealarmrate_template, 0, 1);
                    dici_template(compi, zind) = hr_far_template;
                    hr_template(compi, zind) = hitrate_template;
                    far_template(compi, zind) = falsealarmrate_template;
                    noise_template(compi, zind) = sum(falsealarm_template(:))/sum(component(:));
                    hr_far_diff_template(compi, zind) = (hitrate_template - falsealarmrate_template);
                    
                    zind = zind+1;
                end
            end
            [~,i] = sort(dici_template, 'descend');
            [comp, ~,vals] = mode(i(:,1:length(thresholds)-1),2);
            
            if(length(vals(1))>1)
                maxdici = -Inf;
                maxcomp = 0;
                for vali = 1:length(vals{1})
                    valcomp = vals{1}(vali);
                    dici196 = dici_template(valcomp, length(thresholds));
                    if dici196>maxdici
                        maxdici = dici196;
                        maxcomp = valcomp;
                    end
                end
                topdicis(topdici_ind, 1) = maxcomp;
                topdicis(topdici_ind, 2) = maxdici;
                topdicis(topdici_ind, 3) = hr_template(maxcomp, length(thresholds));
                topdicis(topdici_ind, 4) = far_template(maxcomp, length(thresholds));
                topdicis(topdici_ind, 5) = noise_template(maxcomp, length(thresholds));
                topdicis(topdici_ind, 6) = hr_far_diff_template(maxcomp, length(thresholds));
                
                
            else
                topdicis(topdici_ind, 1) = comp(1);
                topdicis(topdici_ind, 2) = dici_template(comp(1), length(thresholds));
                topdicis(topdici_ind, 3) = hr_template(comp(1), length(thresholds));
                topdicis(topdici_ind, 4) = far_template(comp(1), length(thresholds));
                topdicis(topdici_ind, 5) = noise_template(comp(1), length(thresholds));
                topdicis(topdici_ind, 6) = hr_far_diff_template(comp(1), length(thresholds));
                
            end
            topdici_thresh(topdici_ind, :) = dici_template(topdicis(topdici_ind, 1), :);
            
            topdici_ind = topdici_ind+1;
        end

        ncoltopdicis = size(topdicis,2);
        for topi = 1:size(topdicis,1)
            tncmelodic = spm_read_vols(spm_vol([fmriprepdir subid '/restmelodicsmooth_' char(numcomps(topi)) '/melodic_IC.nii']));
            tnctopcomp = squeeze(tncmelodic(:,:,:,topdicis(topi,1)));
            [~, idx] = max(tnctopcomp(:));
            [r,c,p] = ind2sub(size(tnctopcomp), idx);
            if topi ==1
                topdicis(topi,ncoltopdicis+1) = 0;
            else
                topdicis(topi,ncoltopdicis+1) = str2double(char(numcomps(topi)));
            end
            
            topdicis(topi,ncoltopdicis+2) = 0;
            for ri = max(r-1,1):min(r+1,size(tncmelodic,1))
                for ci = max(c-1,1):min(c+1,size(tncmelodic,2))
                    for pi = max(p-1,1):min(p+1,size(tncmelodic,3))
                        if template(ri,ci,pi)==1
                            topdicis(topi,ncoltopdicis+2) = 1;
                        end
                    end
                end
            end
        end
        topdicis = topdicis(topdicis(:,ncoltopdicis+2)==1,:);
        [~,I] = max(topdicis(:,2));
        if topdicis(I,ncoltopdicis+1)==0
            outdir = [fmriprepdir subid '/restmelodicsmooth_auto'];
        else
            outdir = [fmriprepdir subid '/restmelodicsmooth_' num2str(topdicis(I,ncoltopdicis+1))];
        end
        nummelodic = spm_read_vols(spm_vol([outdir '/melodic_IC.nii']));
        dici_templatefinal=zeros(size(nummelodic,4),1);
        for compi = 1:size(nummelodic,4)
            zind = 1;
            for zi = thresholds(1:length(thresholds)-1)
                component=squeeze(nummelodic(:,:,:,compi));
                component = component>zi;
                intersect_templatefinal = component & template;
                hitrate_templatefinal = sum(intersect_templatefinal(:))/sum(template(:));
                if hitrate_templatefinal==1
                    hitrate_templatefinal = (sum(template(:))-1)/sum(template(:));
                elseif hitrate_templatefinal==0
                    hitrate_templatefinal = 1/sum(template(:));
                end
                falsealarm_templatefinal = component & ~template & brainmask;
                outsidetemplate_insidebrainmask = ~template & brainmask;
                falsealarmrate_templatefinal = (sum(falsealarm_templatefinal(:))/sum(outsidetemplate_insidebrainmask(:)));
                if falsealarmrate_templatefinal==1
                    falsealarmrate_templatefinal = (sum(outsidetemplate_insidebrainmask(:))-1)/sum(outsidetemplate_insidebrainmask(:));
                elseif falsealarmrate_templatefinal==0
                    falsealarmrate_templatefinal = 1/sum(outsidetemplate_insidebrainmask(:));
                end
                hr_far_templatefinal = norminv(hitrate_templatefinal, 0, 1)-norminv(falsealarmrate_templatefinal, 0, 1);
                dici_templatefinal(compi, zind) = hr_far_templatefinal;
                zind = zind+1;
            end
            
        end
        
        meandici = transpose(mean(dici_templatefinal,2));
        [~, maxmeandiciind] = max(meandici);
        [sortvals, ~] = sort(meandici,'descend');
        sortmeandiffs = abs(diff(sortvals));
        maxdiff = max(sortmeandiffs);
        [~, data2cluster, cluster2dataCell] = MeanShiftCluster(meandici,maxdiff);
        
        maxclust = cell2mat(cluster2dataCell(data2cluster(maxmeandiciind)));
        finalnumcomps = 0;
        finalnetwork = zeros(size(brainmask));
        for clustcompi = 1:length(maxclust)
            onevoxpeakcheck = 0;
            clustcomp = squeeze(nummelodic(:,:,:,maxclust(clustcompi)));
            [~, idx] = max(clustcomp(:));
            [r,c,p] = ind2sub(size(clustcomp), idx);
            for ri = max(r-1,1):min(r+1,size(tncmelodic,1))
                for ci = max(c-1,1):min(c+1,size(tncmelodic,2))
                    for pi = max(p-1,1):min(p+1,size(tncmelodic,3))
                        if template(ri,ci,pi)==1
                            onevoxpeakcheck = 1;
                        end
                    end
                end
            end
            if onevoxpeakcheck==1
                finalnumcomps = finalnumcomps+1;
                
                overlapclust = finalnetwork & clustcomp;
                if sum(overlapclust(:)) == 0
                    finalnetwork = finalnetwork+clustcomp;
                else
                    finalnetworkoverlap = finalnetwork;
                    finalnetworkoverlap(~overlapclust) = 0;
                    clustcompoverlap = clustcomp;
                    clustcompoverlap(~overlapclust) = 0;
                    clustcomp(overlapclust) = 0;
                    finalnetwork(overlapclust) = 0;
                    overlapclust = double(overlapclust);
                    for i = 1:size(overlapclust,1)
                        for j = 1:size(overlapclust,2)
                            for k = 1:size(overlapclust,3)
                                if overlapclust(i,j,k)==1
                                    overlapclust(i,j,k) = max(finalnetworkoverlap(i,j,k), clustcompoverlap(i,j,k));
                                end
                            end
                        end
                    end
                    finalnetwork = finalnetwork+clustcomp;
                    finalnetwork = finalnetwork+overlapclust;
                end
                
            end
        end
        
        if ~isfile([outdir '/mean.nii'])
            gunzip([outdir '/mean.nii.gz' ]);
        end
        compheader = spm_vol([outdir '/mean.nii']);
        compheader.fname = [outdir '/' subid '_wIMM.nii'];
        compheader.private.dat.fname = compheader.fname;
        spm_write_vol(compheader, finalnetwork);

    end
    toc;
end
end
