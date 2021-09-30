function mMCM_creation(fmriprepdir, templatepath)
clear;
subs = dir(fmriprepdir);
subs = {subs.name};
cd(fmriprepdir);
for subi=1:length(subs)
    subid = char(subs(subi));
    disp(subid);
    tic;
    maskICAfile = [subid '/restmelodicsmoothmask_auto/melodic_IC.nii'];
    if isfile([maskICAfile '.gz'])
        
        if ~isfile(maskICAfile)
            gunzip ([maskICAfile '.gz']);
        end
        
        template = spm_read_vols(spm_vol(templatepath));
        template = template>0;
        
        
        meanICAfile = [subid '/restmelodicsmoothmask_auto/mean.nii'];
        if ~isfile(meanICAfile)
            gunzip ([meanICAfile '.gz']);
        end
        
        maskICAheader = spm_vol(maskICAfile);
        maskICA = spm_read_vols(maskICAheader);
        meanICAheader = spm_vol(meanICAfile);
        meanICA = spm_read_vols(spm_vol(meanICAheader));
        
        finalnetwork = zeros(size(meanICA));
        numcomps = size(maskICA, 4);
        for numi = 1:numcomps
            component = squeeze(maskICA(:,:,:,numi));
            component(component<0) = 0;
            
            [~, idx] = max(component(:));
            [r,c,p] = ind2sub(size(component), idx);
            
            if template(r,c,p)==1
                overlapclust = finalnetwork & component;
                if sum(overlapclust(:)) == 0
                    finalnetwork = finalnetwork+component;
                else
                    finalnetworkoverlap = finalnetwork;
                    finalnetworkoverlap(~overlapclust) = 0;
                    componentoverlap = component;
                    componentoverlap(~overlapclust) = 0;
                    component(overlapclust) = 0;
                    finalnetwork(overlapclust) = 0;
                    overlapclust = double(overlapclust);
                    for i = 1:size(overlapclust,1)
                        for j = 1:size(overlapclust,2)
                            for k = 1:size(overlapclust,3)
                                if overlapclust(i,j,k)==1
                                    overlapclust(i,j,k) = max(finalnetworkoverlap(i,j,k), componentoverlap(i,j,k));
                                end
                            end
                        end
                    end
                    finalnetwork = finalnetwork+component;
                    finalnetwork = finalnetwork+overlapclust;
                end
            end
        end
        meanICAheader.fname = [subid '/restmelodicsmoothmask_auto/' subid '_mMCM.nii'];
        meanICAheader.private.dat.fname = meanICAheader.fname;
        spm_write_vol(meanICAheader, finalnetwork);
        toc;
    end
end
end
