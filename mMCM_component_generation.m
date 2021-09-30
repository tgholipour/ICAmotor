function component_generation_mMCM(fslpath, afnipath, antspath, fmriprepdir, mMCM_maskpath, smoothkernel)

clear;
cd(fmriprepdir)
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};

for subi = 1:length(subs)
    tic;
    disp(char(subs(subi)))
    subid = char(subs(subi));
    inputfilename = dir([char(subs(subi)) '/func/' char(subs(subi)) '_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    if ~isempty(inputfilename)
        inputfilename = [inputfilename.folder '/' inputfilename.name];
        smoothfilename = [char(subs(subi)) '/func/' char(subs(subi)) '_task-rest_space-T1w_desc-preproc_bold_smooth.nii.gz'];
        if ~isfile(smoothfilename)
            smoothstring = [fslpath '/fslmaths ' inputfilename ' -kernel gauss ' smoothkernel ' -fmean ' smoothfilename];
            system(smoothstring);
        end
        native_mMCMmask_file = [subid '/template/mMCMmask_native.nii'];
        if isfile(native_mMCMmask_file)
            delete(native_mMCMmask_file);
        end
        native_transform_mat = [char(subs(subi)) '/anat/' subid '_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5'];
        native_anat = [char(subs(subi)) '/anat/' subid '_desc-preproc_T1w.nii.gz'];
        nativize_mMCMmask_string = [antspath ' -i ' mMCM_maskpath ' -r ' native_anat ' -t ' native_transform_mat ' -o ' native_mMCMmask_file];
        system(nativize_mMCMmask_string);
        deob_mMCMmask_file = [subid '/template/mMCMmask_native_deob.nii'];
        if isfile(deob_mMCMmask_file)
            delete(deob_mMCMmask_file);
        end
        deoblique_mMCMmask_string = [afnipath '/3dWarp -deoblique -NN -gridset ' smoothfilename ' -prefix ' deob_mMCMmask_file ' ' native_mMCMmask_file];
        system(deoblique_mMCMmask_string);
        if ~isempty(smoothfilename)
            outdir = [char(subs(subi)) '/restmelodicsmoothmask_auto/'];
            mICAstring = [fslpath '/melodic -i ' smoothfilename ' -o ' outdir ' --mask=' deob_mMCMmask_file];
            if ~exist([outdir '/melodic_IC.nii.gz'], 'file')
                system(mICAstring);
            end
        end
        toc;
    end
end