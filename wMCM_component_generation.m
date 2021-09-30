function wMCM_component_generation(fslpath, fmriprepdir, numcomps, smoothkernel)
cd(fmriprepdir)
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};
for subi = 1:length(subs)
    tic;
    disp(char(subs(subi)))
    inputfilename = dir([char(subs(subi)) '/func/' char(subs(subi)) '_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    if ~isempty(inputfilename)            
        inputfilename = [inputfilename.folder '/' inputfilename.name];
        smoothfilename = [char(subs(subi)) '/func/' char(subs(subi)) '_task-rest_space-T1w_desc-preproc_bold_smooth.nii.gz'];
        if ~isfile(smoothfilename)
            smoothstring = [fslpath '/fslmaths ' inputfilename ' -kernel gauss ' smoothkernel ' -fmean ' smoothfilename];
            system(smoothstring);
        end
        for numi = numcomps
            if strcmp(char(numi), 'auto')
                outdir = [char(subs(subi)) '/restmelodicsmooth_' char(numi)];
                numICAstring = [fslpath '/melodic -i ' smoothfilename ' -o ' outdir];
            else
                outdir = [char(subs(subi)) '/restmelodicsmooth_' char(numi)];
                numICAstring = [fslpath '/melodic -i ' smoothfilename ' --dim=' char(numi) ' -o ' outdir];
            end
            if ~exist([outdir '/melodic_IC.nii.gz'], 'file')
                system(numICAstring);
            end
        end
    end
    toc;
end
end
