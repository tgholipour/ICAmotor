function wMCM_component_generation(fslpath, fmriprepdir, numcomps, smoothkernel)
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};
subs = subs(3:length(subs));
for subi = 1:length(subs)
    tic;
    disp(char(subs(subi)))
    inputfilename = dir([fmriprepdir '/' char(subs(subi)) '/func/' char(subs(subi)) '_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    if ~isempty(inputfilename)            
        inputfilename = [inputfilename.folder '/' inputfilename.name];
        smoothfilename = [fmriprepdir '/' char(subs(subi)) '/func/' char(subs(subi)) '_task-rest_space-T1w_desc-preproc_bold_smooth.nii.gz'];
        if ~isfile(smoothfilename)
            smoothstring = ['export FSLDIR=' fslpath '; . ${FSLDIR}/etc/fslconf/fsl.sh; ' fslpath '/bin/fslmaths ' inputfilename ' -kernel gauss ' num2str(smoothkernel) ' -fmean ' smoothfilename];
            system(smoothstring);
        end
        for numi = numcomps
            if strcmp(char(numi), 'auto')
                outdir = [fmriprepdir '/' char(subs(subi)) '/restmelodicsmooth_' char(numi)];
                numICAstring = ['export FSLDIR=' fslpath '; . ${FSLDIR}/etc/fslconf/fsl.sh; ' fslpath '/bin/melodic -i ' smoothfilename ' -o ' outdir];
            else
                outdir = [fmriprepdir '/' char(subs(subi)) '/restmelodicsmooth_' char(numi)];
                numICAstring = ['export FSLDIR=' fslpath '; . ${FSLDIR}/etc/fslconf/fsl.sh; ' fslpath '/bin/melodic -i ' smoothfilename ' --dim=' char(numi) ' -o ' outdir];
            end
            if ~exist([outdir '/melodic_IC.nii.gz'], 'file')
                system(numICAstring);
            end
        end
    end
    toc;
end
end
