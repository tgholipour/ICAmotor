function component_generation(fslpath, fmriprepdir, numcomps)
if nargin==0
    fslpath = '';
    fmriprepdir = 'Data/';
    numcomps = {'auto', '20', '30', '40', '50', '60'};
end
cd(fmriprepdir)
inputsubs=dir(fmriprepdir);
subs = {inputsubs.name};
for subi = 1:length(subs)
    disp(char(subs(subi)))
    inputfilename = dir([[char(subs(subi)) '/rest/'], '*_task-rest*_space-T1w_desc-preproc_bold.nii.gz']);
    if ~isempty(inputfilename)
        for numi = numcomps
            outdir = [char(subs(subi)) '/restmelodic_' char(numi)];
            numICAstring = [fslpath '/bin/melodic -i ' char(subs(subi)) '/rest/' inputfilename ' --dim=' char(numi) ' -o ' outdir];
            if ~exist([outdir '/melodic_IC.nii.gz'], 'file')
                system(numICAstring);
            end
        end
    end
end
end