function job = spm_config_iTT
% automatic configuration of the TT toolbox
% FORMAT job = spm_config_TT;
%
% During evry startup of SPM, this program checks whether the modified
% spm_getSPM.m file is to be used. If this check fails, the program renames
% the located spm_getSPM.m files into spm_getSPM.m.org.
%______________________________________________________________________
% Copyright (C) 2011 Biomedizinische NMR Forschungs GmbH

% Tibor Auer
% $Id: spm_config_TT.m 2011-07-06 $ 

%-add toolbox directory into the MATLAB path
%--------------------------------------------------------------------------
addpath(fullfile(spm('dir'),'toolbox','iTT'));

%-checking MATLAB path for existing spm_getSPM.m files
%-renaming them if they are located upper in the path than the modified one
%-saving changes into restore_getSPM.m if there was any
%--------------------------------------------------------------------------
iTT = IniFile(fullfile(spm('dir'),'toolbox','iTT','config.ini'));
if iTT.config.check_getSPM
    fprintf('\nChecking path for spm_getSPM.m\n');
    p = fileparts(mfilename('fullpath'));
    isMod = false;
    fid = fopen(fullfile(p, 'restore_getSPM.m'),'a');
    w = path_check('spm_getSPM');
    for i = 1:numel(w)
        if isempty(strfind(help(fullfile(w{i},'spm_getSPM.m')),'Tibor Auer'))
            movefile(fullfile(w{i},'spm_getSPM.m'), fullfile(w{i},'spm_getSPM.m.org'));
            copyfile(fullfile(p,'spm_getSPM.m'), fullfile(w{i},'spm_getSPM.m'));
            isMod = true;
            fprintf(fid,['movefile(fullfile(''%s'', ''spm_getSPM.m.org''), fullfile(''%s'', ''spm_getSPM.m''));\n'],w{i},w{i});
        end
    end
    fclose(fid);
    
    %-informing user
    %--------------------------------------------------------------------------
    fprintf('Initialization of iterative Two-Threshold (iTT) toolbox is completed. Settings were loaded from config.ini.\n');
    fprintf('Ignore ''Loading of toolbox %s\\spm_config_iTT.m failed.'' message!\n',p);
    if isMod
        fprintf('Previous spm_getSPM.m files were renamed! Use restore_getSPM to roll back changes!\n');
    end
end
end

function w = path_check(fn)
% looks for a file fn in the MATLAB path and lists every location of it
% FORMAT w = path_check(fn);
%
% Required input:
%   fn  - file of interest
%
% Provided output:
%   w   - cell array of full paths of every loctation of the file fn
%

%-reading MATLAB path
%--------------------------------------------------------------------------
p = path;

%-turning pathstring into a cell array of paths
%--------------------------------------------------------------------------
ind = [0 find(path == pathsep)];
for i = 1:numel(ind)-1
    P{i,1} = p(ind(i)+1:ind(i+1)-1);
end

%-putting path into the output cell array, if it is a valid location of
%file fn
%--------------------------------------------------------------------------
w = {};
for i = 1:numel(P)
    if exist(fullfile(P{i},fn),'file')
        w{end+1} = P{i};
    end
end
end