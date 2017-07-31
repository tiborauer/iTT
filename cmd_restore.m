function o = cmd_restore
% restore changes made by spm_config_TT
% FORMAT o = cmd_restore;
%
% This function is usually called from TT Configuration GUI, but in case of
% need it is possible to call from command line.
%_______________________________________________________________________
% Copyright (C) 2009 Biomedizinische NMR Forschungs GmbH

% Tibor Auer
% $Id: spm_config_TT.m 2009-06-28 $ 

p = which('restore_getSPM.m');
o = 0;
if ~isempty(p)
    restore_getSPM;
    delete(p);
    o = 1;
end