# iTT
iTT toolbox for SPM

iTT toolbox for SPM implements the itertive Two-Threshold (iTT) approach
to determine height thresholds in SPM.

Reference:
Auer T, Schweizer R, Frahm J. An iterative two-threshold analysis 
for single-subject functional MRI of the human brain.
[DOI: 10.1007/s00330-011-2184-5](http://dx.doi.org/10.1007/s00330-011-2184-5)

This file contains information about Install Configuration Usage of the
toolbox.
___________________________________________________________________________
Copyright (C) 2011 Biomedizinische NMR Forschungs GmbH

Tibor Auer
$Id: README 2011-07-06 $ 


System requirements
---------------------------------------------------------------------------
MATLAB (The MathWorks) 
MATLAB 7.6 or later is required.

Curve Fitting Toolbox - required for Gaussian fitting.
Optimization Toolbox - required for iterative optimization.

SPM

Note: You need write permission - during installation and running - in the
'spm' directory (SPMDIR) and in its subdirectories where other spm_getSPM.m
files might be stored (e.g. SPMDIR\toolbox\vbm5).


Testing environment
---------------------------------------------------------------------------
Ubuntu 11.04. 64-bit, Windows 7 64-bit
MATLAB 7.10 (R2010a)
SPM8 v4290


Install
---------------------------------------------------------------------------
Simply copy the iTT folder into the SPMDIR\toolbox dircetory.
spm_config_TT.m will arrange everything during startup of SPM. 

Note: Previous spm_getSPM.m files (e.g. in SPMDIR) will be renamed
to spm_getSPM.m.org. restore_getSPM.m file will be created in the TT folder
to roll back these changes, using cmd_restore.m.


Configure
---------------------------------------------------------------------------
By default, the lower threshold is set to p=0.05, itertive optimization is 
employed (iTT) to determined the lower and upper cuts as well as the upper 
threshold, and results of the fitting will be displayed. Without the itertive
optimization (TT) the central portion of the histogram determined by the
20% and 80% maximum height will be fitted and the upper threshold of p=0.0001 
will be used.
Launch iTT toolbox in the SPM's toolboxes to modify these settings.
This Configuration Toolbox can be also used to disable automatic check for
spm_get_SPM.m files during startup (not recommended), or to roll back changes
during Install (see above).
After accepting changes, they will have an immediate effect and changes 
will be also written into the configuration file (config.ini) to store them 
for the next session. 


Use
---------------------------------------------------------------------------
After sucessfull install, iTT and TT can be used as easily as other approaches 
(FWE or FDR). During displaying Result, simlpy select iTT or TT when asked for
"p-value adjustment". For TT, give probaility threshold for the lower threshold
(default: 0.05).

Note: For calculation and application of the thresholds are performed in
one step, iTT and TT may need more time than other approaches. 


Uninstall
---------------------------------------------------------------------------
Note: Run cmd_restore.m file in to roll back renaming previous 
spm_getSPM.m files (e.g. in SPMDIR).

Simply remove the TT folder from the SPMDIR\toolbox directory.


Versions
---------------------------------------------------------------------------
1.1 - cell_index.m has been added
1.0	- Initial Release
