%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% SETUP_PATH FUNCTION
%==========================================================================
% This function adds all required directories into the Matlab path list.
% Input: relative path of the problem to be solved
% Output: relative path of the problem to be solved
%==========================================================================

function problem_path = setup_path(problem_path)
    fprintf('%s','Setup path...');

    pathlist = string({...
        'MOASMO';
        'Solver';
        'Utility';
        'Library/hcparula';
        'Library/export_fig';
        problem_path ...
    });

    % Begin
    currentpath = pwd;
    if ispc
        dirsep = '\';
    elseif isunix
        dirsep = '/';
    end

    pathlist = replace(pathlist, {'/','\'}, dirsep);
    for i = 1:length(pathlist)
        path(path,strcat(currentpath, dirsep, pathlist{i}));
    end

    fprintf('%s\n','done');
end