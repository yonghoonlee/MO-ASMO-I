% This function obtains number of available parallel pool.
% Usage: npool = HFF_GETPOOLSIZE()
% Input:
% Output: npool
%   npool: Number of active parallel pool (0, if no active parallel pool)

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function npool = hff_getpoolsize()
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj)
        npool = 0;
    else
        npool = poolobj.NumWorkers;
    end
end