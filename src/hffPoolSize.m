%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Returns parallel pool size if exists; otherwise, returns 0
%===============================================================================
function npool = hffPoolSize()
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj)
        npool = 0;
    else
        npool = poolobj.NumWorkers;
    end
end
%===============================================================================
