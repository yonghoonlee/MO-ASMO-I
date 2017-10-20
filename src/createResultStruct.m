%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Automatically generates the result struct used for initializing the problem
%===============================================================================
function result = createResultStruct(prob)
    result = struct('x',[],'fval',[],'exitflag',[],'output',[],...
        'population',[],'scores',[],'data',[],'time',[],'prob',prob);
    result.data = table;
    result.time = table;
end
%===============================================================================
