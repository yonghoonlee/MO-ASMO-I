%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Generates a randomized samples of Latin Hypercube design.
%===============================================================================
function x = samplingLHS(number,dimension)
    % Random permutation to create Latin Hypercube
    x = zeros(number,dimension);
    for i = 1:dimension
        x(:,i) = randperm(number)';
    end
    %---------------------------------------------------------------------------
    % Randomize within Latin Hypercube voxel.
    x = (x-0.5)/number;
    rn = rand(size(x));
    rn = (rn-0.5)/number;
    x = x + rn;
end
%===============================================================================
