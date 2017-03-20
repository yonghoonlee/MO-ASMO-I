% This function generates a randomized samples of Latin Hypercube design.
% Usage: x = LHS(number,dimension)
% Input: number, dimension
% Output: x
%   number: Number of samples to be generated
%   dimension: Number of design variables
%   x: Generated sample points

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function x = LHS(number,dimension)
    % Random permutation to create Latin Hypercube
    x = zeros(number,dimension);
    for i = 1:dimension
        x(:,i) = randperm(number)';
    end
    % Randomize within Latin Hypercube voxel.
    x = (x-0.5)/number;
    rn = rand(size(x));
    rn = (rn-0.5)/number;
    x = x + rn;
end