% This function converts values in [0,1] space to [lb,ub] space.
% Usage: x = SAMPLING_DESCALING(xs,lb,ub)
% Input: xs,lb,ub
% Output: x
%   xs: Matrix of size n*m - n number of m dimension values in [0,1] space
%   lb: Lower bound
%   ub: Upper bound
%   x: Matrix of size n*m - n number of m dimension values in [lb,ub] space

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function x = sampling_descaling(xs,lb,ub)
    lb = lb';
    ub = ub';
    repeat = size(xs,1);
    lb = repmat(lb,repeat,1);
    ub = repmat(ub,repeat,1);
    x = lb + (ub - lb).*xs;
end