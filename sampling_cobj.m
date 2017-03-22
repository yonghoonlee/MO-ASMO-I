% This function computes obj fn to make samples satisfying constraints.
% Usage: f = SAMPLING_COBJ(x,x0,wei1,wei2,prevPoints)
% Input: x,x0,wei1,wei2,prevPoints
% Output: f
%   x: relocated position (design variable) by optimization algorithm
%   x0: original sample location (may not satisfy constraint)
%   wei1: weight 1
%   wei2: weight 2
%   prevPoints: existing points
%   f: objective function

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function f = sampling_cobj(x,x0,wei1,wei2,prevPoints)
    xm = repmat(x,size(prevPoints,1),1);
    f1 = wei1*sum((x - x0).^2);
    if (size(prevPoints,1) == 0)
        f2 = 0;
    else
        f2 = wei2*sum(-log(max(sum((xm - prevPoints).^2,2),1e-12)))/size(prevPoints,1);
    end
    f = f1 + f2;
end