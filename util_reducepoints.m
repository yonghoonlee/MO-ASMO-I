% This function reduces data point by k-mean clustering or random sample.
% Usage: xf = UTIL_REDUCEPOINTS(xin,n,method)
% Input: xin,n,method
% Output: xf
%   xin: Input data set
%   n: Number of desired data points
%   method: Sampling methods 'KMC' or 'RND'
%   xf: Sampled data set

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function xf = util_reducepoints(xin,n,method)
    if (size(xin,1) <= n)
        xf = xin;
    else
        switch lower(method)
            case 'kmc'
                [~,xf] = kmeans(xin,n);
            case 'rnd'
                [xf,~] = datasample(xin,n,'Replace',false);
        end
    end
end