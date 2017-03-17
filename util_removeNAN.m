%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% UTIL_REMOVENAN FUNCTION
%==========================================================================
% This function removes infeasible (NaN) sample points
% Input: x, f, problem, xx
%   x: designs points
%   f: result of function evaluation
%   xx: additional data to be treated
% Output: [X, F, XX]
%==========================================================================

function [X,F,XX] = util_removeNAN(x,f,problem,varargin)
    if (problem.control.verbose > 0)
        fprintf('%s','Removing infeasible results...');
    end
    X = x;
    F = f;
    if (nargin == 4)
        XX = varargin{1};
    else
        XX = X;
    end
    for i = 1:size(X,2)
        clear idx;
        [idx,~] = find(~isnan(X(:,i)));
        X = X(idx,:);
        F = F(idx,:);
        XX = XX(idx,:);
    end
    for i = 1:size(F,2)
        clear idx;
        [idx,~] = find(~isnan(F(:,i)));
        X = X(idx,:);
        F = F(idx,:);
        XX = XX(idx,:);
    end
    for i = 1:size(XX,2)
        clear idx;
        [idx,~] = find(~isnan(XX(:,i)));
        X = X(idx,:);
        F = F(idx,:);
        XX = XX(idx,:);
    end
    n_removed = size(x,1) - size(X,1);
    if n_removed > 0
        if (problem.control.verbose > 0)
            fprintf('%s\n',strcat(num2str(n_removed),...
                '-points removed...done'));
        end
    else
        if (problem.control.verbose > 0)
            fprintf('%s\n','nothing removed...done');
        end
    end
end