%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Separates feasible and infeasible sample points by high fidelity results
%===============================================================================
function [sampleFea, sampleInfea, objHFFFea, objHFFInfea, ...
          additionFea, additionInfea] = hffSeparateNaN( ...
        sample, objHFF, prob, varargin)
    if (prob.control.verbose > 0)
    	fprintf('%s','Separating infeasible results...');
    end
    %---------------------------------------------------------------------------
    % if additional data exists
    if (nargin == 4)
    	addition = varargin{1};
    else
    	addition = sample;
    end
    %---------------------------------------------------------------------------
    % indices of NaN containing rows (idxInfea) / feasible rows (idxFea)
    idxInfea = find(any(isnan([sample,objHFF,addition]),2));
    idxFea = find(~any(isnan([sample,objHFF,addition]),2));
    %---------------------------------------------------------------------------
    % results
    sampleFea = sample(idxFea,:);
    sampleInfea = sample(idxInfea,:);
    objHFFFea = objHFF(idxFea,:);
    objHFFInfea = objHFF(idxInfea,:);
    additionFea = addition(idxFea,:);
    additionInfea = addition(idxInfea,:);
end
%===============================================================================
