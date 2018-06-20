%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Error evaluation
%===============================================================================
function [errorVec, errorAvg] = errorEvaluation(hffF, surF, prob, k, eH, surr)
    [nFea, mFea] = size(hffF);
    
    if (k == 1) && (nFea == 0)
        errorVec = 1;
        errorAvg = 1;
    elseif (k > 1) && (nFea == 0)
        maxError = 0;
        for idx = 1:(k-1)
            data = cell2mat(eH(idx));
            prevError = max(data);
            maxError = max(maxError, prevError);
        end
        errorVec = maxError;
        errorAvg = maxError;
    else
        if prob.bound.adaptive
            flb = surr.scale.flb;
            fub = surr.scale.fub;
        else
            flb = prob.bound.flb;
            fub = prob.bound.fub;
        end
        flbm = repmat(reshape(flb,1,mFea),nFea,1);
        fubm = repmat(reshape(fub,1,mFea),nFea,1);
        normHffFFea = (hffF - flbm)./(fubm - flbm);
        normSurFFea = (surF - flbm)./(fubm - flbm);
        errorVec = sqrt(sum((normHffFFea - normSurFFea).^2,2));
        errorAvg = sum(errorVec)/size(errorVec,1);
    end
end
%===============================================================================
