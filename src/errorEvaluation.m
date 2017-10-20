function [errorVec, errorAvg] = errorEvaluation(hffFFea, surFFea, prob, k, eH)
    [nFea, mFea] = size(hffFFea);
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
        flbm = repmat(reshape(prob.bound.flb,1,mFea),nFea,1);
        fubm = repmat(reshape(prob.bound.fub,1,mFea),nFea,1);
        normHffFFea = (hffFFea - flbm)./(fubm - flbm);
        normSurFFea = (surFFea - flbm)./(fubm - flbm);
        errorVec = sqrt(sum((normHffFFea - normSurFFea).^2,2));
        errorAvg = sum(errorVec)/size(errorVec,1);
    end
end
