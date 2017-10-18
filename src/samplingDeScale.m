%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% De-scale training points from [0,1] to [lb,ub].
%===============================================================================
function x = samplingDeScale(xs,lb,ub)
	lb = reshape(lb,1,numel(lb));
	ub = reshape(ub,1,numel(ub));
    repeat = size(xs,1);
    lb = repmat(lb,repeat,1);
    ub = repmat(ub,repeat,1);
    x = lb + (ub - lb).*xs;
end
%===============================================================================
