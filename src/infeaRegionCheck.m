%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Check if each point is located in the infeasible domain
% Reference: Malak Jr., R.J. and Paredis, C.J.J., "Using Support Vector Machines
% to Formalize the Valid Input Domain of Predictive Models in Systems Design
% Problems," Journal of Mechanical Design, 132, Oct 2010, pp.101001,
% doi: 10.1115/1.4002151
%===============================================================================
function feasibility = infeaRegionCheck(infeaRegion,xtest)
    feasibility = zeros(size(xtest,1),1);
    if ~isstruct(infeaRegion)
        return;
    end
    %---------------------------------------------------------------------------
    dopt = infeaRegion.dopt;
    p.C = infeaRegion.C;
    p.q = infeaRegion.q;
    p.epsilon = infeaRegion.epsilon;
    p.x = infeaRegion.x;
    nx = size(p.x,1);
    nxtest = size(xtest,1);
    for k = 1:nxtest
        z = xtest(k,:);
        Kxz = sum(dopt.*exp(-p.q*sum((p.x-repmat(z,nx,1)).^2,2)));
        Kxx = sum(dopt.*exp(-p.q*sum((p.x-repmat(p.x(1,:),nx,1)).^2,2)));
        feasibility(k,1) = 2*(Kxz - Kxx) + p.epsilon;
    end
end
%===============================================================================
