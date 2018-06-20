%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Approximate ratio of non-dominated hypervolume to hyper-cube volume between
% utopia and anti-utopia using Monte Carlo sampling method
% Original method discussed in:
% Zitzler, E. and Thiele, L.,
% "Multiobjective evolutionary algorithms: a comparative study and the strength
% Pareto approach," IEEE Transactions on Evolutionary Computation, 3(4), 1999,
% pp. 257?271. doi: 10.1109/4235.797969
% Further discussion and summary appears in:
% (1) van Veldhuizen, David A.,
% "Multiobjective evolutionary algorithms: classifications, analyses, and new
% innovations,"	Doctoral Dissertation, Air Force Institute of Technology,
% Wright-Patterson Air Force Base, OH, USA, 1999. isbn: 0-599-28316-5
% (2) Janssens, Gerrit K. and Pangilinan, José M.,
% "Multiple Criteria Performance Analysis of Non-dominated Sets Obtained by
% Multi-objective Evolutionary Algorithms for Optimisation,"
% In Artificial Intelligence Applications and Innovations, Springer,
% Berlin, 2010, pp. 94-103. doi: 10.1007/978-3-642-16239-8_15
%===============================================================================
function [hv,rhv] = approxNDHypervolume(f,n)
    [nf,mf] = size(f);
    uf = min(f);
    auf = max(f);
    samples = uf + (auf - uf).*rand(n,mf);
    dominated = 0;
    for ii = 1:nf
        idx = (sum((f(ii,:) <= samples),2) == mf);
        dominated = dominated + sum(idx);
        samples(idx,:) = [];
    end
    hvtotal = prod(auf - uf);
    rhv = 1 - dominated/n;
    hv = rhv*hvtotal;
end