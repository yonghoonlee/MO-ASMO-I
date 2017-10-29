%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Non-dominated sorting algorithm
%===============================================================================
function [xsort,fsort,ndidx] = ndSort(xin,fin)
% ndSort code adopted by non_domination_sort_mod function of NSGA-II
% implemented by Aravind Seshadri in 2009 (with modification by Yong Hoon Lee)
% from https://www.mathworks.com/matlabcentral/fileexchange/10429
% See the license/NSGA-II.License for detail.
    x = xin;
    f = fin;
    [nx,mx] = size(x);
    [nf,mf] = size(f);
    if (nx ~= nf)
        error('dimension not match');
    end
    %---------------------------------------------------------------------------
    X = [x, f, zeros(nx,1)];
    frtID = 1;
    frtStr(frtID).frt = [];
    idvl = [];
    %---------------------------------------------------------------------------
    for idx1 = 1:nx
        idvl(idx1).n = 0;
        idvl(idx1).p = [];
        for idx2 = 1:nx
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for idx3 = 1:mf
                if (X(idx1,mx+idx3) < X(idx2,mx+idx3))
                    dom_less = dom_less + 1;
                elseif (X(idx1,mx+idx3) == X(idx2,mx+idx3))
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if ((dom_less == 0) && (dom_equal ~= mf))
                idvl(idx1).n = idvl(idx1).n + 1;
            elseif ((dom_more == 0) && (dom_equal ~= mf))
                idvl(idx1).p = [idvl(idx1).p, idx2];
            end
        end
        if (idvl(idx1).n == 0)
            X(idx1,mx+mf+1) = 1;
            frtStr(frtID).frt = [frtStr(frtID).frt, idx1];
        end
    end
    %---------------------------------------------------------------------------
    while ~isempty(frtStr(frtID).frt)
        Q = [];
        for idx1 = 1:length(frtStr(frtID).frt)
            if ~isempty(idvl(frtStr(frtID).frt(idx1)).p)
                for idx2 = 1:length(idvl(frtStr(frtID).frt(idx1)).p)
                    idvl(idvl(frtStr(frtID).frt(idx1)).p(idx2)).n ...
                        = idvl(idvl(frtStr(frtID).frt(idx1)).p(idx2)).n - 1;
                    if idvl(idvl(frtStr(frtID).frt(idx1)).p(idx2)).n == 0
                        X(idvl(frtStr(frtID).frt(idx1)).p(idx2),mf+mx+1) ...
                            = frtID + 1;
                        Q = [Q, idvl(frtStr(frtID).frt(idx1)).p(idx2)];
                    end
                end
            end
        end
        frtID = frtID + 1;
        frtStr(frtID).frt = Q;
    end
    %---------------------------------------------------------------------------
    [~,idxv1] = sort(X(:,mf+mx+1));
    for idx1 = 1:length(idxv1)
        Xsrt(idx1,:) = X(idxv1(idx1),:);
    end
    %---------------------------------------------------------------------------
    xsort = Xsrt(:,1:mx);
    fsort = Xsrt(:,(mx+1):(mx+mf));
    ndidx = Xsrt(:,(mx+mf+1));
end
%===============================================================================
