% This fn generates samples for exploitation by Cartesian domain division.
% Usage: xf = SAMPLING_CDD(xP,number,problem)
% Input: xP,number,problem
% Output: xf
%   xP: predicted Pareto set
%   number: number of samples to be created
%   problem: problem definition structure
%   xf: sample points for exploitation of Pareto-vicinity region

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function xf = sampling_CDD(xP,number,problem)
    if (problem.control.verbose > 0)
        fprintf('Cartesian domain division method...');
    end
    
    
    minxP = min(xP);
    maxxP = max(xP);
    minD = minxP - 0.5*(maxxP - minxP);         % expand 50%
    maxD = maxxP + 0.5*(maxxP - minxP);         % expand 50%
    minD = max([minD; problem.xlb']);           % limit to lb
    maxD = min([maxD; problem.xub']);           % limit to ub
    sizeD = maxD - minD;
    
    i = 0;
    flg = true;
    while (flg)
        i = i + 1;
        
        % Generate list of all point combinations
        clear blist;
        clear clist;
        clear dlist;
        xlist = [];
        for j = 1:i+1
            blist(j,:) = [minD + (j-1)*sizeD/i];
        end
        for k = 1:size(blist,2)
            for j = 1:size(blist,1)
                clist{k}(1,j) = blist(j,k);
            end
        end
        dlist = [combvec(clist{:})]';
        % Check the distance from existing points
        for j = 1:size(dlist,1)
            flg = false;
            for k = 1:size(xP,1)
                distance = abs(dlist(j,:)-xP(k,:));
                if (all(distance < 0.05*sizeD)) % 5% threshold
                    flg = true;
                end
            end
            if (flg == true)
                xlist = [xlist; dlist(j,:)];
            end
        end
        if (size(xlist,1) >= number)
            flg = false;
        else
            flg = true;
        end
    end
    xf = util_reducepoints(xlist,n,'RND');
    % 1% distance randomize
    rx = (2*(rand(size(xf))-0.5)).*(ones(size(xf,1),1)*(0.01*sizeD));
    xf = xf + rx;
    
end
    