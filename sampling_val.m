% This function generates validation samples from the predicted Pareto set.
% Usage: [xV,fV] = SAMPLING_VAL(xP,fP,problem)
% Input: xP,fP,problem
% Output: xV,fV
%   xP: Predicted Pareto set in design space
%   fP: Predicted Pareto set in objective function space
%   problem: Problem definition structure
%   xV: Validation sample points in design space
%   fV: Validation sample points in objective function space (predicted)

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [xV,fV] = sampling_val(xP,fP,problem)
    % Size of predicted Pareto set
    [nx,mx] = size(xP);
    [nf,mf] = size(fP);
    number = problem.sampling.valnumber;
    xV = [];
    fV = [];
    
    % Scale
    lbx = ones(nx,1)*min(xP)-1e-10;
    ubx = ones(nx,1)*max(xP)+1e-10;
    x = (xP - lbx)./(ubx - lbx);
    lbf = ones(nf,1)*min(fP)-1e-10;
    ubf = ones(nf,1)*max(fP)+1e-10;
    f = (fP - lbf)./(ubf - lbf);
    
    % Choose extreme points from the predicted Pareto set
    x_minmaxlist = [];
    for i = 1:mf
        [~,ixmin] = min(f(:,i));
        [~,ixmax] = max(f(:,i));
        x_minmaxlist = [x_minmaxlist; ixmin; ixmax];
    end
    for i = 1:size(x_minmaxlist,1)
        xV = [xV; x(x_minmaxlist(i),:)];
        fV = [fV; f(x_minmaxlist(i),:)];
    end
    % Eliminate duplicates
    try clear tmp; end
    tmp = unique([xV,fV],'rows');
    xV = tmp(:,1:mx);
    fV = tmp(:,(mx+1):end);
    try clear tmp; end
    
    % Adjust numbers
    number = number - size(xV,1);
    fnumber = round(sqrt(mf)/(sqrt(mx)+sqrt(mf))*number);
    number = number - fnumber;
    
    % Remaining x is stored separately
    xrem = x;
    frem = f;
    
    % Searching validation points in the design space
    while (size(xrem,1)>number)
        % Compute distance between each point
        clear dist;
        [nxrem,~] = size(xrem);
        dist = zeros(nxrem,nxrem);    % distance between each point
        for i = 1:nxrem
            for j = 1:nxrem
                dist(i,j) = sqrt(sum((xrem(i,:) - xrem(j,:)).^2));
            end
        end
        % Sort distance
        [~,distmk] = sort(dist,2);
        % Find tentative shortest distance points
        sdlist = zeros(nxrem-1,2);
        for i = 1:nxrem-1
            k = 2;
            fnd = false;
            while ((k <= nxrem) && (fnd == false))
                if (distmk(i,k) > i)
                    sdlist(i,1) = i;
                    sdlist(i,2) = distmk(i,k);
                    sdlist(i,3) = dist(i,round(sdlist(i,2)));
                    fnd = true;
                end
                k = k + 1;
            end
        end
        % Sort in the shortest distance point pairs
        [~,distst] = sort(sdlist(:,3),1);
        % If too much points left for deletion remove 10 points at a time
        rmlist = [];
        if (nxrem - number) > 25
            Nrm = 10;
        else
            Nrm = 1;
        end
        for i = 1:Nrm
            idx = distst(i);
            rmlist = [rmlist; round(sdlist(idx,round(rand)+1))]; % Choose one of two randomly
        end
        rmlist = unique(rmlist);
        % Removing points
        indx = true(nxrem,1);
        indx(rmlist) = false;
        xrem = xrem(indx,:);
        frem = frem(indx,:);
    end
    
    xV = [xV; xrem];
    fV = [fV; frem];
    clear xrem;
    clear frem;
    
    xrem = x;
    frem = f;
    
    % Searching validation points in the objective function space
    for k = 1:fnumber
        try clear dist; end
        try clear idx; end
        [nfrem,~] = size(frem);
        [nfV,~] = size(fV);
        dist = zeros(nfrem,nfV);    % distance between each point
        for i = 1:nfrem
            for j = 1:nfV
                dist(i,j) = sqrt(sum((frem(i,:) - fV(j,:)).^2));
            end
        end
        mindist = min(dist,[],2);
        [~,idx] = max(mindist);
        xV = [xV; x(idx,:)];
        fV = [fV; f(idx,:)];
    end
    
    try clear tmp; end
    tmp = unique([xV,fV],'rows');
    xV = tmp(:,1:mx);
    fV = tmp(:,(mx+1):end);
    nxV = size(xV,1);
    try clear tmp; end
    
    % Descale
    lbx = lbx(1:nxV,:);
    ubx = ubx(1:nxV,:);
    xV = lbx + (ubx - lbx).*xV;
    lbf = lbf(1:nxV,:);
    ubf = ubf(1:nxV,:);
    fV = lbf + (ubf - lbf).*fV;
    
end