% This fn generates samples for exploitation by force-directed layout.
% Usage: xf = SAMPLING_FDL(xP,number,problem)
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

function xf = sampling_FDL(xP,number,problem)
    if (problem.control.verbose > 0)
        fprintf('Force-directed layout method...');
    end
    
    % Scale [0 1]
    nxP = size(xP,1);
    mxP = size(xP,2);
    xPlb = min(xP)-1e-10;
    xPub = max(xP)+1e-10;
    xPS = (xP - repmat(xPlb,nxP,1))./repmat((xPub - xPlb),nxP,1);
    
    xB = [];
    nbounding = ceil(0.2*number); % 20% for bounding points
    nothers = number - nbounding; % rest for interior points
    
    % Find bounding points
    xC = 0.5*ones(size(xPlb)); % Scaled center should be 0.5
    distC = sqrt(sum((xPS - repmat(xC,nxP,1)).^2,2));
    [~,idistC] = sort(distC,'descend'); % Get index in larger first
    nb = min(nbounding,length(distC));
    for i = 1:nb
        xB = [xB; xPS(idistC(i),:)];
    end
    
    % Find clusterred centers
    if (nothers >= size(xPS,1))
        nothers = size(xPS,1)-1;
    end
    if ((nothers > 0) && (size(xPS,1) > nothers))
        [~,xtmp,~] = kmeans(xPS,nothers);
    else
        xtmp = xPS;
    end
    xB = [xB; xtmp];
    number = size(xB,1);
    
    % Generate nearby points to the xB
    pm = rand(size(xB)); pm(pm<0.5) = -1; pm(pm>=0.5) = 1;
    xM = xB + 0.1*pm;
    vM = zeros(size(xM));
    
    if (problem.control.verbose == 2)
        fgFDL = figure();
        set(fgFDL,'Position',[150 100 540 840]);
    end
    for iloop = 1:1000

%         % For debug
%         if (problem.control.verbose == 2)
%             subplot(2,1,1); hold off; plot3(xPS(:,1),xPS(:,2),xPS(:,3),'rx'); hold on; plot3(xB(:,1),xB(:,2),xB(:,3),'bo'); plot3(xM(:,1),xM(:,2),xM(:,3),'k.'); axis([-0.2 1.2 -0.2 1.2 -0.2 1.2]); xlabel('x1'); ylabel('x2'); zlabel('x3');
%             for kk = 1:size(xB,1); plot3([xB(kk,1);xM(kk,1)],[xB(kk,2);xM(kk,2)],[xB(kk,3);xM(kk,3)],'k-'); end;
%             %for kk = 1:size(xM,1); for ll = 1:size(xM,1); if (kk ~= ll); plot3([xM(kk,1);xM(ll,1)],[xM(kk,2);xM(ll,2)],[xM(kk,3);xM(ll,3)],'c:'); end; end; end;
%             subplot(2,1,2); hold off; plot3(xPS(:,4),xPS(:,5),xPS(:,6),'rx'); hold on; plot3(xB(:,4),xB(:,5),xB(:,6),'bo'); plot3(xM(:,4),xM(:,5),xM(:,6),'k.'); axis([-0.2 1.2 -0.2 1.2 -0.2 1.2]); xlabel('x4'); ylabel('x5'); zlabel('x6');
%             for kk = 1:size(xB,1); plot3([xB(kk,4);xM(kk,4)],[xB(kk,5);xM(kk,5)],[xB(kk,6);xM(kk,6)],'k-'); end;
%         end
        
        % Run FDL
        % Peter Eades. A heuristic for graph drawing. Congressus Numerantium, 42:149–160, 1984.

        ForceSpring = zeros(number,mxP);

        % Attractive and repulsive forces between Base point (xB) and Moving point (xM)
        % Log spring force
        C1BM = 2; C2BM = 1.25/exp(1); % force balances if distance is 1.25
        distBM = sqrt(sum((xB - xM).^2,2));
        forceBM = C1BM * log(distBM./C2BM);
        vectorBM = (xB - xM)./repmat(distBM,1,mxP);
        ForceSpring = ForceSpring + vectorBM.*repmat(forceBM,1,mxP);
        
        % Repulsive forces between Pareto points (xPS) and Moving points (xM)
        % Inverse square law force
        C3PM = 1e-6;
        for i = 1:number
            for j = 1:nxP
                distPM = sqrt(sum((xM(i,:) - xPS(j,:)).^2));
                vectorPM = (xM(i,:) - xPS(j,:))./repmat(distPM,1,mxP);
                forcePM = C3PM/distPM^2;
                ForceSpring(i,:) = ForceSpring(i,:) + vectorPM * forcePM;
            end
        end
        
        % Dynamics
        mass = 1;
        timestep = 1/(1+log(iloop));
        vM = vM + ForceSpring/mass*timestep;
        vM(vM>0.1) = 0.1;
        vM(vM<-0.1) = -0.1;
        vM = vM / 1.1; % Velocity diminishing
        xM = xM + vM*timestep;
        
%         % For debug
%         if (problem.control.verbose == 2)
%             mindist = 1e3*ones(number,1);
%             minapart = 1e3*ones(number,1);
%             for i = 1:number
%                 for j = 1:nxP
%                     mindist(i) = min(mindist(i),sqrt(sum((xM(i,:) - xPS(j,:)).^2)));
%                 end
%             end
%             for i = 1:number
%                 for j = 1:number
%                     if ~(i == j)
%                         minapart(i) = min(minapart(i),sqrt(sum((xM(i,:)-xM(j,:)).^2)));
%                     end
%                 end
%             end
%             disp('minapart / mindist');
%             disp([minapart, mindist]);
%             if (iloop == 1)
%                 stdminapart = std(minapart);
%                 stdmindist = std(mindist);
%             end
%             drawnow;
%         end
        
        % Stopping
        if max(sqrt(sum(vM.^2,2))) < 1e-2
            break;
        end
    end
    
    % For debug
    if (problem.control.verbose == 2)
%         v = VideoWriter('FDL_Animation.mp4','MPEG-4');
%         v.FrameRate = 10;
%         open(v);
%         writeVideo(v,M);
%         close(v);
        close(fgFDL);
    end
    
    nxM = size(xM,1);
    xf = repmat(xPlb,nxM,1) + (repmat(xPub,nxM,1) - repmat(xPlb,nxM,1)) .* xM;
    
end