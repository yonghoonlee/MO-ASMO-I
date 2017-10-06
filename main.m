% The MAIN script contains the outer loop of MO-ASMO algorithm.
% Users can run this script after modifying the problem.name variable.
% Usage: MAIN

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

clear; clc; close all; restoredefaultpath;

% Problem path is determined as "./Problems/{problem.name}"
problem.name = 'VESuspensionContSpectra';
                        % Options:  problem.name = 'OsyczkaKundu';
                        %           problem.name = 'Lee2017a';
                        %           problem.name = 'VESuspension';
                        %           problem.name = 'VESuspensionVarFreq';
                        %           problem.name = 'GiesekusModel';
                        %           problem.name = 'Giesekus3Dv03';
                        %           problem.name = 'VESuspensionContSpectra';
                        %           or user-defined problems can be solved.
% Process problem structure
problem = setup_problem(problem);               % Load problem definition
problem = setup_default_parameter(problem);     % Load default params
problem = feval(problem.settingsfun,problem);   % Problem-specific params

% Random seed control
randomseedfrontfactor = 1;                      % non-negative integers
rng(randomseedfrontfactor*100+0);               % Controlled random stream

% Prepare plotting
plot_preparation;                               % Plot preparation script

% Description of DATA structure
% Column  1:  Samples generated at each iteration
% Column  2:  High fidelity function evaluation results of (Col. 1)
% Column  3:  Merged samples
% Column  4:  Merged high fidelity function evaluation results of (Col. 3)
% Column  5:  Predicted Pareto set in design space
% Column  6:  Predicted Pareto set in objective function space
% Column  7:  High fidelity function evaluation results of (Col. 5)
% Column  8:  Error between (Col. 6) and (Col. 7)
% Column  9:  Validation samples in design space
% Column 10:  Validation samples in objective function space
% Column 11:  High fidelity function evaluation results of (Col. 9)
% Column 12:  Error between (Col. 10) and (Col. 11)
% Column 13:  Time measurement data
%   sub-column  1:  High-fidelity function evaluation for generated samples
%   sub-column  2:  Surrogate model SM construction using data (Col. 3,4)
%   sub-column  3:  Initial population generation using data (Col. 3,4,5)
%   sub-column  4:  Multiobjective optimization on the surrogate model SM
%   sub-column  5:  (High-fidelity fn evaluation for predicted Pareto set)
%   sub-column  6:  Sampling for validation
%   sub-column  7:  Sampling for update
% Column 14:  Initial population given to the genetic algorithm

% Initialize
k = 0;
DATA{k+1,1} = sampling_exploration(problem,k);  % Sampling initial designs
                                                % (constraints enforced)

while(k<problem.control.maxiter)
    % Iteration begins
    k = k + 1;
    disp(strcat('[[Iteration:',num2str(k),']]'));
    
    % High fidelity simulation for the first objective computation
    tic;
    DATA{k,2} = hff_eval(DATA{k,1},problem);    % High-fidelity fn eval
    [DATA{k,1},DATA{k,2}] = util_removeNAN(DATA{k,1},DATA{k,2},problem);
    DATA{k,13} = toc;
    
    % Merging pervious simulation results and current simulation results
    % for generating training point pool
    if (k>1)
        DATA{k,3} = [DATA{k-1,3}; DATA{k-1,9}; DATA{k,1}];
        DATA{k,4} = [DATA{k-1,4}; DATA{k-1,11}; DATA{k,2}];
    else
        DATA{k,3} = DATA{k,1};
        DATA{k,4} = DATA{k,2};
    end
    
    % Surrogate model construction
    tic;
    SM{k,1} = surrogate_construct(DATA{k,3},DATA{k,4},problem);
    DATA{k,13} = [DATA{k,13}, toc];
    
    % Generate initial population using previous results
    tic;
    if (k>1)
        % Initial population sampling from previous approximated Pareto set
        initpop = DATA{k-1,5};
        % Initial population sampling from previous high-fidelity results
        ldx = length(DATA{k,4}(:,1));
        for i = 1:size(DATA{k,4},2)             % For each obj variable
            [~,idx] = sortrows(DATA{k,4},i);
            idx = idx(1:ceil(0.1*ldx));
            initpop = vertcat(initpop, DATA{k,3}(idx,:));
        end
        [initpop,~,~] = unique(initpop,'rows'); % Make each item unique
    else
        initpop = [];
    end
    DATA{k,14} = initpop;                       % Save initial populations
    DATA{k,13} = [DATA{k,13}, toc];
    
    % Multiobjective optimization running on the surrogate model
    tic;
    MO{k} = surrogate_gamultiobj(SM{k,1},problem,initpop);
    tmpdat = [MO{k}.f, MO{k}.x];
    tmpdat = sortrows(tmpdat,1:problem.nfvar);  % Sort data in obj fn value
    DATA{k,5} = tmpdat(:,(problem.nfvar + 1):end);
    DATA{k,6} = tmpdat(:,1:(problem.nfvar));
    clear tmpdat;
    DATA{k,13} = [DATA{k,13}, toc];
    
    % (Optional) high fidelity fn evaluation for the predicted Pareto set
    tic;
    if (~problem.highfidelity.expensive)
        DATA{k,7} = hff_eval(DATA{k,5},problem);% High fidelity fn eval
        [DATA{k,5},DATA{k,7},DATA{k,6}] ...
            = util_removeNAN(DATA{k,5},DATA{k,7},problem,DATA{k,6});
        DATA{k,8} = sum(sqrt(sum((DATA{k,7}-DATA{k,6}).^2,2)))...
            /size(DATA{k,7},1);                 % Error computation
    end
    DATA{k,13} = [DATA{k,13}, toc];
    
    % Generate figures of predicted Pareto set
    if (problem.control.plot ~= 0)
        plot_pareto;
    end
    
    % Validation sampling
    tic;
    [DATA{k,9},DATA{k,10}] = sampling_val(DATA{k,5},DATA{k,6},problem);
    DATA{k,11} = hff_eval(DATA{k,9},problem);   % High fidelity fn eval
    [DATA{k,9},DATA{k,11},DATA{k,10}] ...
        = util_removeNAN(DATA{k,9},DATA{k,11},problem,DATA{k,10});
    DATA{k,12} = sum(sqrt(sum((DATA{k,11}-DATA{k,10}).^2,2)))...
        /size(DATA{k,11},1);                    % Error computation
    DATA{k,13} = [DATA{k,13}, toc];
    
    % Generate figures of validation sample points
    if (problem.control.plot ~= 0)
        plot_validation;
        if isfield(problem,'plotcustom')
            for iplotcustom = 1:length(problem.plotcustom)
                figure(fgcustom{iplotcustom}); clf;
                run(fullfile(problem.probpath,problem.plotcustom{iplotcustom}));
                if (problem.control.plotexport ~= 0)
                    eval(['export_fig ',fullfile(problem.probpath,...
                        ['fig_custom',num2str(iplotcustom),'_',num2str(k),'.pdf']),...
                        ' -pdf']);
                end
            end
        end
    end
    
    % Update sampling
    tic;
    prevPoints = [DATA{k,3};DATA{k,9}];
    DATA{k+1,1} = [sampling_exploitation(DATA{k,5},problem);
                   sampling_exploration(problem,k,[DATA{k,3};DATA{k,9}])];
    DATA{k,13} = [DATA{k,13}, toc];
    
    if (problem.control.savedata ~= 0)
        for i = 1:14
            D{1,i} = DATA{k,i};
        end
        save([fullfile(problem.probpath,'DATA_'),num2str(k)],'D','k');
    end
    
end