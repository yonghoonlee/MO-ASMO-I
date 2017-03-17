%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% MAIN SCRIPT
%==========================================================================
% Please run this script in the root path of MO-ASMO-I framework.
%==========================================================================

clear; clc; restoredefaultpath;

% Problem path is determined as "./Problems/{problem.name}"
problem.name = 'OsyczkaKundu';  % 'OsyczkaKundu', 'Lee2017a', 'Corman2016a'
                                % or user-defined problems can be solved.
problem = setup_problem(problem);               % Load problem definition
problem = setup_default_parameter(problem);     % Load default params
problem = feval(problem.settingsfun,problem);   % Problem-specific params

% Random seed control
randomseedfrontfactor = 0;                      % non-negative integers
rng(randomseedfrontfactor*100+0);               % Controlled random stream

% Description of DATA structure
% Column  1:  Samples generated at each iteration
% Column  2:  High fidelity function evaluation results of (Col. 1)
% Column  3:  Merged samples
% Column  4:  Merged high fidelity function evaluation results of (Col. 3)

% Column 13:  Time measurement data
%   sub-column  1:  High-fidelity function evaluation for generated samples
%   sub-column  2:  Surrogate model SM construction using data (Col. 3,4)
%   sub-column  3:  Initial population generation using data (Col. 3,4,5)
%   sub-column  4:  Multiobjective optimization on the surrogate model SM
% Column 14:  Initial population given to the genetic algorithm

% Initialize
k = 0;
disp(strcat('[[Iteration:',num2str(k),']]'));
DATA{k+1,1} = sampling_exploration(problem,k);  % Sampling initial designs
                                                % (constraints enforced)

while(k<problem.control.maxiter)
    % Iteration begins
    k = k + 1;
    
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
        DATA{k,14} = initpop;                   % Save initial populations
    else
        initpop = [];
    end
    DATA{k,13} = [DATA{k,13}, toc];
    
    % Multiobjective optimization running on the surrogate model
    tic;
    MO{k} = surrogate_gamultiobj(SM{k,1},problem,initpop);
    DATA{k,5} = MO{k}.x;
    DATA{k,6} = MO{k}.f;
    DATA{k,13} = [DATA{k,13}, toc];
    
    % (Optional) high-fidelity fn evaluation for the predicted Pareto set
    if (~problem.highfidelity.expensive)
    end
end