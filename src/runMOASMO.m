%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Main routine for running MO-ASMO algorithm.
%===============================================================================
function R = runMOASMO(settingsfun, objfun, nonlconfun, casefile)
    prob = createProblemStruct(settingsfun,objfun,nonlconfun,casefile);
    R = createResultStruct(prob); % create result struct
    rng(prob.random.seed,prob.random.generator); % set random seed
    %---------------------------------------------------------------------------
    % Plot preparation
    if (prob.control.plot)
        plotPreparation;
    end
    %---------------------------------------------------------------------------
    % sampling initial designs
    tic;
    k = 0;
    R.data(k+1,1).c01_smpX = samplingExploration(prob,k);
    R.data(k+1,1).time.c01_01 = toc;
    %---------------------------------------------------------------------------
    % main loop of adaptive refinement
    while (k<prob.control.maxiter)
    	k = k + 1;
        disp(strcat('[[Iteration:',num2str(k),']]'));

        % high fidelity function evaluation for initial/update training samples
        tic;
        R.data(k,1).c02_smpHffF = hffEvaluation(R.data(k,1).c01_smpX, prob);
        [R.data(k,1).c03_smpXFea, R.data(k,1).c05_smpXInfea, ...
         R.data(k,1).c04_smpHffFFea, R.data(k,1).c06_smpHffFInfea, ~, ~] ...
            = hffSeparateNaN( ...
            R.data(k,1).c01_smpX, R.data(k,1).c02_smpHffF, prob);
        R.data(k,1).time.c02_06 = toc;
        %-----------------------------------------------------------------------
        % Merging previous and current results to generate training point pool
        tic;
        if (k>1)
            R.data(k,1).c07_PoolXFea = [R.data(k-1,1).c07_PoolXFea;
                R.data(k-1,1).c27_valSurXFea; R.data(k,1).c03_smpXFea];
            R.data(k,1).c08_PoolHffFFea = [R.data(k-1,1).c08_PoolHffFFea;
                R.data(k-1,1).c29_valHffFFea; R.data(k,1).c04_smpHffFFea];
            R.data(k,1).c09_PoolXInfea = [R.data(k-1,1).c09_PoolXInfea;
                R.data(k-1,1).c30_valSurXInfea; R.data(k,1).c05_smpXInfea];
            R.data(k,1).c10_PoolHffFInfea = [R.data(k-1,1).c10_PoolHffFInfea;
                R.data(k-1,1).c32_valHffFInfea; R.data(k,1).c06_smpHffFInfea];
        else
            R.data(k,1).c07_PoolXFea = R.data(k,1).c03_smpXFea;
            R.data(k,1).c08_PoolHffFFea = R.data(k,1).c04_smpHffFFea;
            R.data(k,1).c09_PoolXInfea = R.data(k,1).c05_smpXInfea;
            R.data(k,1).c10_PoolHffFInfea = R.data(k,1).c06_smpHffFInfea;
        end
        R.data(k,1).time.c07_10 = toc;
        %-----------------------------------------------------------------------
        % Surrogate model construction
        tic;
        R.data(k,1).c11_surrogate = surrogateConstruct( ...
            R.data(k,1).c07_PoolXFea, R.data(k,1).c08_PoolHffFFea, prob);
        R.data(k,1).time.c11_11 = toc;
        %-----------------------------------------------------------------------
        % Infeasible region model construction
        tic;
        if (size(R.data(k,1).c09_PoolXInfea,1) ~= 0)
            R.data(k,1).c12_infeaRegion ...
                = infeaRegionConstruct(R.data(k,1).c09_PoolXInfea, prob);
        else
            R.data(k,1).c12_infeaRegion = [];
        end
        R.data(k,1).time.c12_12 = toc;
        %-----------------------------------------------------------------------
        % Generate initial population using previous results
        tic;
        if (k>1)
            % Previous approximated Pareto set
            tmp = R.data(k-1,1).c14_parSurX;
            % Sample from accumulated high fidelity results
            lx = size(R.data(k,1).c08_PoolHffFFea,1);
            for idx1 = 1:(prob.nfvar) % for each objective function variable
                [~,idx2] = sortrows(R.data(k,1).c08_PoolHffFFea,idx1);
                idx2 = idx2(1:(min(5,ceil(0.1*lx)))); % Try to get anchor points
                tmp = vertcat(tmp, R.data(k,1).c07_PoolXFea(idx2,:));
            end
            % Make each sample unique
            [R.data(k,1).c13_initpop,~,~] = unique(tmp,'rows');
            vars = {'tmp','lx','idx1','idx2'}; clear(vars{:}); clear vars;
        else
            R.data(k,1).c13_initpop = [];
        end
        R.data(k,1).time.c13_13 = toc;
        %-----------------------------------------------------------------------
        % Surrogate model-based optimization
        tic;
        [R.data(k,1).c14_parSurX, R.data(k,1).c15_parSurF, ...
         R.data(k,1).c16_parSurOut] = surrogateMOO( ...
            R.data(k,1).c11_surrogate, R.data(k,1).c12_infeaRegion, ...
            R.data(k,1).c13_initpop, prob);
        R.data(k,1).time.c14_16 = toc;
        %-----------------------------------------------------------------------
        % (Optional) High fidelity function evaluation
        % if the objective function evaluation is not expensive
        % We do not use these results for MO-ASMO computation,
        % but only for users to get information about exact errors
        tic;
        if (~prob.highfidelity.expensive)
            R.data(k,1).c17_parHffF = hffEvaluation( ...
                R.data(k,1).c14_parSurX, prob);
        else
            R.data(k,1).c17_parHffF = R.data(k,1).c15_parSurF;
        end
        [R.data(k,1).c18_parSurXFea, R.data(k,1).c21_parSurXInfea, ...
         R.data(k,1).c20_parHffFFea, R.data(k,1).c23_parHffFInfea, ...
         R.data(k,1).c19_parSurFFea, R.data(k,1).c22_parSurFInfea] ...
            = hffSeparateNaN(R.data(k,1).c14_parSurX, ...
            R.data(k,1).c17_parHffF, prob, R.data(k,1).c15_parSurF);
        R.data(k,1).time.c17_23 = toc;
        %-----------------------------------------------------------------------
        % Sample validation points, high fidelity function evaluation,
        % and compute error metrices.
        tic;
        [R.data(k,1).c24_valSurX, R.data(k,1).c25_valSurF] ...
            = samplingValidation( ...
            R.data(k,1).c14_parSurX, R.data(k,1).c15_parSurF, prob);
        R.data(k,1).c26_valHffF = hffEvaluation( ...
            R.data(k,1).c24_valSurX, prob);
        [R.data(k,1).c27_valSurXFea, R.data(k,1).c30_valSurXInfea, ...
         R.data(k,1).c29_valHffFFea, R.data(k,1).c32_valHffFInfea, ...
         R.data(k,1).c28_valSurFFea, R.data(k,1).c31_valSurFInfea] ...
            = hffSeparateNaN(R.data(k,1).c24_valSurX, ...
            R.data(k,1).c26_valHffF, prob, R.data(k,1).c25_valSurF);
        R.data(k,1).c33_valErrorVec ...
            = sqrt(sum((R.data(k,1).c29_valHffFFea ...
            - R.data(k,1).c28_valSurFFea).^2,2));
        R.data(k,1).c34_valErrorAvg ...
            = sum(R.data(k,1).c33_valErrorVec,1) ...
            / size(R.data(k,1).c33_valErrorVec,1);
        if (~prob.highfidelity.expensive)
            R.data(k,1).c35_parErrorVec ...
                = sqrt(sum((R.data(k,1).c20_parHffFFea ...
                - R.data(k,1).c19_parSurFFea).^2,2));
            R.data(k,1).c36_parErrorAvg ...
                = sum(R.data(k,1).c35_parErrorVec,1) ...
                / size(R.data(k,1).c35_parErrorVec,1);
        else
            R.data(k,1).c35_parErrorVec = [];
            R.data(k,1).c36_parErrorAvg = [];
        end
        R.data(k,1).time.c24_36 = toc;
        %-----------------------------------------------------------------------
        % Plot intermediate solution
        if (prob.control.plot)
            plotFig01Variables;
            plotFig02ObjFun;
            plotFig03ParetoEvol;
            plotFig04Convergence;
        end
        %-----------------------------------------------------------------------
        % Save intermediate result
        clear result; result = R.data(k,1);
        save(fullfile(prob.control.solpath, ...
            [prob.control.case, '_iter', num2str(k,'%04d'), '.mat']), 'result');
        %-----------------------------------------------------------------------
        % Termination condition check
        if ((R.data(k,1).c34_valErrorAvg < prob.control.maxerror) ...
                || (k == prob.control.maxiter))
            break;
        end
        %-----------------------------------------------------------------------
        % Sample update points
        tic;
        R.data(k+1,1).c01_smpX = [ ...
            samplingExploitation( ...
                R.data(k,1).c14_parSurX, prob);
            samplingExploration(prob, k, ...
                [R.data(k,1).c07_PoolXFea; R.data(k,1).c09_PoolXInfea]) ...
        ];
        R.data(k+1,1).time.c01_01 = toc;
    end
    %---------------------------------------------------------------------------
    % Save final result
    clear result; result = R;
    save(fullfile(prob.control.solpath, ...
        [prob.control.case, '_final.mat']), 'result');
end
%===============================================================================
