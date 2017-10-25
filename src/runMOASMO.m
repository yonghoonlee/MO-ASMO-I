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
    c01_smpX = samplingExploration(prob,k);
    R.data.c01_smpX(k+1) = {c01_smpX};
    R.time.c01_01(k+1) = toc;
    %---------------------------------------------------------------------------
    % main loop of adaptive refinement
    while (k<prob.control.maxiter)
    	k = k + 1;
        disp(strcat('[[Iteration:',num2str(k),']]'));

        % high fidelity function evaluation for initial/update training samples
        tic;
        c02_smpHffF = hffEvaluation(c01_smpX, prob);
        [c03_smpXFea, c05_smpXInfea, c04_smpHffFFea, c06_smpHffFInfea, ~, ~] ...
            = hffSeparateNaN(c01_smpX, c02_smpHffF, prob);
        R.data.c02_smpHffF(k) = {c02_smpHffF};
        R.data.c03_smpXFea(k) = {c03_smpXFea};
        R.data.c04_smpHffFFea(k) = {c04_smpHffFFea};
        R.data.c05_smpXInfea(k) = {c05_smpXInfea};
        R.data.c06_smpHffFInfea(k) = {c06_smpHffFInfea};
        R.time.c02_06(k) = toc;
        %-----------------------------------------------------------------------
        % Merging previous and current results to generate training point pool
        tic;
        if (k>1)
            c07_old = cell2mat(R.data.c07_PoolXFea(k-1));
            c08_old = cell2mat(R.data.c08_PoolHffFFea(k-1));
            c09_old = cell2mat(R.data.c09_PoolXInfea(k-1));
            c10_old = cell2mat(R.data.c10_PoolHffFInfea(k-1));
            c27_old = cell2mat(R.data.c27_valSurXFea(k-1));
            c29_old = cell2mat(R.data.c29_valHffFFea(k-1));
            c30_old = cell2mat(R.data.c30_valSurXInfea(k-1));
            c32_old = cell2mat(R.data.c32_valHffFInfea(k-1));
            c07_PoolXFea = [c07_old; c27_old; c03_smpXFea];
            c08_PoolHffFFea = [c08_old; c29_old; c04_smpHffFFea];
            c09_PoolXInfea = [c09_old; c30_old; c05_smpXInfea];
            c10_PoolHffFInfea = [c10_old; c32_old; c06_smpHffFInfea];
        else
            c07_PoolXFea = c03_smpXFea;
            c08_PoolHffFFea = c04_smpHffFFea;
            c09_PoolXInfea = c05_smpXInfea;
            c10_PoolHffFInfea = c06_smpHffFInfea;
        end
        R.data.c07_PoolXFea(k) = {c07_PoolXFea};
        R.data.c08_PoolHffFFea(k) = {c08_PoolHffFFea};
        R.data.c09_PoolXInfea(k) = {c09_PoolXInfea};
        R.data.c10_PoolHffFInfea(k) = {c10_PoolHffFInfea};
        R.time.c07_10(k) = toc;
        %-----------------------------------------------------------------------
        % Surrogate model construction
        tic;
        c11_surrogate = surrogateConstruct( ...
            c07_PoolXFea, c08_PoolHffFFea, prob, k);
        R.data.c11_surrogate(k) = {c11_surrogate};
        R.time.c11_11(k) = toc;
        %-----------------------------------------------------------------------
        % Infeasible region model construction
        tic;
        if (size(c09_PoolXInfea,1) ~= 0)
            c12_infeaRegion = infeaRegionConstruct(c09_PoolXInfea, prob);
        else
            c12_infeaRegion = [];
        end
        R.data.c12_infeaRegion(k) = {c12_infeaRegion};
        R.time.c12_12(k) = toc;
        %-----------------------------------------------------------------------
        % Generate initial population using previous results
        tic;
        if (k>1)
            % Previous approximated Pareto set
            tmp = cell2mat(R.data.c14_parSurX(k-1));
            % Sample from accumulated high fidelity results
            lx = size(c08_PoolHffFFea,1);
            for idx1 = 1:(prob.nfvar) % for each objective function variable
                [~,idx2] = sortrows(c08_PoolHffFFea,idx1);
                idx2 = idx2(1:(min(5,ceil(0.1*lx)))); % Try to get anchor points
                tmp = vertcat(tmp, c07_PoolXFea(idx2,:));
            end
            % Make each sample unique
            [c13_initpop,~,~] = unique(tmp,'rows');
            vars = {'tmp','lx','idx1','idx2'}; clear(vars{:}); clear vars;
        else
            c13_initpop = [];
        end
        R.data.c13_initpop(k) = {c13_initpop};
        R.time.c13_13(k) = toc;
        %-----------------------------------------------------------------------
        % Surrogate model-based optimization
        tic;
        [c14_parSurX, c15_parSurF, c16_parSurOut] ...
            = surrogateMOO(c11_surrogate, c12_infeaRegion, c13_initpop, prob);
        R.data.c14_parSurX(k) = {c14_parSurX};
        R.data.c15_parSurF(k) = {c15_parSurF};
        R.data.c16_parSurOut(k) = {c16_parSurOut};
        R.time.c14_16(k) = toc;
        %-----------------------------------------------------------------------
        % (Optional) High fidelity function evaluation
        % if the objective function evaluation is not expensive
        % We do not use these results for MO-ASMO computation,
        % but only for users to get information about exact errors
        tic;
        if (~prob.highfidelity.expensive)
            c17_parHffF = hffEvaluation(c14_parSurX, prob);
        else
            c17_parHffF = c15_parSurF;
        end
        [c18_parSurXFea, c21_parSurXInfea, ...
         c20_parHffFFea, c23_parHffFInfea, ...
         c19_parSurFFea, c22_parSurFInfea] ...
            = hffSeparateNaN(c14_parSurX, c17_parHffF, prob, c15_parSurF);
        R.data.c17_parHffF(k) = {c17_parHffF};
        R.data.c18_parSurXFea(k) = {c18_parSurXFea};
        R.data.c19_parSurFFea(k) = {c19_parSurFFea};
        R.data.c20_parHffFFea(k) = {c20_parHffFFea};
        R.data.c21_parSurXInfea(k) = {c21_parSurXInfea};
        R.data.c22_parSurFInfea(k) = {c22_parSurFInfea};
        R.data.c23_parHffFInfea(k) = {c23_parHffFInfea};
        R.time.c17_23(k) = toc;
        %-----------------------------------------------------------------------
        % Sample validation points, high fidelity function evaluation
        tic;
        [c24_valSurX, c25_valSurF] ...
            = samplingValidation(c14_parSurX, c15_parSurF, prob);
        c26_valHffF = hffEvaluation(c24_valSurX, prob);
        [c27_valSurXFea, c30_valSurXInfea, ...
         c29_valHffFFea, c32_valHffFInfea, ...
         c28_valSurFFea, c31_valSurFInfea] ...
            = hffSeparateNaN(c24_valSurX, c26_valHffF, prob, c25_valSurF);
        R.data.c24_valSurX(k) = {c24_valSurX};
        R.data.c25_valSurF(k) = {c25_valSurF};
        R.data.c26_valHffF(k) = {c26_valHffF};
        R.data.c27_valSurXFea(k) = {c27_valSurXFea};
        R.data.c28_valSurFFea(k) = {c28_valSurFFea};
        R.data.c29_valHffFFea(k) = {c29_valHffFFea};
        R.data.c30_valSurXInfea(k) = {c30_valSurXInfea};
        R.data.c31_valSurFInfea(k) = {c31_valSurFInfea};
        R.data.c32_valHffFInfea(k) = {c32_valHffFInfea};
        R.time.c24_32(k) = toc;
        %-----------------------------------------------------------------------
        % Compute error metrices
        tic;
        if (k>1); eH = R.data.c33_valErrorVec; else; eH = {1}; end
        [c33_valErrorVec, c34_valErrorAvg] ...
            = errorEvaluation(c29_valHffFFea, c28_valSurFFea, ...
            prob, k, eH, c11_surrogate);
        % If not expensive
        if (~prob.highfidelity.expensive)
            if (k>1); eH = R.data.c35_parErrorVec; else; eH = {1}; end
            [c35_parErrorVec, c36_parErrorAvg] ...
                = errorEvaluation(c20_parHffFFea, c19_parSurFFea, ...
                prob, k, eH, c11_surrogate);
        else
            c35_parErrorVec = 0;
            c36_parErrorAvg = 0;
        end
        R.data.c33_valErrorVec(k) = {c33_valErrorVec};
        R.data.c34_valErrorAvg(k) = c34_valErrorAvg;
        R.data.c35_parErrorVec(k) = {c35_parErrorVec};
        R.data.c36_parErrorAvg(k) = c36_parErrorAvg;
        R.time.c33_36(k) = toc;
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
        clear result; result = R.data(k,:);
        save(fullfile(prob.control.solpath, ...
            [prob.control.case, '_iter', num2str(k,'%04d'), '.mat']), 'result');
        %-----------------------------------------------------------------------
        % Termination condition check
        if ((c34_valErrorAvg < prob.control.maxerror) ...
                || (k == prob.control.maxiter))
            break;
        end
        %-----------------------------------------------------------------------
        % Sample update points
        tic;
        c01_smpX = [ ...
            samplingExploitation(c14_parSurX, prob);
            samplingExploration(prob, k, [c07_PoolXFea; c09_PoolXInfea]) ...
        ];
        R.data.c01_smpX(k+1) = {c01_smpX};
        R.time.c01_01(k+1) = toc;
    end
    %---------------------------------------------------------------------------
    % Save final result
    clear result; result = R;
    save(fullfile(prob.control.solpath, ...
        [prob.control.case, '_final.mat']), 'result');
end
%===============================================================================
