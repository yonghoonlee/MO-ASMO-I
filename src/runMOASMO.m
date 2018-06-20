%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Main routine for running MO-ASMO algorithm.
%===============================================================================
function R = runMOASMO(settingsfun, objfun, nonlconfun, casefile, varargin)
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
    R.exitflag = 1;
    R.output.termination = 'Stopping criteria not met.';
    while (k<prob.control.maxiter)
    	k = k + 1;
        disp(strcat('[[Iteration:',num2str(k),']]'));

        % high fidelity function evaluation for initial/update training samples
        tic;
        c02_smpHffF = hffEvaluation(c01_smpX, prob);
        [c03_smpXFea, c05_smpXInfea, c04_smpHffFFea, c06_smpHffFInfea, ~, ~] ...
            = hffSeparateNaN(c01_smpX, c02_smpHffF, prob);
        R.data.c02_smpHffF(k,1) = {c02_smpHffF};
        R.data.c03_smpXFea(k,1) = {c03_smpXFea};
        R.data.c04_smpHffFFea(k,1) = {c04_smpHffFFea};
        R.data.c05_smpXInfea(k,1) = {c05_smpXInfea};
        R.data.c06_smpHffFInfea(k,1) = {c06_smpHffFInfea};
        R.time.c02_06(k,1) = toc;
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
        R.data.c07_PoolXFea(k,1) = {c07_PoolXFea};
        R.data.c08_PoolHffFFea(k,1) = {c08_PoolHffFFea};
        R.data.c09_PoolXInfea(k,1) = {c09_PoolXInfea};
        R.data.c10_PoolHffFInfea(k,1) = {c10_PoolHffFInfea};
        R.time.c07_10(k,1) = toc;
        %-----------------------------------------------------------------------
        % Surrogate model construction
        tic;
        c11_surrogate = surrogateConstruct( ...
            c07_PoolXFea, c08_PoolHffFFea, prob, k);
        R.data.c11_surrogate(k,1) = {c11_surrogate};
        R.time.c11_11(k,1) = toc;
        %-----------------------------------------------------------------------
        % Infeasible region model construction
        tic;
        if (size(c09_PoolXInfea,1) ~= 0)
            c12_infeaRegion = infeaRegionConstruct(c09_PoolXInfea, prob);
        else
            c12_infeaRegion = [];
        end
        R.data.c12_infeaRegion(k,1) = {c12_infeaRegion};
        R.time.c12_12(k,1) = toc;
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
        R.data.c13_initpop(k,1) = {c13_initpop};
        R.time.c13_13(k,1) = toc;
        %-----------------------------------------------------------------------
        % Surrogate model-based optimization
        tic;
        [c14_parSurX, c15_parSurF, c16_parSurOut] ...
            = surrogateMOO(c11_surrogate, c12_infeaRegion, c13_initpop, prob);
        R.data.c14_parSurX(k,1) = {c14_parSurX};
        R.data.c15_parSurF(k,1) = {c15_parSurF};
        R.data.c16_parSurOut(k,1) = {c16_parSurOut};
        R.time.c14_16(k,1) = toc;
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
        R.data.c17_parHffF(k,1) = {c17_parHffF};
        R.data.c18_parSurXFea(k,1) = {c18_parSurXFea};
        R.data.c19_parSurFFea(k,1) = {c19_parSurFFea};
        R.data.c20_parHffFFea(k,1) = {c20_parHffFFea};
        R.data.c21_parSurXInfea(k,1) = {c21_parSurXInfea};
        R.data.c22_parSurFInfea(k,1) = {c22_parSurFInfea};
        R.data.c23_parHffFInfea(k,1) = {c23_parHffFInfea};
        R.time.c17_23(k,1) = toc;
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
        R.data.c24_valSurX(k,1) = {c24_valSurX};
        R.data.c25_valSurF(k,1) = {c25_valSurF};
        R.data.c26_valHffF(k,1) = {c26_valHffF};
        R.data.c27_valSurXFea(k,1) = {c27_valSurXFea};
        R.data.c28_valSurFFea(k,1) = {c28_valSurFFea};
        R.data.c29_valHffFFea(k,1) = {c29_valHffFFea};
        R.data.c30_valSurXInfea(k,1) = {c30_valSurXInfea};
        R.data.c31_valSurFInfea(k,1) = {c31_valSurFInfea};
        R.data.c32_valHffFInfea(k,1) = {c32_valHffFInfea};
        R.time.c24_32(k,1) = toc;
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
        R.data.c33_valErrorVec(k,1) = {c33_valErrorVec};
        R.data.c34_valErrorAvg(k,1) = c34_valErrorAvg;
        R.data.c35_parErrorVec(k,1) = {c35_parErrorVec};
        R.data.c36_parErrorAvg(k,1) = c36_parErrorAvg;
        R.time.c33_36(k,1) = toc;
        %-----------------------------------------------------------------------
        % Compute convergence metrices
        tic;
        % 01-1: Hypervolume, Hypervolume Ratio in Predicted Pareto Set
        [R.data.c37_hypervolume(k,1), R.data.c38_hypervolumeRatio(k,1)] ...
            = approxNDHypervolume(R.data.c19_parSurFFea{k,1}, ...
                min(10*100^size(R.data.c19_parSurFFea{k,1},2),1000000));
        % 01-2: Normalized Hypervolume Change in Predicted Pareto Set
        if (k > 1)
            R.data.c39_changeHypervolume(k,1) = abs( ...
                R.data.c37_hypervolume(k-1,1) - R.data.c37_hypervolume(k,1)) ...
                / abs(max(R.data.c37_hypervolume(1:k)) ...
                - min(R.data.c37_hypervolume(1:k)));
            R.data.c40_changeHypervolumeRatio(k,1) = abs( ...
                R.data.c38_hypervolumeRatio(k-1,1) ...
                - R.data.c38_hypervolumeRatio(k,1)) ...
                / abs(max(R.data.c38_hypervolumeRatio(1:k)) ...
                - min(R.data.c38_hypervolumeRatio(1:k)));
        else
            R.data.c39_changeHypervolume(k,1) = 1;
            R.data.c40_changeHypervolumeRatio(k,1) = 1;
        end
        % 02-1: Hypervolume, Hypervolume Ratio in Pareto Set of HF Results
        [~,ndFHF,ndiHF] = ndSort(c07_PoolXFea,c08_PoolHffFFea);
        [R.data.c41_hypervolumeHF(k,1), R.data.c42_hypervolumeRatioHF(k,1)] ...
            = approxNDHypervolume(ndFHF(ndiHF==1,:), ...
            	min(10*100^size(R.data.c19_parSurFFea{k,1},2),1000000));
        % 02-2: Normalized Hypervolume Change in Pareto Set of HF Results
        if (k > 1)
            R.data.c43_changeHypervolumeHF(k,1) = abs( ...
                R.data.c41_hypervolumeHF(k-1,1) - R.data.c41_hypervolumeHF(k,1)) ...
                / abs(max(R.data.c41_hypervolumeHF(1:k)) ...
                - min(R.data.c41_hypervolumeHF(1:k)));
            R.data.c44_changeHypervolumeRatioHF(k,1) = abs( ...
                R.data.c42_hypervolumeRatioHF(k-1,1) ...
                - R.data.c42_hypervolumeRatioHF(k,1)) ...
                / abs(max(R.data.c42_hypervolumeRatioHF(1:k)) ...
                - min(R.data.c42_hypervolumeRatioHF(1:k)));
        else
        	R.data.c43_changeHypervolumeHF(k,1) = 1;
        	R.data.c44_changeHypervolumeRatioHF(k,1) = 1;
        end
        clear ndFHF ndiHF
        R.time.c37_44(k,1) = toc;
        %-----------------------------------------------------------------------
        % Plot intermediate solution
        if (prob.control.plot)
            plot_Fig01_02_Variables;
            plot_Fig03_03_ObjFun;
            plot_Fig04_04_ParetoEvol;
            plot_Fig05_05_EulerianErrorConvergence;
            plot_Fig06_07_Hypervolume;
            plot_Fig08_HFParetoSet;
            plot_Fig09_10_HypervolumeHF;
        end
        %-----------------------------------------------------------------------
        % Save intermediate result
        clear result; result = R.data(k,:);
        if (nargin == 6)
            outeriter = varargin{1};
            innercase = varargin{2};
            save(fullfile(prob.control.solpath, ...
                [prob.control.case, '_outeriter', num2str(outeriter,'%04d'), ...
                '_innercase', num2str(innercase,'%04d'), ...
                '_iter', num2str(k,'%04d'), '.mat']), ...
                'result','-v7.3');
        else
            save(fullfile(prob.control.solpath, ...
                [prob.control.case, '_iter', num2str(k,'%04d'), '.mat']), ...
                'result','-v7.3');
        end
        %-----------------------------------------------------------------------
        % Termination condition check
        if ((c34_valErrorAvg < prob.control.maxerror) ...
                || (k == prob.control.maxiter))
            R.exitflag = 0;
            R.output.termination = 'Stopping criteria met.';
            break;
        end
        %-----------------------------------------------------------------------
        % Sample update points
        tic;
        c01_smpX = [ ...
            samplingExploitation(c14_parSurX, prob);
            samplingExploration(prob, k, [c07_PoolXFea; c09_PoolXInfea]) ...
        ];
        R.data.c01_smpX(k+1,1) = {c01_smpX};
        R.time.c01_01(k+1,1) = toc;
    end
    %---------------------------------------------------------------------------
    % Plot final solution
    if (prob.control.plot)
        %
    end
    %---------------------------------------------------------------------------
    % Save final result
    clear result; result = R;    
    if (nargin == 6)
        outeriter = varargin{1};
        innercase = varargin{2};
        save(fullfile(prob.control.solpath, ...
            [prob.control.case, '_outeriter', num2str(outeriter,'%04d'), ...
            '_innercase', num2str(innercase,'%04d'), '_final.mat']), ...
            'result','-v7.3');
    else
        save(fullfile(prob.control.solpath, ...
            [prob.control.case, '_final.mat']), 'result','-v7.3');
    end
    [ndpop, ndscr, ndind] = ndSort(c07_PoolXFea, c08_PoolHffFFea);
    R.population = ndpop;
    R.scores = ndscr;
    R.output.iteration = k;
    R.x = ndpop(ndind==1,:);
    R.fval = ndscr(ndind==1,:);
end
%===============================================================================
