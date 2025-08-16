% surfstat_two_step_negonly.m
% -------------------------------------------------------------------------
% Two-step SurfStat pipeline (neg-only):
%   Step 1: regress out BMI -> residuals
%   Step 2: residuals ~ adiposity
% Only negative association is tested and visualized.
% Cluster-forming threshold = 0.01
% Corrected P-masks: voxel-P < 0.001, cluster-P < 0.01
% Outputs: standardized beta, cluster CSV, cluster PNGs
% -------------------------------------------------------------------------

%% Configuration Section
CONFIG.DATA_DIR   = 'example_data';
CONFIG.MORPH_DIR  = fullfile('example_data','T1');
CONFIG.OUT_DIR    = fullfile('outputs','surfstat');
CONFIG.LH_PIAL    = 'lh.pial';
CONFIG.RH_PIAL    = 'rh.pial';  
CONFIG.PARTS      = {'arm','leg','trunk','vat'};
CONFIG.MEASUREMENTS = {'thickness','volume','area'};
CONFIG.AGG_FMT    = 'SUBJ%d.%s.00.mgh';  
CONFIG.VOXEL_P    = 1e-3;
CONFIG.CLUSTER_P  = 1e-2;
CONFIG.CLUSTER_THRESH = 0.01;
CONFIG.CLIM_STD_BETA = [-0.15, 0.15];
CONFIG.CSV_PATH = @(part) fullfile(CONFIG.DATA_DIR, sprintf('%s_with_covariates.csv', part));

%% Load surface and create output directories
avsurf = SurfStatReadSurf({CONFIG.LH_PIAL, CONFIG.RH_PIAL});
if ~exist(CONFIG.OUT_DIR, 'dir'); mkdir(CONFIG.OUT_DIR); end

%% Main analysis loop for each adiposity variable ("part")
for ip = 1:numel(CONFIG.PARTS)
    part = CONFIG.PARTS{ip};
    T = readtable(CONFIG.CSV_PATH(part));
    assert(ismember('BMI', T.Properties.VariableNames), 'BMI column missing.');
    assert(ismember(part, T.Properties.VariableNames), 'Adiposity variable missing.');

    X_adip = T.(part);
    BMI    = T.BMI;
    nSubj  = height(T);

    % Load cortical metrics
    Y = struct();
    for im = 1:numel(CONFIG.MEASUREMENTS)
        meas = CONFIG.MEASUREMENTS{im};
        lh_file = fullfile(CONFIG.MORPH_DIR, sprintf('lh.%s', sprintf(CONFIG.AGG_FMT, nSubj, meas)));
        rh_file = fullfile(CONFIG.MORPH_DIR, sprintf('rh.%s', sprintf(CONFIG.AGG_FMT, nSubj, meas)));
        Y.(meas) = SurfStatReadData({lh_file, rh_file});
    end

    % STEP 1: Regress out BMI -> residuals
    R = struct();
    for im = 1:numel(CONFIG.MEASUREMENTS)
        meas = CONFIG.MEASUREMENTS{im};
        model1 = 1 + term(BMI);
        slm1 = SurfStatLinMod(Y.(meas), model1, avsurf);
        R.(meas) = Y.(meas) - slm1.X * slm1.coef;
    end

    % STEP 2: residual ~ adiposity
    for im = 1:numel(CONFIG.MEASUREMENTS)
        meas = CONFIG.MEASUREMENTS{im};
        
        % Extract covariates (excluding BMI and adiposity variable)
        covariate_names = setdiff(T.Properties.VariableNames, {'eid', part, 'BMI'});
        X_covars = table2array(T(:, covariate_names));
        
        % Convert to SurfStat terms
        model2 = 1 + term(X_adip);
        for k = 1:size(X_covars, 2)
            model2 = model2 + term(X_covars(:, k));
        end
        slm2 = SurfStatLinMod(R.(meas), model2, avsurf);

        % Negative direction only  
        slm_neg = SurfStatT(slm2, -[0 1]);  
        [pval, peak, clus, clusid] = SurfStatP(slm_neg, [], CONFIG.CLUSTER_THRESH);  

        % Dual-threshold mask
        voxel_mask = pval.P < CONFIG.VOXEL_P;
        clust_mask = pval.C < CONFIG.CLUSTER_P;
        sig_mask = voxel_mask & clust_mask;

        % Standardized beta
        sX = std(X_adip);
        sY = std(R.(meas), 0, 1);
        beta_unstd = slm_neg.coef(2,:);
        beta_std = (beta_unstd .* sX) ./ sY;
        beta_std(~isfinite(beta_std)) = 0;
        beta_std_sig = beta_std .* sig_mask;  

        % Save results
        out_path = fullfile(CONFIG.OUT_DIR, part);
        if ~exist(out_path, 'dir'); mkdir(out_path); end
        save(fullfile(out_path, sprintf('%s_%s_results.mat', part, meas)), ...
            'slm_neg', 'pval', 'clus', 'clusid', 'sig_mask', 'beta_std', 'beta_std_sig');

        % Save cluster CSV
        if ~isempty(clus)
            writetable(struct2table(clus), fullfile(out_path, sprintf('%s_%s_clus.csv', part, meas)));
        end

        % Visualization
        fig = figure('Color','w','Position',[100 100 900 700]);
        SurfStatView(beta_std_sig, avsurf, sprintf('%s - %s (StdBeta)', part, meas));
        colormap(fig, make_diverging_cmap());
        caxis(CONFIG.CLIM_STD_BETA);
        exportgraphics(fig, fullfile(out_path, sprintf('%s_%s_stdBeta.png', part, meas)), 'Resolution', 300);
        close(fig);

        % Save each significant cluster as PNG
        sig_ids = clus(1).clusid(clus(1).P < CONFIG.CLUSTER_P);
        for cid = sig_ids'
            vm = (clusid == cid) & sig_mask;
            if ~any(vm(:)), continue; end
            f = figure('Color','w','Position',[80 60 900 700]);  
            SurfStatView(double(vm), avsurf, sprintf('%s - %s - clus_%d', part, meas, cid));
            colormap(f, gray);
            exportgraphics(f, fullfile(out_path, sprintf('clus_%s_%s_%d.png', part, meas, cid)), 'Resolution', 300);
            close(f);
        end
    end

    fprintf(' Finished: %s\n', part);  
end

fprintf('\n All analysis complete. Results saved in: %s\n', CONFIG.OUT_DIR);  

%% Custom colormap function
function cmap = make_diverging_cmap()
    base = [ ...
        0.4039 0.0000 0.1216
        0.6980 0.0941 0.1686
        0.8392 0.3765 0.3020
        0.9569 0.6471 0.5098
        0.9922 0.8588 0.7804
        1.0000 1.0000 1.0000
        0.8196 0.8980 0.9412
        0.5725 0.7725 0.8706
        0.2627 0.5765 0.7647
        0.1294 0.4000 0.6745
        0.0196 0.1882 0.3804 ];
    n = 256;
    cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    cmap = flipud(cmap);
    mid = round(size(cmap,1)/2);
    cmap(mid,:) = [1 1 1];
end
