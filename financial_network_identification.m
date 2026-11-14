%% Identification of financial network LTV

clear all; close all; clc;
warning('off', 'all');
%% Data loading and preprocessing

dataset = readtable('dataset_MIMO_30.xlsx');
dates = datetime(dataset.Date,'InputFormat','yyyy-MM-dd HH:mm:ss');
u = dataset{:,2};       
y = dataset{:,3:end};   
Ts = 1; 
companies_names = {'AAL', 'ALGT', 'ALK', 'BKNG', 'CCL', 'COP', 'CPA', 'CSX', 'CVX', 'DAL', 'EOG', 'EXPE', 'FDX', 'JBLU', 'LUV', 'MAR', 'MPC', 'NCLH', 'NSC', 'OXY', 'RCL', 'SKYW', 'SLB', 'SNCY', 'UAL', 'UNP', 'UPS', 'VLO', 'XOM'};
N = size(y,1);
chunk_days= 120;
step = 1;
chunk_indices = 1:step:(N-chunk_days+1);
chunk_num = length(chunk_indices);
companies_num = 29; 
%% Global detrending

data = iddata(y, u, Ts);
mu = getTrend(data, 0);
data_d = detrend(data,mu);

%% ARX and analysis

na = 30; 
nb = 30; 
nk = 1;  
orders_siso = [na, nb, nk]; 
network_Y_Y = zeros(companies_num, companies_num, chunk_num);
network_U_Y = zeros(companies_num, 1, chunk_num);
opt_kernel = arxRegulOptions('RegulKernel', 'SS');


for k = 1:chunk_num
    
    y_chunk = y(k:(k+chunk_days-1), :);
    u_chunk = u(k:(k+chunk_days-1), :);
    
    %local scale normalization (variance balancing)
    y_chunk = y_chunk ./ std(y_chunk);
    u_chunk = u_chunk ./ std(u_chunk);
    
    temp_Y_Y = zeros(companies_num, companies_num);
    temp_U_Y = zeros(companies_num, 1);
    
    for i = 1:companies_num
        y_target = y_chunk(:, i);
        
        % exclude target density from exogenous inputs
        y_others = y_chunk;
        y_others(:, i) = []; 
        inputs_miso = [u_chunk, y_others];
        
        % detrending
        data_miso = iddata(y_target, inputs_miso, Ts);
        mu_miso = getTrend(data_miso, 0);
        det_data_miso = detrend(data_miso, mu_miso);
        Y_det = det_data_miso.y;
        U_det = det_data_miso.u;
        
        %kernel SS
        siso_init = iddata(Y_det, U_det(:, 1), Ts);
        [Lambda_opt, R_siso] = arxRegul(siso_init, orders_siso, opt_kernel);
        
        R_A = R_siso(1:na, 1:na);
        R_B = R_siso(na+1:end, na+1:end);
        
        R_base_blocks = cell(1, 1 + companies_num);
        R_base_blocks{1} = R_A;
        for j = 1:companies_num
            R_base_blocks{1+j} = R_B;
        end
        R_full = blkdiag(R_base_blocks{:});
        
        Lambda_miso = Lambda_opt * companies_num * 2;
        
        % construction of Phi
        max_delay = max(na, nb + nk - 1); 
        N_eff = chunk_days - max_delay;
        Phi = zeros(N_eff, na + companies_num * nb);
        Y_est = Y_det(max_delay + 1 : end); 
        
        for t = max_delay + 1 : chunk_days
            phi_A = -Y_det(t-1 : -1 : t-na)';
            phi_B = zeros(1, companies_num * nb);
            for j = 1:companies_num
                idx_start_B = (j-1)*nb + 1;
                idx_end_B = j*nb;
                phi_B(idx_start_B : idx_end_B) = U_det(t-nk : -1 : t-nk-nb+1, j)';
            end
            Phi(t - max_delay, :) = [phi_A, phi_B];
        end
        
        % RELS direct solution
        H = (Phi' * Phi) + Lambda_miso * (R_full' * R_full) + 1e-6 * eye(size(Phi, 2));
        theta_hat = H \ (Phi' * Y_est);
        
        % iterative reweighting
        epsilon = 0.5;   
        R_new_blocks = cell(1, 1 + companies_num); 
        theta_A = theta_hat(1:na); 
        R_new_blocks{1} = R_A / sqrt(epsilon + norm(R_A * theta_A)^2); 
        
        for j = 1:companies_num
            idx_start_B = na + (j-1)*nb + 1;
            idx_end_B = na + j*nb;
            theta_j = theta_hat(idx_start_B : idx_end_B); 
            R_new_blocks{1+j} = R_B / sqrt(epsilon + norm(R_B * theta_j)^2); 
        end
        
        R_new_full = blkdiag(R_new_blocks{:});
        H_iter = (Phi' * Phi) + Lambda_miso * (R_new_full' * R_new_full) + 1e-6 * eye(size(Phi, 2));
        theta_hat = H_iter \ (Phi' * Y_est);
        
        % sparse network extraction
        % auto regression (diagonal)
        temp_Y_Y(i, i) = norm(theta_hat(1:na));
        
        % fuel impact
        temp_U_Y(i, 1) = norm(theta_hat(na+1 : na+nb));
        
        % interactions Y -> Y 
        input_col = 2; 
        for j = 1:companies_num
            if j ~= i
                idx_start_B = na + (input_col-1)*nb + 1;
                idx_end_B = na + input_col*nb;
                temp_Y_Y(i, j) = norm(theta_hat(idx_start_B : idx_end_B));
                input_col = input_col + 1;
            end
        end
    end % <-- end cycle for i = 1:num_y
    
    network_Y_Y(:,:,k) = temp_Y_Y;
    network_U_Y(:,1,k) = temp_U_Y;

end

%% Dynamic thresholding


target_density = 8.5; 
adj_binary = zeros(companies_num, companies_num, chunk_num); 
net_stability = zeros(chunk_num, 1);
eps_evolution = zeros(chunk_num, 1); 

for k = 1:chunk_num
    curr_network = network_Y_Y(:,:,k);
    
    % create a new version of the matrix with zero on the main diagonal
    cross_network = curr_network;
    cross_network(logical(eye(companies_num))) = 0;
    
    % vector with values outside the diagonal
    off_diag_vals = cross_network(~eye(companies_num));
    
    % compute threshold to get the target density
    epsilon_thresh = prctile(off_diag_vals, 100 - target_density);
    epsilon_thresh = max(epsilon_thresh, 1e-6); 
    
   % adjoint binary matrix  
    adj_binary(:,:,k) = (cross_network >= epsilon_thresh);
    eps_evolution(k) = epsilon_thresh;
    
    % compute network stability
    if k > 1    
        prev_cross = network_Y_Y(:,:,k-1);
        prev_cross(logical(eye(companies_num))) = 0; 
        
        diff_matrix = cross_network .* adj_binary(:,:,k) - prev_cross .* adj_binary(:,:,k-1);  
        sum_squares = sum(diff_matrix(:).^2); 
        net_stability(k) = sqrt(sum_squares); 
    end
end


%% Shock date


shock_date = datetime('2022-02-24');
window_dates = dates(chunk_indices + chunk_days - 1);

margin_search = 20; 
date_start_search = shock_date - days(margin_search);
date_end_search   = shock_date + days(margin_search);

valid_indices = find(window_dates >= date_start_search & window_dates <= date_end_search);

[~, idx_max_var] = max(net_stability(valid_indices));
index_shock_day_identified = valid_indices(idx_max_var);

identified_shock_date = window_dates(index_shock_day_identified);
fprintf('The real shock happened on: %s\n', datestr(shock_date));
fprintf('The identified shock happened on: %s\n', datestr(identified_shock_date));

%% colors

sectors = {
    {'AAL', 'DAL', 'UAL', 'LUV', 'ALK', 'ALGT', 'CPA', 'JBLU', 'SKYW', 'SNCY'},...
    {'CCL', 'RCL', 'NCLH', 'BKNG', 'EXPE', 'MAR'},...
    {'XOM', 'CVX', 'COP', 'OXY', 'EOG', 'MPC', 'SLB', 'VLO'},...
    {'FDX', 'UPS', 'UNP', 'CSX', 'NSC'}
};
node_sectors = 5 * ones(companies_num, 1);
for i = 1:companies_num
    for s = 1:length(sectors)
       if ismember(companies_names{i}, sectors{s})
            node_sectors(i) = s; 
            break;
        end
    end
end
color_map_local = [0.85 0.32 0.09; 0.46 0.67 0.18; 0.92 0.69 0.12; 0.49 0.18 0.55; 0.3010 0.7450 0.9330]; 
node_colors_local = color_map_local(node_sectors, :);
%% Networks

offsets_days = [-90, -60, -30, -15, 0, 15, 30, 60, 90]; 
titles_grid = {'90 days before','60 days before','30 days before','15 days before','day of the shock','15 days later','30 days later','60 days later','90 days later'};

% matrix Y -> Y (inter-companies, without diagonal)
figure('Position', [50 50 1400 900], 'Name', 'matrices - network Y_Y');
sgtitle('Dynamical analysis: adjacency matrices (cross-network Y -> Y)', 'FontSize', 14, 'FontWeight', 'bold');
for p = 1:9
    subplot(3, 3, p); 
    k_idx = index_shock_day_identified + round(offsets_days(p) / step);
    k_idx = max(1, min(k_idx, chunk_num)); 
    
     % use cross matrices binary
    cross_net = network_Y_Y(:,:,k_idx);
    cross_net(logical(eye(companies_num))) = 0;
    weighted_matrix = cross_net .* adj_binary(:,:,k_idx);
    
    imagesc(weighted_matrix); 
    colormap('jet'); colorbar;
    title(titles_grid{p});
    xlabel('To company'); ylabel('from company');
    
    % add names on axes
    set(gca, 'XTick', 1:companies_num, 'XTickLabel', companies_names, 'XTickLabelRotation', 90, ...
             'YTick', 1:companies_num, 'YTickLabel', companies_names, 'FontSize', 6);
end

%% intermediate levels to zoom in, for a better understanding
% option 1: from 15 days before to day 0 
offsets_zoom = [-15, -12, -9, -6, -3, 0]; 
title_fig = 'pre-shock analysis (-15 to day 0)';

% option 2: from 0 to +15 days 
%offsets_zoom = [0, 3, 6, 9, 12, 15]; 
%title_fig = 'post shock analysis (day 0 to day 15)';

% option 3: from +30 to +60 days 
%offsets_zoom = [30, 36, 42, 48, 54, 60]; 
%title_fig= 'shock dissipation analysis (30 to 60 days later)';

figure('Position', [100 100 1300 700], 'Name', title_fig);
sgtitle(title_fig, 'FontSize', 14, 'FontWeight', 'bold');

for p = 1:length(offsets_zoom)
    subplot(2, 3, p); 
    
    % compute corresponding chunk
    k_idx = index_shock_day_identified + round(offsets_zoom(p) / step);
    k_idx = max(1, min(k_idx, chunk_num)); 
    
    % matrix extraction
    cross_net = network_Y_Y(:, :, k_idx);
    cross_net(logical(eye(companies_num))) = 0;
    weighted_matrix = cross_net .* adj_binary(:, :, k_idx);
    
    % graph
    imagesc(weighted_matrix); 
    colormap('jet'); 
    colorbar;
    
    % subplot titles
    if offsets_zoom(p) < 0
        title_str = sprintf('%d days before', abs(offsets_zoom(p)));
    elseif offsets_zoom(p) == 0
        title_str = 'shock day';
    else
        title_str = sprintf('+%d days after', offsets_zoom(p));
    end
    
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('to company'); 
    ylabel('from company');
    
    set(gca, 'XTick', 1:companies_num, 'XTickLabel', companies_names, 'XTickLabelRotation', 90, ...
             'YTick', 1:companies_num, 'YTickLabel', companies_names, 'FontSize', 6);
end

% matrices U -> Y (oil impact)
figure('Position', [150 150 1400 900], 'Name', 'Matrices network U_Y');
sgtitle('dynamical Analysis: fuel impact (U -> Y)', 'FontSize', 14, 'FontWeight', 'bold');
for p = 1:9
    subplot(3, 3, p); 
    k_idx = index_shock_day_identified + round(offsets_days(p) / step);
    k_idx = max(1, min(k_idx, chunk_num));
    
    oil_vector = network_U_Y(:, 1, k_idx);
    imagesc(oil_vector); 
    colormap('jet'); colorbar;
    title(titles_grid{p});
    set(gca, 'XTick', 1, 'XTickLabel', {'Oil'}, 'YTick', 1:companies_num, 'YTickLabel', companies_names, 'FontSize', 8);
    pbaspect([0.3 1 1]); 
end

%% intermediate levels to zoom in, for a better understanding
% [-90 to +90 DAYS]

%option 1: from 90 to -60 days before (long time before shock)
% offsets_zoom_U = [-90, -84, -78, -72, -66, -60]; 
% title_fig_U = 'Oil impact: far from shock (-90 to -60 days)';

% option 2: from 60 to 30 days before 
% offsets_zoom_U = [-60, -54, -48, -42, -36, -30]; 
% title_fig_U = 'Oil impact: from 60 to 30 days before ';

%option 3:from 30 to 15 days before
% offsets_zoom_U = [-30, -27, -24, -21, -18, -15]; 
% title_fig_U = 'Oil impact: from one month to two weeks earlier';

%option 4: from 15 days earlier to shock day 
% offsets_zoom_U = [-15, -12, -9, -6, -3, 0]; 
% title_fig_U = 'oil impact: immediate pre-shock (-15 to shock day)';

%option 5: from shock day to 2 days after 
offsets_zoom_U = [0, 3, 6, 9, 12, 15]; 
title_fig_U = 'Oil impact: immediate post-shock reaction (day 0 to +15)';

% option 6: from 15 days later to 30 days later 
% offsets_zoom_U = [15, 18, 21, 24, 27, 30]; 
% title_fig_U = 'Oil impact: +15 to +30 days';

% option 7:  from 30  to 60 days later  
% offsets_zoom_U = [30, 36, 42, 48, 54, 60]; 
% title_fig_U = 'Oil impact: from +30 to +60 days)';

%option 8: from 60 to 90 days after  
% offsets_zoom_U = [60, 66, 72, 78, 84, 90]; 
% title_fig_U = 'Oil impact: long-term stabilization (+60 to +90 Days)';


figure('Position', [150 150 1300 700], 'Name', title_fig_U);
sgtitle(title_fig_U, 'FontSize', 14, 'FontWeight', 'bold');

for p = 1:length(offsets_zoom_U)
    subplot(2, 3, p); 
    
    % chunk index
    k_idx = index_shock_day_identified + round(offsets_zoom_U(p) / step);
    k_idx = max(1, min(k_idx, chunk_num)); 
    
    % vector extraction (U -> Y)
    oil_vector = network_U_Y(:, 1, k_idx);
    
    % colormap
    imagesc(oil_vector); 
    colormap('jet'); 
    colorbar;
    
    % title subplot
    if offsets_zoom_U(p) < 0
        title_str = sprintf('%d Days Before', abs(offsets_zoom_U(p)));
    elseif offsets_zoom_U(p) == 0
        title_str = 'Day 0 (Shock)';
    else
        title_str = sprintf('+%d Days Later', offsets_zoom_U(p));
    end
    
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    % labels and axis
    set(gca, 'XTick', 1, 'XTickLabel', {'Oil'}, ...
             'YTick', 1:companies_num, 'YTickLabel', companies_names, ...
             'FontSize', 7);
         
    % change ratio for better visualisation
    pbaspect([0.3 1 1]); 
end

% autoregression matrices (memory Y_i -> Y_i)
figure('Position', [200 200 1400 900], 'Name', 'matrices - autoregression (Y_i -> Y_i)');
sgtitle('Dynamical analysis: autoregression (self-memory)', 'FontSize', 14, 'FontWeight', 'bold');

for p = 1:9
    subplot(3, 3, p); 
    k_idx = index_shock_day_identified + round(offsets_days(p) / step);
    k_idx = max(1, min(k_idx, chunk_num));
    
    ar_vector = diag(network_Y_Y(:,:,k_idx));
    imagesc(ar_vector); 
    colormap('jet'); colorbar;
    title(titles_grid{p});
    set(gca, 'XTick', 1, 'XTickLabel', {'AR (Memory)'}, 'YTick', 1:companies_num, 'YTickLabel', companies_names, 'FontSize', 8);
    pbaspect([0.3 1 1]); 
end
%% intermediate graphs for a better understanding

%option 1: from -30 to -15 days 
% offsets_zoom_AR = [-30, -27, -24, -21, -18, -15]; 
% title_fig_AR = 'AR self-memory: initial phase (-30 to -15 Days)';

% option 2: from 15 days earlier to shock day 
offsets_zoom_AR = [-15, -12, -9, -6, -3, 0]; 
title_fig_AR = 'AR self-memory: immediately before shock (-15 to shock day)';

%option 3: from shock day to 15 days after
% offsets_zoom_AR = [0, 3, 6, 9, 12, 15]; 
% title_fig_AR = 'AR self-memory: immediately after shock (day 0 to +15)';

% option4: from 15 to 30 days after
% offsets_zoom_AR = [15, 18, 21, 24, 27, 30]; 
% title_fig_AR = 'AR self-memory:  adjustment after shock (+15 to +30 Days)';


figure('Position', [200 200 1300 700], 'Name', title_fig_AR);
sgtitle(title_fig_AR, 'FontSize', 14, 'FontWeight', 'bold');

for p = 1:length(offsets_zoom_AR)
    subplot(2, 3, p); 
    
    %  chunk index
    k_idx = index_shock_day_identified + round(offsets_zoom_AR(p) / step);
    k_idx = max(1, min(k_idx, chunk_num)); 
    
    % diagonal extraction (autoregressive past)
   ar_vector = diag(network_Y_Y(:,:,k_idx)); 
    
    % single column color map
    imagesc(ar_vector); 
    colormap('jet'); 
    colorbar;
    
    % subplot title
    if offsets_zoom_AR(p) < 0
        title_str = sprintf('%d Days before', abs(offsets_zoom_AR(p)));
    elseif offsets_zoom_AR(p) == 0
        title_str = 'day 0 (Shock)';
    else
        title_str = sprintf('+%d Days later', offsets_zoom_AR(p));
    end
    
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    % axes and labels
    set(gca, 'XTick', 1, 'XTickLabel', {'AR (Memory)'}, ...
             'YTick', 1:companies_num, 'YTickLabel', companies_names, ...
             'FontSize', 7);
         
    % change ratio for better visualisation
    pbaspect([0.3 1 1]); 
end

% Bar plot U -> Y
figure('Position', [250 250 1600 900], 'Name', 'Bar plot - oil impact U_Y');
sgtitle('Ranking of oil impact on companies over time', 'FontSize', 14, 'FontWeight', 'bold');
for p = 1:9
    subplot(3, 3, p); 
    k_idx = index_shock_day_identified + round(offsets_days(p) / step);
    k_idx = max(1, min(k_idx, chunk_num));
    
    oil_vector = network_U_Y(:, 1, k_idx);
    b_oil = barh(oil_vector, 'FaceColor', 'flat');
    for i = 1:companies_num
        b_oil.CData(i, :) = node_colors_local(i, :);
    end
    set(gca, 'YTick', 1:companies_num, 'YTickLabel', companies_names, 'FontSize', 7);
    title(titles_grid{p});
    grid on; xlabel('Impact Weight');
end

%% Stability, epsilon and oil impact

figure('Name', 'Market Variation over time', 'Position', [150 150 900 400]); 
plot(window_dates, net_stability, 'b', 'LineWidth', 1.5); 
hold on; 
xline(shock_date, '--r', 'real shock', 'LineWidth', 1.5);
xline(identified_shock_date, '--g', 'Identified shock', 'LineWidth', 1.5);
grid on; title('Market variation over time (network stability)');
ylabel('Norm'); xlabel('Data'); legend('Location', 'best');

figure('Name', 'evolution of epsilon', 'Position', [200 200 900 400]);
plot(window_dates, eps_evolution, 'm', 'LineWidth', 2); 
hold on;
xline(shock_date, '--r', 'Real shock', 'LineWidth', 1.5);
xline(identified_shock_date, '--g', 'Identified shock', 'LineWidth', 1.5);
grid on; title('Evolution of dynamic threshold (Target Density)');
ylabel('Epsilon'); xlabel('Data'); legend('Location', 'best');

fuel_impact = zeros(chunk_num, 1);
for k = 1:chunk_num
    fuel_impact(k) = mean(network_U_Y(:, 1, k)); 
end

window_reaction_days = 60; 
search_limit = min(chunk_num, index_shock_day_identified + round(window_reaction_days / step));
[~, local_peak] = max(fuel_impact(index_shock_day_identified : search_limit));
chunk_revelation = index_shock_day_identified + local_peak - 1;

figure('Name', 'Fuel impact on the market', 'Position', [250 250 900 400]);
plot(window_dates, fuel_impact, 'b', 'LineWidth', 2);
hold on;
xline(shock_date, '--r', 'War shock', 'LineWidth', 1.5);
xline(window_dates(chunk_revelation), '--m', 'Fuel Peak', 'LineWidth', 1.5);
grid on; title('Fuel impact on the market');
ylabel('Mean impact'); xlabel('Data'); legend('Location', 'best');

if chunk_revelation > index_shock_day_identified
    days_for_reaction = (chunk_revelation - index_shock_day_identified) * step;
    fprintf('\nThe peak of oil impact happened in: %s\n', datestr(window_dates(chunk_revelation)));
    fprintf('Reaction time: %d working days.\n', days_for_reaction);
end

%% Hold-out cross validation

N_training = round(N * 0.75);
data_T = data_d(1 : N_training);    
data_V = data_d(N_training + 1 : end); 

na_miso = 30;  
nb_miso = 30; 
nk_miso = 1;
orders_miso_model = [na_miso, repmat(nb_miso, 1, companies_num), repmat(nk_miso, 1, companies_num)];

all_fits_z = zeros(companies_num, 1);

opt_comp = compareOptions('InitialCondition', 'z');

for i = 1:companies_num
    %Training data
    y_tr = data_T.y(:, i);
    y_others_tr = data_T.y;
    y_others_tr(:, i) = []; 
    inputs_tr = [data_T.u, y_others_tr];
    
    miso_train_obj = iddata(y_tr, inputs_tr, Ts);
    
    %kernel regulation
    [Lambda_opt_miso, R_miso_model] = arxRegul(miso_train_obj, orders_miso_model, opt_kernel);
    
    opt_arx_miso = arxOptions;
    opt_arx_miso.Regularization.Lambda = Lambda_opt_miso;
    opt_arx_miso.Regularization.R = R_miso_model;
    
    % model MISO
    m_miso_train = arx(miso_train_obj, orders_miso_model, opt_arx_miso);
    
    % validation data
    y_val = data_V.y(:, i);
    y_others_val = data_V.y;
    y_others_val(:, i) = []; 
    inputs_val = [data_V.u, y_others_val];
    
    miso_val_obj = iddata(y_val, inputs_val, Ts);
    
    % compute Fit  1-step ahead out-of-sample
    [~, fit_val_z, ~] = compare(miso_val_obj, m_miso_train, 1, opt_comp); 
    all_fits_z(i) = fit_val_z;
end

mean_fit_z = mean(all_fits_z);
fprintf('average global miso fit: %.2f%%\n', mean_fit_z);

% plot
figure('Name', 'Hold-out miso fit (1 step ahead)', 'Position', [100 300 1100 500]);
b_fit_z = bar(all_fits_z, 'FaceColor', 'flat');
for i = 1:companies_num
    b_fit_z.CData(i, :) = node_colors_local(i, :);
end
set(gca, 'XTick', 1:companies_num, 'XTickLabel', companies_names, 'XTickLabelRotation', 45, 'FontSize', 8);
grid on; 
title(sprintf('MISO hold out cross validation fit %% (1-step)- mean: %.2f%%', mean_fit_z));
ylabel('fit %'); 
xlabel('Company');
yline(mean_fit_z, '--r', ['mean (', num2str(mean_fit_z, '%.2f'), '%)'], 'LineWidth', 2);