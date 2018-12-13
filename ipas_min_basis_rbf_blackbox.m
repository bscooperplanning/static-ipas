% This blackbox can be used to test our IPAS planner for different number
% of sensors, number of Guassian basis (parameters), and different grid
% sizes (as squared function of grid row size). 
% See 'blackbox_sweep_example.m' for example using IPAS approach given a
% list of # of sensors used.
%
% INPUTS
%   - n_sensors: number of sensors to be used in each iteration
%   - n_params: number of Gaussian basis functions that compose the field,
%   also setup as perfect square to evenly distribute around field.
%   Therefore: Total basis functions = n_params^2
%   - n_grid_dim: Assuming a square grid workspace, this is the row/column
%   size of the grid. So, total grid size = n_grid_row ^ 2
%   - SEED: true or false. If you want repeatble fields SEED = true
%
% OUTPUTS
%   - path_cost_incurred: the resulting path cost from our iterative
%   approach. Should always be >= path_cost_true
%   - path_cost_true: the path cost from planning with full information about
%   the environment. Should be lower bound for path_cost_incurred. 

function [path_cost_incurred, path_cost_expected, path_cost_true, total_measurements, path_cost_var, conv_iters] = ...
    ipas_min_basis_rbf_blackbox(n_sensor, threat_basis_data, grid_world, sensor_noise)

%% TABULA RASA - the clean slate club
% clear variables; close all; clc
addpath('./routine/') % Add greedy sensor placement algorithms
% From J. Ranieri "Near Optimal Sensor Placement for Linear Inverse
% Problems", https://github.com/jranieri/OptimalSensorPlacement

threat_parameters_true = threat_basis_data.threat_parameters_true;
Np = threat_basis_data.n_threat_parameters;
Ng = grid_world.n_grid_points;
%% SENSORS

%===== SENSOR PLACEMENT
sensor_locations = zeros(1, n_sensor);

%----- Random selection of Ns basis from Np^2 basis
if (n_sensor > Np)
    rand_basis = [randperm(Np, Np) randi(Np, 1, n_sensor - Np)];
else
    rand_basis = randperm(Np, n_sensor);
end
for m1 = 1:length(rand_basis)
    % Convert basis center to nearest grid point
    x_loc = threat_basis_data.basis_parameters.mean(1,rand_basis(m1));
    y_loc = threat_basis_data.basis_parameters.mean(2,rand_basis(m1));
    grid_vec = kron([x_loc, y_loc]', ones(1, grid_world.n_grid_points)) - ...
            grid_world.coordinates;
    grid_xy_sep = sqrt(grid_vec(1,:).^2 + grid_vec(2,:).^2);
    [~, grid_id] = min(grid_xy_sep);
    sensor_locations(m1) = grid_id;
end

%----- Arbitary placement
% for m1 = 1:n_sensor
% 	this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
% 	while any(this_location == sensor_locations(1:m1))
% 		this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
% 	end
% 	sensor_locations(m1) = this_location;
% end
% % Or empty initial set to allow for initial optimistic path
% sensor_locations = [1];

%% PATH PLANNING WITH TRUE THREAT (FOR COMPARISON)
v_start	  = 1;
v_goal	  = grid_world.n_grid_points;

threat_value_true = calc_threat_rbf(threat_basis_data, ...
	threat_parameters_true, grid_world.coordinates);

transition_cost = grid_world.adjacency;
for m1 = 1:grid_world.n_grid_points
	nhbrs	= find(grid_world.adjacency(m1, :));
	transition_cost(m1, nhbrs) = max(1E-3, threat_value_true(:, nhbrs)); % If exposure less than 0, take 0
end
search_data.adjacency_matrix= transition_cost;
search_data.adjacency_struct= [];
search_data.fcn_find_nhbr	= [];

search_data.v_start			= v_start; % start/goal verts set above sensor initial config
search_data.v_goal			= v_goal;
search_data.all_goal		= false;
% search_data.heuristic		= zeros(grid_world.n_grid_points, 1);
search_data.heuristic       = @min_threat_manhat_heur;
search_data.grid_world      = grid_world;

% [vertex_data_true, ~]	= astar(search_data);
[vertex_data_true, ~]	= astar_hfunc(search_data); % Uses heurstic function
opt_path_true = greedy_trace(search_data.v_start, search_data.v_goal, vertex_data_true);

path_cost_true       		= sum(threat_value_true(opt_path_true(2:end)));

%% PATH PLANNING WITH ESTIMATED THREAT IN LOOP WITH SENSOR RECONFIGURATION

iters = 100; % max iterations before artificial termination
incurred_rec = zeros(1,1);
expected_rec = zeros(1,1);
path_cost_var_rec = zeros(1,1);

% Binary sensor vector [Ng^2 x iters]
bin_sensor_locations = zeros(grid_world.n_grid_points,iters);
bin_sensor_locations(sensor_locations,1) = 1;

%========================================================================
%============= Initial Batch Estimation of Threat Parameters ============
%========================================================================
% Find basis which can be best estimated 
identifiable_basis = identify_basis_close2sensor(sensor_locations, grid_world, threat_basis_data);
threat_basis_data_subset = threat_basis_data;
threat_basis_data_subset.basis_subset = identifiable_basis;
threat_basis_data_subset.basis_parameters.mean = ...
    threat_basis_data.basis_parameters.mean(:, identifiable_basis);
threat_basis_data_subset.n_threat_parameters = numel(identifiable_basis);

[tpe_mean, tpe_covar] = threat_estimate_rbf(sensor_locations, grid_world, ...
	threat_basis_data_subset, threat_parameters_true(identifiable_basis), sensor_noise);

% Reconstruct Field using only parameters with confidence
threat_value_estimate.mean	= calc_threat_rbf(threat_basis_data_subset, tpe_mean, grid_world.coordinates);

%===== Path planning
transition_cost_estimate = grid_world.adjacency;
for m1 = 1:grid_world.n_grid_points
    nhbrs	= find(grid_world.adjacency(m1, :));
    % intended to exclude negative vals? but erases connectivity for neg vals
    transition_cost_estimate(m1, nhbrs) = max( 1E-3, (threat_value_estimate.mean(:, nhbrs)) );
end
search_data.adjacency_matrix= transition_cost_estimate;

% [vertex_data, ~]= astar(search_data);
[vertex_data, ~]= astar_hfunc(search_data); % Uses heurstic function
opt_path		= greedy_trace(search_data.v_start, search_data.v_goal, vertex_data);

%===== Key numbers
expected_cost_estimate	= sum( max(0, threat_value_estimate.mean(opt_path(2:end))) );
incurred_cost			= sum(threat_value_true(opt_path(2:end)));

incurred_rec(1)= incurred_cost;
expected_rec(1)= expected_cost_estimate;
path_cost_var_rec(1) = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar);

% tpe_covar_full is the minimum covariance matrix possible if every grid
% point was measured. Serves as the lower bound on paremeter estimation and
% resulting path cost variance.
[~, tpe_covar_full] = threat_estimate_rbf(1:grid_world.n_grid_points,...
    grid_world, threat_basis_data, threat_parameters_true, sensor_noise);
min_discrete_path_variance = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar_full);
k = 5; % constant which limits the sub-optimality of IPAS path variance
m00 = 1;
while(abs(path_cost_var_rec(m00)) > k^2*min_discrete_path_variance && m00 < iters) 
    %===== Sensor reconfiguration
    % Place sensors along path, grid points closest to basis of interest
    % [new_locations, tp_oi] = actor_sensor_reconf_01(threat_basis_data, grid_world, n_sensor, opt_path);
    % Places sensors along path, prioritize the unmeasured basis
    % [new_locations, tp_intsct] = sensor_reconf_basis_unmeasured(threat_basis_data,...
    %     threat_basis_data_subset, grid_world, n_sensor, opt_path);
    % Place sensors along path, place minimum # of sensors to cover unmeasured basis
    [new_locations, tp_intsct] = sensor_reconf_min_basis(threat_basis_data,...
        threat_basis_data_subset, grid_world, n_sensor, opt_path);
    % Using Frame Potential and greedy selection for basis subdomain chosen by
    % cross-correlation
    %     [sensor_locations, tp_oi] = xcorr_FP_sensor_reconf(threat_basis_data, grid_world, n_sensor, opt_path);
    % Using proximity to basis subset, where basis subset chosen by cross-corr
    % 	[sensor_locations, tp_oi] = xcorr_sensor_reconf_v2_rbf(threat_basis_data, grid_world, n_sensor, opt_path);
    % Using fft based xcorr2
    %     [new_locations, tp_oi] = xcorr_fft_sensor_reconf_rbf(threat_basis_data, grid_world, n_sensor, opt_path);
    sensor_locations = [sensor_locations new_locations];
    bin_sensor_locations(sensor_locations,m00 + 1) = 1;
    norm_sensor_diff(m00) = norm(bin_sensor_locations(:,m00+1) - bin_sensor_locations(:,m00),1);
	%===== Threat estimation using Batch weighted LS
    % Find basis which can be best estimated 
    identifiable_basis = identify_basis_close2sensor(sensor_locations, grid_world, threat_basis_data);
    threat_basis_data_subset = threat_basis_data;
    threat_basis_data_subset.basis_subset = identifiable_basis;
    threat_basis_data_subset.basis_parameters.mean = ...
        threat_basis_data.basis_parameters.mean(:, identifiable_basis);
    threat_basis_data_subset.n_threat_parameters = numel(identifiable_basis);
    
    [tpe_mean, tpe_covar] = threat_estimate_rbf(sensor_locations, grid_world, ...
        threat_basis_data_subset, threat_parameters_true(identifiable_basis), sensor_noise);
    
	threat_value_estimate.mean	= calc_threat_rbf(threat_basis_data_subset, tpe_mean, grid_world.coordinates);       
	
	%===== Path planning
	transition_cost_estimate = grid_world.adjacency;
	for m1 = 1:grid_world.n_grid_points
		nhbrs	= find(grid_world.adjacency(m1, :));
		transition_cost_estimate(m1, nhbrs) = max( 1E-3, (threat_value_estimate.mean(:, nhbrs)) ); % intended to exclude negative vals? but erases connectivity for neg vals
	end
	search_data.adjacency_matrix= transition_cost_estimate;

% 	[vertex_data, ~]= astar(search_data);
    [vertex_data, ~]	= astar_hfunc(search_data); % Uses heurstic function
    ghost_path      = opt_path; % Save previous path for ghost visual 
	opt_path		= greedy_trace(search_data.v_start, search_data.v_goal, vertex_data);
    % update min_discrete_path_variance for possibly different length path
    min_discrete_path_variance = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar_full);

	%===== Key numbers
	expected_cost_estimate	= sum( max(0, threat_value_estimate.mean(opt_path(2:end))) );
	incurred_cost			= sum(threat_value_true(opt_path(2:end)));
	
    incurred_rec(m00+1)  = incurred_cost;
    expected_rec(m00+1)  = expected_cost_estimate;
    trace_cov_rec(m00+1) = trace(tpe_covar);
    path_cost_var_rec(m00+1) = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar);


    m00 = m00 + 1;
end
bin_sensor_locations(sensor_locations,m00 + 1) = 1;
norm_sensor_diff(m00) = norm(bin_sensor_locations(:,m00 + 1) - bin_sensor_locations(:,m00),1);
total_measurements = length(sensor_locations);
conv_iters = m00;

path_cost_incurred = incurred_rec(end);
path_cost_expected = expected_rec(end);

path_cost_var = path_cost_variance_rbf(opt_path, threat_basis_data, threat_basis_data_subset,...
                                       grid_world, tpe_covar);                                   

end
