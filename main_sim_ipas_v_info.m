%{
Copyright 2018 Benjamin Cooper

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
%}
% Comparison of rbf IPAS vs. Info approach
% This version contains modifications such that:
%   - IPAS sensor reconf uses minimum necessary sensors 
%   - Info allowed to use same minimum measurements as IPAS
%   - Info also only evaluated using estimable basis (basis in range of any
%   of the sensors placed), basically the same as IPAS reconstruction
% The above levels the playing field for comparing IPAS v INFO

% Also includes min_discrete_variance which is the minimum variance of the
% discretization. If you take a measurement at every grid point (H has
% dimension N_G^2 x N_P), then this is minimum achievable variance. We want
% to be within some k^2 multiple of that variance:
% Var(\theta_{subset})_{IPAS} < k^2*Var_{min-discretization}

% Added on from Cumulative Least Squares
% Cumulative in the sense that rather than
%   a) throw away old measurements at each iteration or
%   b) use recursive LSE
% this version just stacks measurements into larger and larger vectors to
% be evaluated. This should be computationally bad way to do sequential
% estimation, but I'm looking at LSE for underdetermined systems.
% So mainly, experimenting to see what's going on with estimate and
% error covariances.

clc; clear all
% Specify primary parameters for study

n_sensors = 10;
n_params = 6; % Total parameters will be n_params^2
n_grid_dim = 20;
SEED = true;
DRAW_ITERATIONS = true; % Set true if you want to draw all iters of setup

%% TABULA RASA - the clean slate club
% clear variables; close all; clc
addpath('./routine/') % Add greedy sensor placement algorithms
% From J. Ranieri "Near Optimal Sensor Placement for Linear Inverse
% Problems", https://github.com/jranieri/OptimalSensorPlacement

% Specify seed for random field generation
% To give consistent field and true_cost
if SEED
%     rng(666);
    rng(12345678)
end

Ns = n_sensors;
Np = n_params;

folder_name	= ['Data_Pics/RBF_v_INFO_MinBasis_Np',num2str(Np^2),'_Ns',num2str(n_sensors),'_',...
    datestr(clock,'yyyymmdd_HHMMSS')];
% folder_name	= ['Data_Pics/RBF_Cumulative_Xcorr_Np',num2str(Np^2),'_Ns',num2str(n_sensors),'_',...
%     datestr(clock,'yyyymmdd_HHMMSS')];
giffilename = 'min-basis-animation';
mkdir(folder_name)
%% GRIDWORLD
%----- Assume square workspace
half_wksp_size = 1;
grid_world.half_wksp_size = half_wksp_size; 
grid_world.n_grid_row	= n_grid_dim; % Perfect square
grid_world.n_grid_points= grid_world.n_grid_row^2;
grid_world.spacing		= 2*half_wksp_size / (grid_world.n_grid_row - 1);

grid_world.coordinates	= zeros(2, grid_world.n_grid_points);
for m1 = 0:(grid_world.n_grid_points - 1)	
	grid_world.coordinates(:, m1 + 1) = [...
		-half_wksp_size + (mod(m1, grid_world.n_grid_row))*grid_world.spacing; ...
		-half_wksp_size + floor(m1/grid_world.n_grid_row)*grid_world.spacing];
end

%----- Setup adjacency matrix
n_edges		= 0;
n_exp_edges	= grid_world.n_grid_points*4;
edge_list	= zeros(n_exp_edges, 3);
for m1 = 1:grid_world.n_grid_points
	if (m1 + 1 <= grid_world.n_grid_points) && (mod(m1, grid_world.n_grid_row) ~= 0)
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [m1 (m1 + 1) 1];
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [(m1 + 1) m1 1];
	end

	if (m1 + grid_world.n_grid_row) <= grid_world.n_grid_points
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [m1 (m1 + grid_world.n_grid_row) 1];
		n_edges				= n_edges + 1;
		edge_list(n_edges, :) = [(m1 + grid_world.n_grid_row) m1 1];
	end
end
grid_world.adjacency = sparse(edge_list(1:n_edges,1), edge_list(1:n_edges,2), edge_list(1:n_edges,3));

v_start	  = 1;
loc_start = grid_world.coordinates(:,v_start);
v_goal	  = grid_world.n_grid_points;
loc_goal  = grid_world.coordinates(:,v_goal);
%% SENSORS

% -----------------------------------------------------------------------
% ----------------------- NO OF SENSORS ---------------------------------
n_sensor = Ns;																% Less than #grid points
% -----------------------------------------------------------------------
noise_sigma = 0.1;
sensor_noise.n_sensor = Ns;
sensor_noise.noise_sigma = noise_sigma;
sensor_noise.R		= sparse( noise_sigma.^2*diag(ones(1, n_sensor)) );
sensor_noise.Rinv	= sparse( (1/(noise_sigma.^2))*diag(ones(1, n_sensor)) );

%===== SENSOR PLACEMENT EXPERIMENTS
sensor_locations = zeros(1, n_sensor);
%----- Arbitary placement
for m1 = 1:n_sensor
	this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
	while any(this_location == sensor_locations(1:m1))
		this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
	end
	sensor_locations(m1) = this_location;
end
% Or start with no sensors
sensor_locations = [1];


%% THREAT FIELD GENERATION
%{
Spatial basis functions are symmetric 2D Gaussians
	\sqrt{\sigma_\Psi} is the "one sigma" width
	(\mu_{x,i}, \mu_{y,i}) is the "mean"(center)
	\mathsf{x} = [x,y]-[\mu_{x,i},\mu_{y,i}]

	Psi_p(x,y) = exp( \mathsf{x} * diag((1/\sigma), (1/\sigma)) * \mathsf{x}')

	For now, basis functions have centers uniformly spaced in an equal
	number of  rows and columns.
%}
% ------------------------------------------------------------------------
% ------------------- NO OF THREAT PARAMETERS ----------------------------
% ------------------------------------------------------------------------
% Visually code a square map and trasnform to appropriate format ---------
threat_parameters_true_map = [1 1 1 1 0 0; ...  % How the threats look visually
                              1 0.28 0 0 0 0; ... % but, need to flip, transpose 
                              1 0 5 5 0 1; ... % and reshape to get into correct
                              1 0 5 5 0 1; ...  % order for location on field in
                              0 0 0 0 .2 1; ...  % threat_basis_data below.
                              0 0 1 1 1 1];
% threat_parameters_true_map = [5 5 1 1 1 0 0 0 0 0; ...  % How the threats look visually
%                               5 1 1 1 4 2 0 0 8 0; ... % but, need to flip, transpose
%                               1 1 1 5 4 0 0 0 0 0; ... % and reshape to get into correct
%                               5 5 5 1 0 0 0 0 0 0; ...  % order for location on field in
%                               5 5 5 0 0 8 0 0 3 3; ...  % threat_basis_data below.
%                               5 5 0 0 0 3 3 3 3 3; ...  
%                               5 5 0 0 3 3 3 3 3 3; ...  
%                               5 5 0 0 3 3 3 3 0 0; ...  
%                               5 5 0 0 3 3 3 0 0 0; ...  
%                               0 0 0 1 2 2 0 0 0 0];
if exist('threat_parameters_true_map','var')
    threat_parameters_true = reshape(flipud(threat_parameters_true_map)',[numel(threat_parameters_true_map), 1]);
    threat_basis_data.n_threat_parameters = numel(threat_parameters_true_map);
    threat_basis_data.threat_map = threat_parameters_true_map;
else
    % Perfect squares for now, for equal #rows and columns
    threat_basis_data.n_threat_parameters	= Np^2;								
    % Range of rand threat params (a,b): r = a + (b - a).*rand(N,1);
    min_threat = 0;
    max_threat = 10;
    threat_parameters_true = min_threat + (max_threat - min_threat).*...
        rand(threat_basis_data.n_threat_parameters,1);
end
% ------------------------------------------------------------------------
n_center_rows	= sqrt(threat_basis_data.n_threat_parameters);
% center_spacing	= 2*half_wksp_size / (n_center_rows + 1);
center_spacing	= 2*half_wksp_size / (n_center_rows - 1);
threat_basis_data.basis_parameters.mean = zeros(2, threat_basis_data.n_threat_parameters);
for m1 = 1:n_center_rows
	for m2 = 1:n_center_rows
		threat_basis_data.basis_parameters.mean(:, (m2 - 1)*n_center_rows + m1) = ...
			[-half_wksp_size + (m1 - 1)*center_spacing; -half_wksp_size + (m2 - 1)*center_spacing];
	end
end
    
% basis variance chosen such that the range of influence of each basis just
% touch at the diagonal between two basis centers. Basis centers are 1.4 *
% center_spacing distance along diagonal, therefore meet at 0.7 * center
% Range of influce (x_rng) is governed by threat value and noise level. It
% has similar meaning/use to a signal-to-noise ratio (SNR)
x_rng = 1/2*sqrt(2) * center_spacing;
theta_guess = 5; % theta_guess should be a value in the median range of parameter values
basis_variance = x_rng^2 / (2 * log(theta_guess / noise_sigma));
threat_basis_data.basis_parameters.var	= basis_variance;		% This is \sigma^2_\Psi
threat_basis_data.offset				= 1;

%% PATH PLANNING WITH TRUE THREAT (FOR COMPARISON)

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
search_data.heuristic		= zeros(grid_world.n_grid_points, 1);
search_data.grid_world      = grid_world;

[vertex_data_true, ~]	= astar(search_data);
opt_path_true = greedy_trace(search_data.v_start, search_data.v_goal, vertex_data_true);
% fprintf('Cost of true optimal path:\t%f\n', vertex_data_true(search_data.v_goal).d)

path_cost_true       		= sum(threat_value_true(opt_path_true(2:end)));

%===== 2D VISUALATION OF TRUE FIELD AND TRUE OPTIMAL PATH
%-----Threat field
n_threat_plot	= 101;
[x_plot, y_plot]= meshgrid( linspace(-half_wksp_size, half_wksp_size, n_threat_plot), ...
	linspace(-half_wksp_size, half_wksp_size, n_threat_plot));
x_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
y_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
threat_plot		= threat_basis_data.offset*ones(n_threat_plot);

for m1 = 1:threat_basis_data.n_threat_parameters
	threat_plot = threat_plot + threat_parameters_true(m1) * exp(...
		(-1 / (2 * threat_basis_data.basis_parameters.var)) * ...
		((x_plot - threat_basis_data.basis_parameters.mean(1, m1)).^2 + ...
		(y_plot - threat_basis_data.basis_parameters.mean(2, m1)).^2) );
end
figure()
ax = gca;
set(gca,'FontSize',16)
% surfc(x_plot,y_plot,threat_plot,'LineStyle','none')
imagesc(ax, x_plot_vec, y_plot_vec, threat_plot)
set(ax, 'YDir','normal')
title(['True optimal path cost = ',num2str(path_cost_true)])
colorbar; view(2); axis equal; hold on;
xlim([-half_wksp_size half_wksp_size]);
ylim([-half_wksp_size half_wksp_size]);
currColorLims = [min(threat_plot(:)) max(threat_plot(:))];

%----- Gridworld
plot(grid_world.coordinates(1, :), grid_world.coordinates(2, :), ...
	 'k.', 'MarkerSize',20);

%----- Path
plot(grid_world.coordinates(1, opt_path_true), grid_world.coordinates(2, opt_path_true), ...
	 'wo', 'MarkerSize',10','LineWidth',3);
plot_name = 'true-map-path';
set(gcf, 'InvertHardcopy', 'off')
print('-dpng', fullfile(folder_name,plot_name))
close

%% PATH PLANNING WITH ESTIMATED THREAT IN LOOP WITH SENSOR RECONFIGURATION

 % Iterations should not be larger than total grid points divided by
 % available sensors (reasonably). May be necessary if required Var high.
iters = ceil(grid_world.n_grid_points / n_sensor);
incurred_rec = zeros(1,1);
estimated_rec = zeros(1,1);
trace_cov_rec = zeros(1,1);
path_cost_var_rec = zeros(1,1);

% Binary sensor vector [Ng^2 x iters]
bin_sensor_locations = zeros(grid_world.n_grid_points,iters);
bin_sensor_locations(sensor_locations,1) = 1;
n_sensor_used(1) = length(sensor_locations);

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

[vertex_data, ~]= astar(search_data);
opt_path		= greedy_trace(search_data.v_start, search_data.v_goal, vertex_data);
ghost_path = opt_path; % ghost path is used to visualize previous iteration path
                       % which should cover current sensor placement

%===== Key numbers
estimated_cost_estimate	= sum( max(0, threat_value_estimate.mean(opt_path(2:end))) );
incurred_cost			= sum(threat_value_true(opt_path(2:end)));
% fprintf('\n===== ITERATION 0 =====\n')
% fprintf('Worst-case cost of optimal path with estimated threat: \t%f\n', vertex_data(search_data.v_goal).d)
% fprintf('estimated cost of optimal path with estimated threat: \t%f\n', estimated_cost_estimate)
% fprintf('Incurred cost of optimal path with estimated threat: \t%f\n', incurred_cost)

incurred_rec(1)= incurred_cost;
estimated_rec(1)= estimated_cost_estimate;
path_cost_var_rec(1) = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar);

%==== Drawing first iteration
if DRAW_ITERATIONS
    h = figure();
    axis tight manual % this ensures that getframe() returns a consistent size
    ax = gca;
% 	cla;
	n_threat_plot	= 101;
	[x_plot, y_plot]= meshgrid( linspace(-half_wksp_size, half_wksp_size, n_threat_plot), ...
		linspace(-half_wksp_size, half_wksp_size, n_threat_plot));
    x_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
    y_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
	threat_plot		= threat_basis_data_subset.offset*ones(n_threat_plot);

	for m1 = 1:threat_basis_data_subset.n_threat_parameters
		threat_plot = threat_plot + tpe_mean(m1) * exp(...
			(-1 / (2 * threat_basis_data_subset.basis_parameters.var)) * ...
			((x_plot - threat_basis_data_subset.basis_parameters.mean(1, m1)).^2 + ...
			(y_plot - threat_basis_data_subset.basis_parameters.mean(2, m1)).^2) );
    end
    
    threat_plot_measured = NaN*threat_plot; % Start with unknown background
    for m2 = 1:threat_basis_data_subset.n_threat_parameters
        % Fill in threat map values for regions measured
        x_loc = threat_basis_data_subset.basis_parameters.mean(1, m2);
        y_loc = threat_basis_data_subset.basis_parameters.mean(2, m2);
        max_range = 3*sqrt(threat_basis_data_subset.basis_parameters.var);
        threat_plot_measured((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2) = ...
            threat_plot((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2);
    end
% 	surfc(x_plot,y_plot,threat_plot,'LineStyle','none')
    hp = pcolor(ax, x_plot, y_plot, threat_plot_measured);
    set(hp, 'edgecolor', 'none')
    set(h,'Color',[1 1 1]);
    set(ax,'color',[0.5 0.5 0.5]) % background gray for unknown/unmeasured map
%     set(ax, 'YDir','normal')
    title(['Iter. ',num2str(1),',   Incurred Cost: ',num2str(incurred_cost)],'FontSize',16)
    caxis(currColorLims);
    colorbar; view(2); axis equal; hold on;
    xlim([-half_wksp_size half_wksp_size]);
    ylim([-half_wksp_size half_wksp_size]);

	%--- Gridworld
	plot(grid_world.coordinates(1, :), grid_world.coordinates(2, :), 'k.', 'MarkerSize',20);

	%--- Sensors
	plot(grid_world.coordinates(1, sensor_locations), grid_world.coordinates(2, sensor_locations),...
         'ko', 'MarkerSize',20,'LineWidth',3);
    
%     %--- Basis Locations
%     plot(threat_basis_data.basis_parameters.mean(1,:),...
%          threat_basis_data.basis_parameters.mean(2,:), 'ro', 'MarkerSize',20, 'LineWidth',2);
%     
%     %--- Identifiable Basis Locations
%     draw_basis_range(ax, threat_basis_data, identifiable_basis, 'g')
    
	%--- Path
	plot(grid_world.coordinates(1, opt_path), grid_world.coordinates(2, opt_path),...
         'wo', 'MarkerSize',10','LineWidth',3);
    xlim([-half_wksp_size half_wksp_size]);
    ylim([-half_wksp_size half_wksp_size]);
	drawnow()
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,fullfile(folder_name,giffilename),'gif','Loopcount',inf);
    
    plot_name = ['ipas-iter-' num2str(1)];
    set(h, 'InvertHardcopy', 'off')
	print('-dpng', fullfile(folder_name,plot_name))
    close
end

%=======================================================================
%========= Iterate with Cumulative LS from reconfigured sensors ========
%=======================================================================
% tpe_covar_full is the minimum covariance matrix possible if every grid
% point was measured. Serves as the lower bound on paremeter estimation and
% resulting path cost variance.
[tpe_mean_full, tpe_covar_full] = threat_estimate_rbf(1:grid_world.n_grid_points,...
    grid_world, threat_basis_data, threat_parameters_true, sensor_noise);

min_discrete_path_variance = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar_full);
k = 5; % consant which limits the sub-optimality of IPAS path variance
k_sqr_discrete_path_variance = k^2*min_discrete_path_variance;

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
    bin_sensor_locations(new_locations,m00 + 1) = 1;
    norm_sensor_diff(m00) = norm(bin_sensor_locations(:,m00+1) - bin_sensor_locations(:,m00),1);
    n_sensor_used(m00+1) = length(new_locations);
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

	[vertex_data, ~]= astar(search_data);
    ghost_path      = opt_path; % Save previous path for ghost visual 
	opt_path		= greedy_trace(search_data.v_start, search_data.v_goal, vertex_data);
    % update min_discrete_path_variance for possibly different length path
    min_discrete_path_variance = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar_full);
    k_sqr_discrete_path_variance = k^2*min_discrete_path_variance;

	%===== Key numbers
	estimated_cost_estimate	= sum( max(0, threat_value_estimate.mean(opt_path(2:end))) );
	incurred_cost			= sum(threat_value_true(opt_path(2:end)));
% 	fprintf('\n===== ITERATION %i =====\n', m00)
% 	fprintf('Worst-case cost of optimal path with estimated threat: \t%f\n', vertex_data(search_data.v_goal).d)
% 	fprintf('estimated cost of optimal path with estimated threat: \t%f\n', estimated_cost_estimate)
% 	fprintf('Incurred cost of optimal path with estimated threat: \t%f\n', incurred_cost)
% 	fprintf('Norm of binary sensor vector difference: \t%f\n', norm_sensor_diff(end))
%   fprintf('Path Cost variance: \t%f\n', path_cost_var_rec(end))
	
    incurred_rec(m00+1)  = incurred_cost;
    estimated_rec(m00+1)  = estimated_cost_estimate;
    trace_cov_rec(m00+1) = trace(tpe_covar);
    path_cost_var_rec(m00+1) = path_cost_variance_rbf(opt_path, threat_basis_data,...
                       threat_basis_data_subset, grid_world, tpe_covar);
    
    if DRAW_ITERATIONS
        %===== Drawing
        %---Threat field
        h = figure();
        axis tight manual % this ensures that getframe() returns a consistent size
        ax = gca;
        % 	cla;
        n_threat_plot	= 101;
        [x_plot, y_plot]= meshgrid( linspace(-half_wksp_size, half_wksp_size, n_threat_plot), ...
            linspace(-half_wksp_size, half_wksp_size, n_threat_plot));
        x_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
        y_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
        threat_plot		= threat_basis_data_subset.offset*ones(n_threat_plot);
        
        for m1 = 1:threat_basis_data_subset.n_threat_parameters
            threat_plot = threat_plot + tpe_mean(m1) * exp(...
                (-1 / (2 * threat_basis_data_subset.basis_parameters.var)) * ...
                ((x_plot - threat_basis_data_subset.basis_parameters.mean(1, m1)).^2 + ...
                (y_plot - threat_basis_data_subset.basis_parameters.mean(2, m1)).^2) );
        end
        
        threat_plot_measured = NaN*threat_plot; % Start with unknown background
        for m2 = 1:threat_basis_data_subset.n_threat_parameters
            % Fill in threat map values for regions measured
            x_loc = threat_basis_data_subset.basis_parameters.mean(1, m2);
            y_loc = threat_basis_data_subset.basis_parameters.mean(2, m2);
            max_range = 3*sqrt(threat_basis_data_subset.basis_parameters.var);
            threat_plot_measured((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2) = ...
                threat_plot((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2);
        end
        % 	surfc(x_plot,y_plot,threat_plot,'LineStyle','none')
        hp = pcolor(ax, x_plot, y_plot, threat_plot_measured);
        set(hp, 'edgecolor', 'none')
        set(ax,'color',[0.5 0.5 0.5]) % background gray for unknown/unmeasured map
        %     set(ax, 'YDir','normal')
        title(['Iter. ',num2str(m00+1),',   Incurred Cost: ',num2str(incurred_cost)],'FontSize',16)
        caxis(currColorLims);
        colorbar; view(2); axis equal; hold on;
        xlim([-half_wksp_size half_wksp_size]);
        ylim([-half_wksp_size half_wksp_size]);
        
        %--- Gridworld
        plot(grid_world.coordinates(1, :), grid_world.coordinates(2, :), ...
             'k.', 'MarkerSize',20);
        
        %--- Sensors
        plot(grid_world.coordinates(1, new_locations), grid_world.coordinates(2, new_locations), ...
             'ko', 'MarkerSize',20,'LineWidth',3);
         
%          %--- Vehicle Position
%         plot(grid_world.coordinates(1, v_current), grid_world.coordinates(2, v_current), ...
%              'go', 'MarkerSize',25,'LineWidth',4);
         
         %--- Basis Locations
%          plot(threat_basis_data.basis_parameters.mean(1,:),...
%              threat_basis_data.basis_parameters.mean(2,:), 'ro', 'MarkerSize',20, 'LineWidth',2);
%          
%          %--- Identifiable Basis Locations
%          draw_basis_range(ax, threat_basis_data, identifiable_basis, 'g')
        
        %--- Ghost Path
%         plot(grid_world.coordinates(1, ghost_path), grid_world.coordinates(2, ghost_path),...
%               'o','Color',[0.8 0.8 0.8], 'MarkerSize',10','LineWidth',3);
        
        %--- Path
        plot(grid_world.coordinates(1, opt_path), grid_world.coordinates(2, opt_path),...
             'wo', 'MarkerSize',10','LineWidth',3);
        
        xlim([-half_wksp_size half_wksp_size]);
        ylim([-half_wksp_size half_wksp_size]);
        drawnow()
        
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,fullfile(folder_name,giffilename),'gif','WriteMode','append');
        
        plot_name = ['ipas-iter-' num2str(m00+1)];
        set(h, 'InvertHardcopy', 'off')
        print('-dpng', fullfile(folder_name,plot_name))
        close
    end

    m00 = m00 + 1;
end
% bin_sensor_locations(sensor_locations,m00 + 1) = 1;
% norm_sensor_diff(m00) = norm(bin_sensor_locations(:,m00 + 1) - bin_sensor_locations(:,m00),1);

%========================================================================
%=========== Greedy Sensor Placement Using Frame Potential ===========
% From J. Ranieri "Near Optimal Sensor Placement for Linear Inverse
% Problems", https://github.com/jranieri/OptimalSensorPlacement
%========================================================================
tic % Start timer for FP path plan
Psi = calc_rbf_value(threat_basis_data, grid_world.coordinates);
% let INFO approach use as many measurements as IPAS
sensor_locations_FP = SP_greedyFP(Psi,length(sensor_locations)); 

%=============     Batch Estimation of Threat Parameters     ============
identifiable_basis_FP = identify_basis_close2sensor(sensor_locations_FP, grid_world, threat_basis_data);
threat_basis_data_subset_FP = threat_basis_data;
threat_basis_data_subset_FP.basis_subset = identifiable_basis_FP;
threat_basis_data_subset_FP.basis_parameters.mean = ...
    threat_basis_data.basis_parameters.mean(:, identifiable_basis_FP);
threat_basis_data_subset_FP.n_threat_parameters = numel(identifiable_basis_FP);
[tpe_mean_FP, tpe_covar_FP] = threat_estimate_rbf(sensor_locations_FP, grid_world, ...
    threat_basis_data_subset_FP, threat_parameters_true(identifiable_basis_FP), sensor_noise);

% Reconstruct Field using only parameters with confidence
threat_value_estimate_FP.mean	= calc_threat_rbf(threat_basis_data_subset_FP, tpe_mean_FP, grid_world.coordinates);

%===== Path planning
transition_cost_estimate_FP = grid_world.adjacency;
for m1 = 1:grid_world.n_grid_points
    nhbrs	= find(grid_world.adjacency(m1, :));
    % intended to exclude negative vals? but erases connectivity for neg vals
    transition_cost_estimate_FP(m1, nhbrs) = max( 1E-3, (threat_value_estimate_FP.mean(:, nhbrs)) );
end
search_data.adjacency_matrix = transition_cost_estimate_FP;

[vertex_data_FP, ~]= astar(search_data);
opt_path_FP		= greedy_trace(search_data.v_start, search_data.v_goal, vertex_data_FP);

info_estimated_cost_estimate	= sum( max(0, threat_value_estimate_FP.mean(opt_path_FP(2:end))) );
info_incurred_cost			= sum(threat_value_true(opt_path_FP(2:end)));

%% Draw INFO plot and path
if DRAW_ITERATIONS
    %===== Drawing
    %---Threat field
    h = figure();
    axis tight manual % this ensures that getframe() returns a consistent size
    ax = gca;
    % 	cla;
    n_threat_plot	= 101;
    [x_plot, y_plot]= meshgrid( linspace(-half_wksp_size, half_wksp_size, n_threat_plot), ...
        linspace(-half_wksp_size, half_wksp_size, n_threat_plot));
    x_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
    y_plot_vec = linspace(-half_wksp_size, half_wksp_size, n_threat_plot);
    threat_plot		= threat_basis_data_subset_FP.offset*ones(n_threat_plot);
    
    for m1 = 1:threat_basis_data_subset_FP.n_threat_parameters
        threat_plot = threat_plot + tpe_mean_FP(m1) * exp(...
            (-1 / (2 * threat_basis_data_subset_FP.basis_parameters.var)) * ...
            ((x_plot - threat_basis_data_subset_FP.basis_parameters.mean(1, m1)).^2 + ...
            (y_plot - threat_basis_data_subset_FP.basis_parameters.mean(2, m1)).^2) );
    end
    
    threat_plot_measured = NaN*threat_plot; % Start with unknown background
    for m2 = 1:threat_basis_data_subset_FP.n_threat_parameters
        % Fill in threat map values for regions measured
        x_loc = threat_basis_data_subset_FP.basis_parameters.mean(1, m2);
        y_loc = threat_basis_data_subset_FP.basis_parameters.mean(2, m2);
        max_range = 3*sqrt(threat_basis_data_subset_FP.basis_parameters.var);
        threat_plot_measured((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2) = ...
            threat_plot((x_plot - x_loc).^2 + (y_plot - y_loc).^2 < max_range.^2);
    end
    % 	surfc(x_plot,y_plot,threat_plot,'LineStyle','none')
    hp = pcolor(ax, x_plot, y_plot, threat_plot_measured);
    set(hp, 'edgecolor', 'none')
    set(h,'Color',[1 1 1]);
    set(ax,'color',[0.5 0.5 0.5]) % background gray for unknown/unmeasured map
    %     set(ax, 'YDir','normal')
    title(['Info-driven placement, Incurred Cost: ',num2str(info_incurred_cost)],'FontSize',16)
    caxis(currColorLims);
    colorbar; view(2); axis equal; hold on;
    xlim([-half_wksp_size half_wksp_size]);
    ylim([-half_wksp_size half_wksp_size]);
    
    %--- Gridworld
    plot(grid_world.coordinates(1, :), grid_world.coordinates(2, :), ...
        'k.', 'MarkerSize',20);
    
    %--- Sensors
    plot(grid_world.coordinates(1, sensor_locations_FP), grid_world.coordinates(2, sensor_locations_FP), ...
        'ko', 'MarkerSize',20,'LineWidth',3);
    
    %          %--- Vehicle Position
    %         plot(grid_world.coordinates(1, v_current), grid_world.coordinates(2, v_current), ...
    %              'go', 'MarkerSize',25,'LineWidth',4);
    
    %--- Basis Locations
    %          plot(threat_basis_data.basis_parameters.mean(1,:),...
    %              threat_basis_data.basis_parameters.mean(2,:), 'ro', 'MarkerSize',20, 'LineWidth',2);
    %
    %          %--- Identifiable Basis Locations
    %          draw_basis_range(ax, threat_basis_data, identifiable_basis, 'g')
    
    %--- Ghost Path
    %         plot(grid_world.coordinates(1, ghost_path), grid_world.coordinates(2, ghost_path),...
    %               'o','Color',[0.8 0.8 0.8], 'MarkerSize',10','LineWidth',3);
    
    %--- Path
    plot(grid_world.coordinates(1, opt_path_FP), grid_world.coordinates(2, opt_path_FP),...
        'wo', 'MarkerSize',10','LineWidth',3);
    
    xlim([-half_wksp_size half_wksp_size]);
    ylim([-half_wksp_size half_wksp_size]);
    drawnow()
    
    plot_name = 'info-map-path';
    set(h, 'InvertHardcopy', 'off')
    print('-dpng', fullfile(folder_name,plot_name))
    close
end


cost_diff = abs(info_incurred_cost - path_cost_true);
percent_diff = cost_diff/path_cost_true*100;
info_theta_var = trace(tpe_covar_FP);
info_path_cost_var = path_cost_variance_rbf(opt_path_FP, threat_basis_data,...
                        threat_basis_data_subset_FP, grid_world, tpe_covar_FP);
%% Plotting
h = figure('position', [50, 50, 650, 800]);

subplot(311)
semilogy(1:m00, info_path_cost_var*ones(1,m00), '--r', 1:m00, path_cost_var_rec, 'b','LineWidth',2)
set(gca,'FontSize',14)
title({['N_P = ' num2str(n_params^2),', grid points = ', num2str(n_grid_dim^2)],...
       ['N_S = ',num2str(n_sensor), ', # measurements = ',num2str(length(sensor_locations))]}) 
% xlabel('Number of iterations')
set(gca,'xtick',1:m00)
YTick = logspace(round(log10(min_discrete_path_variance)),...
                         round(log10(max(path_cost_var_rec))),3);
YTickLabels = cellstr(num2str(round(log10(YTick(:))), '10^%d'));
set(gca,'ytick',YTick)
set(gca,'yticklabels',YTickLabels)
ylabel('Var(J(v))')
ylim([0.1*YTick(1) inf])
legend({'Info - single placement', 'IPAS - sequential placement'}, 'Location','Best', 'FontSize', 14, 'Box', 'off')

subplot(312)
plot(1:m00, incurred_rec, ':b', 1:m00, info_incurred_cost*ones(1,m00), '--r', ...
    1:m00, path_cost_true*ones(1,m00), 'g', 'LineWidth', 2)
set(gca,'FontSize',14)
% xlabel('Number of iterations')
set(gca,'xtick',1:m00)
ylabel('Path Cost');
ylim([path_cost_true*0.9 inf])
legend({'IPAS', 'Info' ,'True optimal'}, 'Location','Best',...
    'FontSize', 14, 'Box', 'off','Orientation','horizontal')

subplot(313)
plot(n_sensor_used,'b','LineWidth',2)
set(gca,'FontSize',14)
% title('Number of sensors used')
xlabel('Iteration')
set(gca,'xtick',1:m00)
ylabel('No. of sensors used')

plot_name = 'info-v-ipas-plot';
% set(h, 'InvertHardcopy', 'off')
print(fullfile(folder_name,plot_name),'-dpng','-r300')

% Plotting only for IPAS, no info stuff ----------------------------------
h = figure('position', [50, 50, 650, 800]);

subplot(311)
semilogy(1:m00, path_cost_var_rec, 'b', 1:m00, k_sqr_discrete_path_variance*ones(1,m00),...
    ':r', 1:m00, min_discrete_path_variance*ones(1,m00), '--g')
set(gca,'FontSize',14)
title({['N_P = ' num2str(n_params^2),', grid points = ', num2str(n_grid_dim^2)],...
       ['N_S = ',num2str(n_sensor), ', # measurements = ',num2str(length(sensor_locations))]}) 
% xlabel('Number of iterations')
set(gca,'xtick',1:m00)
YTick = logspace(round(log10(min_discrete_path_variance)),...
                         round(log10(max(path_cost_var_rec))),3);
YTickLabels = cellstr(num2str(round(log10(YTick(:))), '10^%d'));
set(gca,'ytick',YTick)
set(gca,'yticklabels',YTickLabels)
ylabel('Var(J(v))')
ylim([0.1*YTick(1) inf])
legend({'IPAS Var(J(v))', 'k^2 Var(J(v))', 'Grid Var(J(v))'},'Location', 'Best', 'FontSize', 12, 'Box', 'off')

subplot(312)
plot(1:m00, incurred_rec, ':b', 1:m00, estimated_rec, '--r' ,1:m00, path_cost_true*ones(1,m00), 'g')
set(gca,'FontSize',14)
% xlabel('Number of iterations')
set(gca,'xtick',1:m00)
ylabel('Path Cost');
legend({'IPAS incurred cost', 'IPAS estimated cost' ,'True optimal cost'}, 'Location','Best', 'FontSize', 12, 'Box', 'off')

subplot(313)
plot(n_sensor_used,'b','LineWidth',2)
set(gca,'FontSize',14)
% title('Number of sensors used')
xlabel('Iteration')
set(gca,'xtick',1:m00)
ylabel('sensors used')

plot_name = 'variance-cost-plot';
% set(h, 'InvertHardcopy', 'off')
print(fullfile(folder_name,plot_name),'-dpng','-r300')

%% Save problem data, variables to reproduce later
ipas_vars.k = k;
ipas_vars.theta_guess = theta_guess;
ipas_vars.v_start = v_start;
ipas_vars.v_goal = v_goal;

ipas_result.n_sensors_used = n_sensor_used;
ipas_result.sensor_locations = sensor_locations;

problem_data_name = 'problem_data.mat';
save(fullfile(folder_name,problem_data_name),'grid_world','sensor_noise',...
    'threat_basis_data', 'ipas_vars', 'ipas_result')

