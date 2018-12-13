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
% Generate data from a set of simulations to compare IPAS with various
% information driven sensor placement approaches.
% Sensor Placement methods:
%   1) IPAS
%   2) Random
%   3) greedy MI
%   4) greedy FP
%   5) greedy MSE
% Data to collect
%   1) Np: Number of parameters
%   2) Ns: Number of sensors
%   3) true_cost: True path cost
%   4) inc_cost: Incurred path cost
%   5) exp_cost: Expected path cost
%   6) iters: # of iterations to converge
%   7) comp_time: End-to-End computation time
%   8) sensor_movement: Total sensors moved (norm_sensor_diff)
%   9) var_path_cost: The final variance of path cost
%  10) grid_size: Total number of nodes in grid (n_grid_row^2)
%  11) noise_sig: the noise level sigma applied to random 

clc; clear all;

dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
fprintf('Beginning simulation runs at: %s\n', dt)

PRINT_THINGS = false;

addpath('./routine/') % Add greedy sensor placement algorithms
% From J. Ranieri "Near Optimal Sensor Placement for Linear Inverse
% Problems", https://github.com/jranieri/OptimalSensorPlacement

%% DEFINE INDEPENDENT PARAMETERS: Np, Ns, grid_size, noise

% Define a sensor range and noise
min_sensors = 1;
max_sensors = 90;
sensor_list = min_sensors:max_sensors;

noise_sigma = 0.1;
sensor_noise.sigma = noise_sigma;
sensor_noise.R		= sparse( noise_sigma.^2*diag(ones(1, min_sensors)) );
sensor_noise.Rinv	= sparse( (1/(noise_sigma.^2))*diag(ones(1, min_sensors)) );

% Define parameter size
min_params = 2; % Perfect squares of parameters enforced in threat creation
max_params = 9; % eg. [1 2 3 4 5 6 7] => [1 4 9 16 25 36 49]
param_list = min_params:max_params;   % total number = n_params^2

% Define grid size
min_n_grid_row = 20;
max_n_grid_row = 20;
n_grid_row_step = 20;
n_grid_row_list = min_n_grid_row:n_grid_row_step:max_n_grid_row;  % total grid size = n_grid_row^2
% n_grid_row_list = [20 50 80 100 120]; % Or provide a list of grid sizes

fprintf('Grid sizes: ')
disp(n_grid_row_list)
fprintf('Np:         ')
disp(param_list.^2)
fprintf('Ns:         %i ... %i\n', sensor_list(1), sensor_list(end))


% Make generated fields repeated, so that path_cost_true is the same in
% this round of tests
SEED = false;

folder_name	= ['Sim_Data/IPAS-v-INFO-PARALLEL-',...
    'NG2_',num2str(n_grid_row_list(1)^2),'-',num2str(n_grid_row_list(end)^2),'_',...
    'Np_',num2str(param_list(1)^2),'-',num2str(param_list(end)^2),'_',...
    'Ns_',num2str(sensor_list(1)),'-',num2str(sensor_list(end)),'_',...
    datestr(clock,'yyyymmdd_HHMMSS')];
sim_data_file_name = 'parallel_sim_record.mat';


mkdir(folder_name)

%% DEFINE DATA STRUCTURE TO STORE ALL SIMULATIONS

sim_meta_data.grid_list = n_grid_row_list;
sim_meta_data.param_list = param_list;
sim_meta_data.sensor_list = sensor_list;

full_sims = [];
N_sim = 30; % Number of simulations per Ng-Np-Ns tuple
fprintf('No. of full simulations: %i\n',N_sim)

parfor_progress(N_sim); % Initialize progress bar
parfor sim_i = 1:N_sim
    sim_data_record = [];
%     fprintf('\nStart Full Simulation %i/%i \n', sim_i, N_sim);
    grid_num = 1; % struct index for curr_grid
    for curr_grid = n_grid_row_list  % ITERATE GRIDS
        % ------------------------------------------------------------------------
        % -------------------------- GRIDWORLD -----------------------------------
        % PICK SENSOR LOCATIONS FROM THIS DATA STRUCTURE.
        % ------------------------------------------------------------------------
        % This is the data structure that contains all available locations in the
        % environment.
        % grid_world.coordinates is a (2 x n_grid_row^2) matrix
        % a sensor at the first coordinate would be:
        % first_sensor_location = grid_world.coordinate(:, 1)
        % last_sensor_location = grid_world.coordinates(:, end)
        %----- Assume square workspace
        grid_world = []; % necessary for parfor
        half_wksp_size = 1;
        grid_world.n_grid_row	= curr_grid;
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
        
        % ------------------------------------------------------------------------
        Np_num = 1;   % struct index for current parameter Np
        for Np = param_list  % ITERATE No. OF PARAMETERS
            % ------------------------------------------------------------------------
            % ------------------- NO OF THREAT PARAMETERS ----------------------------
            threat_basis_data = []; % necessary for parfor
            threat_basis_data.n_threat_parameters	= Np^2;								% Perfect squares for now, for equal #rows and columns
            % ------------------------------------------------------------------------
            % Range of rand threat params (a,b): r = a + (b - a).*rand(N,1);
            min_threat = 0;
            max_threat = 10;
            threat_parameters_true = min_threat + (max_threat - min_threat).*...
                rand(threat_basis_data.n_threat_parameters,1);
            n_center_rows	= sqrt(threat_basis_data.n_threat_parameters);
            center_spacing	= 2*half_wksp_size / (n_center_rows - 1);
            threat_basis_data.basis_parameters.mean = zeros(2, threat_basis_data.n_threat_parameters);
            for m1 = 1:n_center_rows
                for m2 = 1:n_center_rows
                    threat_basis_data.basis_parameters.mean(:, (m2 - 1)*n_center_rows + m1) = ...
                        [-half_wksp_size + (m1 - 1)*center_spacing; -half_wksp_size + (m2 - 1)*center_spacing];
                end
            end
            threat_basis_data.threat_parameters_true = threat_parameters_true;
            % basis variance chosen such that the range of influence of each basis just
            % touch at the diagonal between two basis centers. Basis centers are 1.4 *
            % center_spacing distance along diagonal, therefore meet at 0.7 * center
            % Range of influence (x_rng) is governed by threat value and noise level. It
            % has similar meaning/use to a signal-to-noise ratio (SNR)
            % This also ensures the whole environment is covered by basis.
            x_rng = 1/2*sqrt(2) * center_spacing;
            theta_guess = 5; % theta_guess should be a value in the median range of parameter values
            basis_variance = x_rng^2 / (2 * log(theta_guess / noise_sigma));
            threat_basis_data.basis_parameters.var	= basis_variance;		% This is \sigma^2_\Psi
            threat_basis_data.offset				= 1;
            % -----------------------------------------------------------------------
            
            Ns_num = 1;   % struct index for current sensor Ns
            for Ns = sensor_list  % ITERATE No. OF SENSORS
                if PRINT_THINGS
                    fprintf('============================================\n')
                    fprintf('Running with:\n')
                    fprintf('Number of sensors = %d\n', Ns)
                    fprintf('Total Guass basis = %d\n', Np^2)
                    fprintf('Total grid size = %d\n\n', curr_grid^2)
                end
                
                % IPAS sensor placement approach ------------------------------
                if PRINT_THINGS; fprintf('Running IPAS...\n'); end
                tic
                [ipas_path_cost_incurred, ipas_path_cost_expected, ipas_path_cost_true,...
                    ipas_total_measurements, ipas_path_cost_var, ipas_iters] = ...
                    ipas_min_basis_rbf_blackbox(Ns, threat_basis_data, grid_world,...
                    sensor_noise);
                ipas_time = toc;
                if PRINT_THINGS
                    fprintf('No. of measurements = %i\n',ipas_total_measurements)
                    fprintf('Finished in %.3f seconds\n\n',ipas_time)
                end
                % -------------------------------------------------------------
                
                % For the following information approaches: use ipas_total_measurements
                Psi = calc_rbf_value(threat_basis_data, grid_world.coordinates);
                
                % greedy Frame Potential (FP) sensor placement ----------------
                if PRINT_THINGS; fprintf('Running greedy Frame Potential...\n'); end
                tic
                sensor_posns_FP = SP_greedyFP(Psi, ipas_total_measurements);
                [incurred_cost_FP, expected_cost_FP, path_cost_true_FP, path_cost_var_FP] = ...
                    sensor_pos_to_info_rbf_blackbox(sensor_posns_FP, threat_basis_data, ...
                    grid_world, sensor_noise);
                time_FP = toc;
                if PRINT_THINGS; fprintf('Finished in %.3f seconds\n\n',time_FP); end
                % -------------------------------------------------------------
                
                % greedy Mutual Information (MI)-------------------------------
                %             fprintf('Running greedy Mutual Information...\n')
                %             tic
                %             sensor_posns_MI = SP_greedyMI(Psi*Psi', ipas_total_measurements);
                %             [incurred_cost_MI, expected_cost_MI, path_cost_true_MI, path_cost_var_MI] = ...
                %         sensor_pos_to_info_rbf_blackbox(sensor_posns_MI, threat_basis_data, ...
                %             grid_world, sensor_noise);
                %             time_MI = toc;
                %             fprintf('Finished in %.3f seconds\n\n',time_MI)
                % -------------------------------------------------------------
                
                % greedy Mean Squared Error (MSE) -----------------------------
%                 if PRINT_THINGS; fprintf('Running greedy MSE...\n'); end
%                 tic
%                 sensor_posns_MSE = SP_greedyMSE(Psi, ipas_total_measurements);
%                 [incurred_cost_MSE, expected_cost_MSE, path_cost_true_MSE, path_cost_var_MSE] = ...
%                     sensor_pos_to_info_rbf_blackbox(sensor_posns_MSE, threat_basis_data, ...
%                     grid_world, sensor_noise);
%                 time_MSE = toc;
%                 if PRINT_THINGS; fprintf('Finished in %.3f seconds\n\n',time_MSE); end
                % -------------------------------------------------------------
                
                % RANDOM sensor placement
                if PRINT_THINGS; fprintf('Running Random without replacement...\n'); end
                tic
                sensor_posns_RANDOM = sensor_reconf_random_wo_replace(grid_world, ipas_total_measurements);
                [incurred_cost_RAND, expected_cost_RAND, path_cost_true_RAND, path_cost_var_RAND] = ...
                    sensor_pos_to_info_rbf_blackbox(sensor_posns_RANDOM, threat_basis_data, ...
                    grid_world, sensor_noise);
                time_RAND = toc;
                if PRINT_THINGS; fprintf('Finished in %.3f seconds\n\n',time_RAND); end
                % -------------------------------------------------------------
                
                % ====================== STORE ALL DATA =======================
                sim_data_record(grid_num, Np_num, Ns_num).threat_basis_data = threat_basis_data;
                
                sim_data_record(grid_num, Np_num, Ns_num).ipas.Np = Np^2;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.Ns = Ns;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.grid_size = curr_grid^2;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.noise_sigma = noise_sigma;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.iterations = ipas_iters;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.n_measurements = ipas_total_measurements;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.path_cost.true = ipas_path_cost_true;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.path_cost.expected = ipas_path_cost_expected;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.path_cost.incurred = ipas_path_cost_incurred;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.path_cost_var = ipas_path_cost_var;
                sim_data_record(grid_num, Np_num, Ns_num).ipas.comp_time = ipas_time;
                
                sim_data_record(grid_num, Np_num, Ns_num).FP.Np = Np^2;
                sim_data_record(grid_num, Np_num, Ns_num).FP.Ns = Ns;
                sim_data_record(grid_num, Np_num, Ns_num).FP.grid_size = curr_grid^2;
                sim_data_record(grid_num, Np_num, Ns_num).FP.noise_sigma = noise_sigma;
                sim_data_record(grid_num, Np_num, Ns_num).FP.iterations = ceil(ipas_total_measurements/Ns);
                sim_data_record(grid_num, Np_num, Ns_num).FP.n_measurements = ipas_total_measurements;
                sim_data_record(grid_num, Np_num, Ns_num).FP.path_cost.true = path_cost_true_FP;
                sim_data_record(grid_num, Np_num, Ns_num).FP.path_cost.expected = expected_cost_FP;
                sim_data_record(grid_num, Np_num, Ns_num).FP.path_cost.incurred = incurred_cost_FP;
                sim_data_record(grid_num, Np_num, Ns_num).FP.path_cost_var = path_cost_var_FP;
                sim_data_record(grid_num, Np_num, Ns_num).FP.comp_time = time_FP;
                
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.Np = Np^2;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.Ns = Ns;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.grid_size = curr_grid^2;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.noise_sigma = noise_sigma;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.iterations = ceil(ipas_total_measurements/Ns);
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.n_measurements = ipas_total_measurements;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.path_cost.true = path_cost_true_MI;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.path_cost.expected = expected_cost_MI;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.path_cost.incurred = incurred_cost_MI;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.path_cost_var = path_cost_var_MI;
                %             sim_data_record(grid_num, Np_num, Ns_num).MI.comp_time = time_MI;
                
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.Np = Np^2;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.Ns = Ns;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.grid_size = curr_grid^2;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.noise_sigma = noise_sigma;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.iterations = ceil(ipas_total_measurements/Ns);
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.n_measurements = ipas_total_measurements;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.path_cost.true = path_cost_true_MSE;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.path_cost.expected = expected_cost_MSE;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.path_cost.incurred = incurred_cost_MSE;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.path_cost_var = path_cost_var_MSE;
%                 sim_data_record(grid_num, Np_num, Ns_num).MSE.comp_time = time_MSE;
                
                sim_data_record(grid_num, Np_num, Ns_num).RAND.Np = Np^2;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.Ns = Ns;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.grid_size = curr_grid^2;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.noise_sigma = noise_sigma;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.iterations = ceil(ipas_total_measurements/Ns);
                sim_data_record(grid_num, Np_num, Ns_num).RAND.n_measurements = ipas_total_measurements;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.path_cost.true = path_cost_true_RAND;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.path_cost.expected = expected_cost_RAND;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.path_cost.incurred = incurred_cost_RAND;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.path_cost_var = path_cost_var_RAND;
                sim_data_record(grid_num, Np_num, Ns_num).RAND.comp_time = time_RAND;
                
                % Overwrite .mat file each n_sim iterations so that if it
                % crashes we still have some data.
%                 save(fullfile(folder_name,sim_data_file_name),'sim_data_record', 'sim_meta_data')
                Ns_num = Ns_num + 1;
            end % End of current sensor Ns iteration
            Np_num = Np_num + 1;
        end % End of current parameter Np iteration
        grid_num = grid_num + 1;
    end % End of current grid size curr_grid iteration
    full_sims(sim_i).data = sim_data_record;
    parfor_progress; % Count
%     fprintf('\nEnd of Full Simulation %i/%i \n', sim_i, N_sim);
end % End of Full simulation sim_i
parfor_progress(0); % Clean up

save(fullfile(folder_name,sim_data_file_name),'full_sims', 'sim_meta_data')

dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
fprintf('\nSimulation runs ended at: %s\n', dt)

