% The sensor reconfiguration method only returns enough sensors to measure
% the unknown relevant basis at each iteration. This reduces the overall
% movement of sensors which may be unnecessary.

function [sensor_location, basis_intersect] = ...
    sensor_reconf_min_basis(threat_basis_data_orig, threat_basis_data_subset,...
    grid_world, n_sensor, opt_path)

% Places sensors at intersection of basis in range of path and basis that
% are unmeasured/unidentified. With a consequence that in later iterations
% when path has been mostly measured, the sensors will cluster on unknown
% basis

%----- Spatial basis functions of interest
basis_oi = false(1, threat_basis_data_orig.n_threat_parameters);
% max_seperation should have distance units, param .var is [L]^2 since
% variance is sigma^2, and std dev. is sigma. Changinge to sqrt(.var)
% max_separation = 2*threat_basis_data.basis_parameters.var; % If any coordinate of path is within 1-sigma of basis center
max_separation = 3*sqrt(threat_basis_data_orig.basis_parameters.var);
for m1 = 1:threat_basis_data_orig.n_threat_parameters
	basis_path_vec			= kron( threat_basis_data_orig.basis_parameters.mean(:, m1), ...
		ones(1, numel(opt_path)) ) - grid_world.coordinates(:, opt_path);
	basis_path_separation	= basis_path_vec(1, :).^2 + basis_path_vec(2, :).^2;
	
	if ( min(basis_path_separation) <= max_separation^2 )
		basis_oi(m1) = true;
	end
end

% basis_oi
tmp00 = 1:threat_basis_data_orig.n_threat_parameters;
reduced_parameters = tmp00(basis_oi);

% Identify basis with high variance
uibs = 1:threat_basis_data_orig.n_threat_parameters; % Un-Identified Basis Set (uibs)
uibs(threat_basis_data_subset.basis_subset) = [];
basis_intersect = intersect(reduced_parameters, uibs);
n_basis_intersect = length(basis_intersect);

% If we identify more unknown basis along path than sensors, limit to
% n_sensors.
if n_basis_intersect > n_sensor
   n_basis_intersect = n_sensor; 
end

% If all relevant basis have been measured, intersection is []
% Default to full set of reduced parameters. IPAS should terminate, so we
% shouldn't need to move any sensors. What is best strategy in this case?
if length(basis_intersect) < 1
    basis_intersect = reduced_parameters;
    n_basis_intersect = n_sensor;
end

%----- Rank grid points by closeness to spatial bases of interest
n_basis_i		= length(basis_intersect);
basis_grid_sep	= zeros(grid_world.n_grid_points, 1);

for m1 = 1:grid_world.n_grid_points
	basis_grid_vec		= kron( grid_world.coordinates(:, m1), ones(1, n_basis_i) ) - ...
		threat_basis_data_orig.basis_parameters.mean(:, basis_intersect);
	basis_grid_sep(m1)	= min( basis_grid_vec(1, :).^2 + basis_grid_vec(2, :).^2 );
	
end

ranked_basis_sep= sortrows([(1:grid_world.n_grid_points)', basis_grid_sep], 2);
sensor_location = ranked_basis_sep(1:n_basis_intersect, 1)';