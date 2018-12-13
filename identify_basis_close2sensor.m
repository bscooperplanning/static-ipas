function [ identified_bases ] = identify_basis_close2sensor(sensor_location, ...
	grid_world, threat_basis_data)

identified_bases = [];
if ( isempty(sensor_location) )
    return;
end

% Find identifiable basis functions
basis_oi = false(1, threat_basis_data.n_threat_parameters);
% Identifiability defined as within area of significant support: 3*sigma
max_separation = 3*sqrt(threat_basis_data.basis_parameters.var);
for m1 = 1:threat_basis_data.n_threat_parameters
    basis_sensor_vec = kron( threat_basis_data.basis_parameters.mean(:, m1), ...
		ones(1, numel(sensor_location)) ) - grid_world.coordinates(:, sensor_location);
    % Euclidean distance of basis from sensor locations
    basis_sensor_separation = sqrt(basis_sensor_vec(1,:).^2 + basis_sensor_vec(2,:).^2);
    
    % If the closest sensor to basis is less than 3*sigma, it's identifiable
    if ( min(basis_sensor_separation) <= max_separation )
        basis_oi(m1) = true;
    end
end
tmp00 = 1:threat_basis_data.n_threat_parameters;
identifiable_bases = tmp00(basis_oi);
n_basis_oi = sum(basis_oi);

%------ Rank basis by closeness to sensors
% for m1 = 1:numel(identifiable_bases)
%     basis_sensor_vec = kron( threat_basis_data.basis_parameters.mean(:, m1), ...
% 		ones(1, numel(sensor_location)) ) - grid_world.coordinates(:, sensor_location);
%     basis_sensor_sep(m1) = min(sqrt(basis_sensor_vec(1,:).^2 + basis_sensor_vec(2,:).^2));
% end

%----- Find closest basis to each sensor
for m1 = 1:numel(sensor_location)
    sensor_basis_vec = kron(grid_world.coordinates(:, sensor_location(m1)) , ...
		ones(1, n_basis_oi) ) - threat_basis_data.basis_parameters.mean(:, identifiable_bases);
    [sensor_basis_sep(m1), nearest_bases(m1)] = min(sqrt(sensor_basis_vec(1,:).^2 + sensor_basis_vec(2,:).^2));
end
identified_bases = unique(identifiable_bases(nearest_bases),'stable');
end