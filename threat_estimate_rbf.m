function [tpe_mean, tpe_covar] = threat_estimate_rbf(sensor_locations, ...
	grid_world, threat_basis_data, threat_parameters_true, sensor_noise)

% "sensor_location" is 1 x n vector of sensor grid point locations
% true parameters are needed only for sensor simulation, so no cheating

sensor_posn		= grid_world.coordinates(:, sensor_locations);
measured_threat = calc_sensor_measurement(threat_basis_data, threat_parameters_true, sensor_posn, sensor_noise);
measured_threat_wo_offset = measured_threat - threat_basis_data.offset;

H_measurement	= calc_rbf_value(threat_basis_data, sensor_posn);
noise_inv = sensor_noise.Rinv(1,1);
noise_inv_extended = diag(noise_inv*ones(1,size(sensor_posn,2)));

tpe_mean = H_measurement \ measured_threat_wo_offset;						% Least squares estimate of parameters
tpe_covar= inv( H_measurement' * noise_inv_extended * H_measurement );