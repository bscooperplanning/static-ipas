function threat_value_measured = calc_sensor_measurement(threat_basis_data, threat_parameters_true, sensor_posn, sensor_noise)

noise = sensor_noise.R(1,1);
noise_extended = diag(noise*ones(1,size(sensor_posn,2)));

threat_value_measured = (calc_threat_rbf(threat_basis_data, threat_parameters_true, sensor_posn))' ...
	+ sqrt(noise_extended)*randn(size(sensor_posn, 2), 1);

