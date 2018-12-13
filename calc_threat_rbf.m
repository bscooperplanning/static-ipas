function threat_value = calc_threat_rbf(threat_basis_data, threat_parameters, posn)

% "posn" is 2 x n vector, where each column has 2D position coordinates

H_measurement	= calc_rbf_value(threat_basis_data, posn);
threat_value	= threat_basis_data.offset + (H_measurement * threat_parameters)';
