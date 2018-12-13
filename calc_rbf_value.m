function H_measurement = calc_rbf_value(threat_basis_data, posn)

H_measurement = zeros(size(posn, 2), threat_basis_data.n_threat_parameters);

for m1 = 1:size(posn, 2)
	posn_vec = [posn(1, m1)*ones(1, threat_basis_data.n_threat_parameters); ...
		posn(2, m1)*ones(1, threat_basis_data.n_threat_parameters)];

	H_measurement(m1, :) = exp((-1 / (2 * threat_basis_data.basis_parameters.var)) .* ...
		((posn_vec(1, :) - threat_basis_data.basis_parameters.mean(1, :)).^2 + ...
		 (posn_vec(2, :) - threat_basis_data.basis_parameters.mean(2, :)).^2) );
end