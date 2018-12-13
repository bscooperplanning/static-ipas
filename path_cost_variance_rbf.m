function path_var = path_cost_variance_rbf(opt_path, threat_basis_data_orig,...
                    threat_basis_data_subset, grid_world, tpe_covar)

% If tpe_covar is from a subset of basis smaller than Np, then augment the
% covar with extra row&col with large variance in diagonal
if size(tpe_covar,2) < threat_basis_data_orig.n_threat_parameters
    uibs = 1:threat_basis_data_orig.n_threat_parameters; % Un-Identified Basis Set (uibs)
    uibs(threat_basis_data_subset.basis_subset) = [];
    for bi = 1:length(uibs)
        uib = uibs(bi);
        bcol = zeros(size(tpe_covar,2),1);
        brow = zeros(1, size(tpe_covar,1) + 1);
        brow(uib) = 10e2; % some arbitrary large variance value
        tpe_covar = [tpe_covar(:, 1:uib-1) bcol tpe_covar(:, uib:end)];
        tpe_covar = [tpe_covar(1:uib-1, :); brow; tpe_covar(uib:end, :)];
    end
end

path_var = 0;
for p = 1:length(opt_path)
    H_m = calc_rbf_value(threat_basis_data_orig, grid_world.coordinates(:, opt_path(p)));
    path_var = path_var + H_m*tpe_covar*H_m';
%     path_var = path_var + H_m*noise_extended*H_m';
%     fprintf('path step: %i, path_var = %f\n',p, path_var)
end

end