function h = min_threat_manhat_heur(id_current, search_data, grid_world)

id_goal		= search_data.v_goal;
grid_space = grid_world.spacing;
% min_threat_value = min(search_data.adjacency_matrix); % min(A(find(A)))
A = search_data.adjacency_matrix;
min_threat_value = full(min(A(find(A))));
% disp(min_threat_value)

dx = abs(grid_world.coordinates(1, id_current) - grid_world.coordinates(1, id_goal));
dy = abs(grid_world.coordinates(2, id_current) - grid_world.coordinates(2, id_goal));
min_steps = round(dx/grid_space) + round(dy/grid_space);
h = min_steps * min_threat_value;

% fprintf('grid_space: %d\n', grid_space)
% fprintf('Current Position: (%d, %d)\n',  grid_world.coordinates(1, id_current), grid_world.coordinates(2, id_current))
% fprintf('Goal Position: (%d, %d)\n',  grid_world.coordinates(1, id_goal), grid_world.coordinates(2, id_goal))
% fprintf('min_threat_value: %d\n', min_threat_value)
% fprintf('min_steps: %d\n', min_steps)
% fprintf('h = %d\n', h)
end