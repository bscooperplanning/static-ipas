function sensor_locations = sensor_reconf_random_wo_replace(grid_world, n_sensor)
% random sensor locations without replacement (can't put 2 sensors in the
% same grid point.
sensor_locations = zeros(1, n_sensor);
%----- Arbitary placement
for m1 = 1:n_sensor
    if m1 > grid_world.n_grid_points
        sensor_locations = sensor_locations(1:grid_world.n_grid_points);
        break
    end
	this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
	while any(this_location == sensor_locations(1:m1))
		this_location = 1 + round((grid_world.n_grid_points - 1)*rand);
	end
	sensor_locations(m1) = this_location;
end

end