function optimal_path = greedy_trace(v_start, v_goal, vertex_data, known_verts)

if nargin == 3
	v_current	= v_goal;
	optimal_path= v_goal;

	while (v_current ~= v_start)
		v_current	= vertex_data(v_current).b;
		optimal_path= cat(2, v_current, optimal_path);
	end
	return
end

if nargin > 3
	[is_goal_known, id_current] = ismember(v_goal, known_verts);
	if ~is_goal_known, error('knownNodes does not contain goal!'); end
end

optimal_path = v_goal;
while (id_current ~= 1)
	id_current	= vertex_data(id_current).b;
	optimal_path= cat(2, known_verts(id_current), optimal_path);
end
