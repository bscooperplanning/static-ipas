%{
Copyright (c) 2015 Raghvendra V. Cowlagi. All rights reserved.

Copyright notice: 
=================
No part of this work may be reproduced without the written permission of
the copyright holder, except for non-profit and educational purposes under
the provisions of Title 17, USC Section 107 of the United States Copyright
Act of 1976. Reproduction of this work for commercial use is a violation of
copyright.


Disclaimer:
===========
This software program is intended for educational and research purposes.
The author and the institution with which the author is affiliated are not
liable for damages resulting the application of this program, or any
section thereof, which may be attributed to any errors that may exist in
this program.


Author information:
===================
Raghvendra V. Cowlagi, Ph.D,
Assistant Professor, Aerospace Engineering Program,
Department of Mechanical Engineering, Worcester Polytechnic Institute.
 
Higgins Laboratories, 247,
100 Institute Road, Worcester, MA 01609.
Phone: +1-508-831-6405
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi


The author welcomes questions, comments, suggestions for improvements, and
reports of errors in this program.

Program Description
===================
The mother of all A* versions.
	* Either struct or matrix will do
	* Even a function specifying neighbors will do, for search as you go
	* Multiple goals when any one will do
	* Multiple goals when cost-to-come for all is required
	* Single start vertex, for multiple, adjust calling fcn
%}

function [vertex_data, known_verts] = astar_hfunc(search_data)
% nodeS		: Start node
% nodeG		: Goal node set
% heur		: Heuristic (N x 1 vector, where N is number of nodes),
%			can also be function handle
% A			: Transition cost as struct
% G			: Transition cost as matrix
% nhbrHandle: Handle to function which returns (forward) adjacent states,
%			and corresponding transition costs
% nhbrhdData: Any additional data that 'nhbrHandle' and/or 'heur may
%			require, a struct that 'nhbrHandle' or 'heur' decode internally
% allGoal	: Whether cost-to-come to any goal or all goals

% nodeData(n).mk= marker, 0 = NEW, 1 = OPEN, 2 = CLOSED
% nodeData(n).d	= cost to come
% nodeData(n).b	= backpointer


id_start	= search_data.v_start;
id_goal		= search_data.v_goal;
heuristic	= search_data.heuristic;
A			= search_data.adjacency_struct;
G			= search_data.adjacency_matrix;
fcn_find_nhbr = search_data.fcn_find_nhbr;
all_goal	= search_data.all_goal;
grid_world  = search_data.grid_world;


is_large = 0;																% If graph is large
if numel(A)
	n_vertices	= numel(A);
elseif numel(G)
	n_vertices	= size(G, 1);
else
	n_vertices		= 1e3;													% Graph too large, initialize nodeStruct with arbitrary size
	is_large		= 1;
	n_known_verts	= 0;
	known_verts		= zeros(n_vertices, 1);
end

vert_struct	= struct('id', 0, 'mk', 0, 'd', Inf, 'b', []);
vertex_data	= repmat(vert_struct, 1, n_vertices);
if is_large
	vertex_data(1).id	= id_start;
	vertex_data(1).mk	= 1;	vertex_data(1).d	= 0;
	
	n_known_verts				= n_known_verts + 1;
	known_verts(n_known_verts)	= id_start;
else
	vertex_data(id_start).mk	= 1;	vertex_data(id_start).d	= 0;
end
	

n_open	= 1;
if is_large
	open_list = [1 heuristic(id_start, search_data)];
else
    if isa(heuristic, 'function_handle')
        open_list = [id_start heuristic(id_start, search_data, grid_world)];
    else
        open_list = [id_start heuristic(id_start)];
    end
end
goal_closed	= 0;

n_iter	= 0;
while (n_open ~= 0) && (~goal_closed)
% 	clc;
% 	fprintf('Number of iterations  : %i\n', nIter);
% 	fprintf('Number of OPEN nodes  : %i\n', nOpen);
% 	fprintf('Number of known nodes : %i\n\n', nKnownNodes);
	
	n_iter		= n_iter + 1;
	v_current	= open_list(1, 1);								% Get node from top of (sorted) open stack
	vertex_data(v_current).mk = 2;											% Mark that node as dead
	
	n_open		= n_open - 1;
	open_list(1, :) = [];	
	
	if numel(A)																% Transition costs in struct
		nhbrs	= A(v_current).nhbrs;
		costs	= A(v_current).costs;
	elseif numel(G)															% Transition costs in matrix
		nhbrs	= find(G(v_current,:));
		costs	= full(G(v_current, nhbrs));
	else																	% Other function handle for nhbrs and costs
		[nhbrs, costs] = fcn_find_nhbr(vertex_data(v_current).id, search_data);
	end
	
	for k = 1:numel(nhbrs)													% For all neighbors
		if is_large
			[known_nhbr, idx_nhbr]	= ismember(nhbrs(k), known_verts);
		else
			known_nhbr				= 0;
		end
		
		if is_large && known_nhbr
			v_new					= idx_nhbr;
		elseif is_large
			v_new					= n_known_verts + 1;
			vertex_data(v_new).id	= nhbrs(k);
			if v_new > n_vertices				
				vertex_data(v_new).mk	= 0;
				vertex_data(v_new).b	= [];
				vertex_data(v_new).d	= Inf;
			end
			
			n_known_verts				= v_new;
			known_verts(n_known_verts)	= nhbrs(k);
		else
			v_new	= nhbrs(k);			
		end
		cost_new	= costs(k);												% Cost to go from act to new
		
% 		[v_current v_new cost_new]
		
		if vertex_data(v_new).mk == 0										% Unvisited
			vertex_data(v_new).mk	= 1;									% Mark open
			vertex_data(v_new).d	= vertex_data(v_current).d + cost_new;	% Update c2come of newly visited state
			vertex_data(v_new).b	= v_current;
			
			if is_large
				tmpOpen = bin_sort(open_list(1:n_open, :), [v_new ...
					 vertex_data(v_new).d + heuristic(nhbrs(k), search_data)], 2);
            else
                if isa(heuristic, 'function_handle')
                    heur_value = heuristic(v_new, search_data, grid_world);
                else
                    heur_value = heuristic(v_new);
                end
				tmpOpen = bin_sort(open_list(1:n_open, :), ...
					[v_new vertex_data(v_new).d + heur_value], 2);
			end
			if numel(tmpOpen) == 0
				n_open	= 0;
				open_list	= [];
			else
				n_open	= size(tmpOpen, 1);
				open_list(1:n_open, :)	= tmpOpen;								% Add [nodeNew cost] to sorted open list
			end			
		elseif vertex_data(v_new).mk == 1									% Already open, update c2come if necessary
			if vertex_data(v_new).d > vertex_data(v_current).d + cost_new
				vertex_data(v_new).d	= vertex_data(v_current).d + cost_new;
				vertex_data(v_new).b	= v_current;
				
				[~, loc] = ismember(v_new, open_list(1:n_open, 1));
				open_list(loc, :)= [];		n_open = n_open - 1;
				
				if is_large
					tmpOpen = bin_sort(open_list(1:n_open, :), [v_new ...
						vertex_data(v_new).d + heuristic(nhbrs(k), search_data)], 2);					
                else
                    if isa(heuristic, 'function_handle')
                        heur_value = heuristic(v_new, search_data, grid_world);
                    else
                        heur_value = heuristic(v_new);
                    end
					tmpOpen = bin_sort(open_list(1:n_open, :), ...
						[v_new vertex_data(v_new).d + heur_value], 2);
				end
				if numel(tmpOpen) == 0
					n_open	= 0;
					open_list	= [];
				else
					n_open	= size(tmpOpen, 1);
					open_list(1:n_open, :)	= tmpOpen;							% Add [nodeNew cost] to sorted open list
				end
			end
		end
	end
	
	if is_large
		if ~all_goal
			for k = 1:numel(id_goal)
				[is_goal_known, idx_goal] = ismember(id_goal(k), known_verts);
				if is_goal_known && vertex_data(idx_goal).mk == 2, goal_closed = 1; break; end
			end
		else			
			if ~all(ismember(id_goal, known_verts)), continue; end
			goal_closed = 1;
			for k = 1:numel(id_goal)
				[~, idx_goal] = ismember(id_goal(k), known_verts);
				if vertex_data(idx_goal).mk == 2, goal_closed = 0; break; end
			end
		end
	else
		if ~all_goal																% Any goal will do
			for k = 1:numel(id_goal)
				if vertex_data(id_goal(k)).mk == 2, goal_closed = 1; break; end
			end
		else																	% All goals must pop
			goal_closed = 1;
			for k = 1:numel(id_goal)
				if vertex_data(id_goal(k)).mk ~= 2, goal_closed = 0; break; end
			end
		end
	end
end
if is_large
	known_verts = known_verts(1:n_known_verts);
else
	known_verts = [];
end


%**************************************************************************
function A = bin_sort(A, B, c)
%--------------------------------------------------------------------------
% The rows of B are inserted into A, while sorting (ascending) according to
% column c. Both A and B have the same number of columns. A is assumed to
% sorted ascending.

[rA, cA] = size(A);
[rB, cB] = size(B);

if numel(A) == 0, A = B; return; end
if numel(B) == 0, return; end
if cB ~= cA, error('A and B must have same number of columns!\n'); end

for count = 1:rB
	thisIns		= B(count, :);
	thisCost	= thisIns(1, c);
	
	if ismember(thisIns, A, 'rows')
		fprintf('This one came back!\t\t'); disp(thisIns)
		redn = redn + 1;
		continue;
	end

	if A(rA, c) <= thisCost													% Insert after last row
		A	= cat(1, A, thisIns);
		rA	= rA + 1;
		continue;
	elseif A(1, c) >= thisCost												% Insert before first row
		A	= cat(1, thisIns, A);
		rA	= rA + 1;
		continue;
	end
	
	nCand	= rA;															% Number of candidate rows in A that can have greater cost
	testRow	= 0;
	dirn	= 1;	
	while nCand > 1
		p		= floor(nCand/2);
		testRow = testRow + dirn*p;
		dirn	= 2*(A(testRow, c) < thisCost) - 1;
		nCand	= nCand - p;
	end	

	insRow = testRow + (dirn + 1)/2;										% Insert before this row in A
	A	= cat(1, A(1:(insRow - 1), :), thisIns, A(insRow:rA, :));
	rA	= rA + 1;
end
%**************************************************************************