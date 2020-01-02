function t = IsLos_discrete(p1, p2, Map, stepsize, x0)
% To test if there is a building that blocks the LOS between two points p1
% and p2 in 3D space based on a evaluvation map, Map, which is a matrix
% that contains the evaluation data.
%
% INPUT 
%   p1, p2      3D positions of the two points
%   Map         Nrow * Ncol evaluation map, where rows represent the x
%               indices, and columns represent the y indices
%   stepsize    Step size (grid size) in meters of the evaluation Map
%   x0          The left-bottom coordinates (2D) of the map

p1 = p1(:).';
p2 = p2(:).';
x0 = x0(:).';

[Nrow, Ncol] = size(Map);

delta = max(abs(p2(1:2) - p1(1:2)));
n_grid_intervals = round(delta / stepsize); % Note that: #(grid_points) = #(grid_interverals) + 1

if n_grid_intervals == 0
    t = 1;  % perpendicular, LOS
    return
end

% Extract the vector of elevations from the elevatio map, Map
map_grid_span = (p2 - p1) / n_grid_intervals;  % could be negative
n_grid_points = n_grid_intervals + 1;
grid_loc = zeros(n_grid_points, 3);
grid_loc(:, 1) = p1(1) + (0:n_grid_intervals).' * map_grid_span(1);
grid_loc(:, 2) = p1(2) + (0:n_grid_intervals).' * map_grid_span(2);
grid_loc(:, 3) = p1(3) + (0:n_grid_intervals).' * map_grid_span(3); 
           % coordinates in meter, (n_grid_points * 3) matrix

map_id = round((grid_loc(:, 1:2) - repmat(x0, n_grid_points, 1)) / stepsize) + ones(n_grid_points, 2);
           % matrix row and column indices for each grid, (n_grid_points * 2) matrix

valid_grid_id = map_id(:, 1) >= 1 & map_id(:, 1) <= Nrow ...
                        & map_id(:, 2) >= 1 & map_id(:, 2) <= Ncol;
                    
map_elevation = MatVal(Map, map_id(valid_grid_id, :));
line_height = grid_loc(valid_grid_id, 3);

if isempty(find(map_elevation > line_height, 1))
    t = 1;
else
    t = 0;
end

end


function V = MatVal(A, index2d)
% Input a set of indices {(x, y)}, output a set of values {A(x,y)}.
    nrow = size(A, 1);
    index1d = (index2d(:, 2)  - 1) * nrow + index2d(:, 1);
    V = A(index1d);
end