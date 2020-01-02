function t = IsLosK_discrete(p1, p2, Maps, stepsize, x0)
% Version 2: Extension to K segment case
% 
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

K = length(Maps) + 1;
[Nrow, Ncol] = size(Maps{1});

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
                    
line_height = grid_loc(valid_grid_id, 3);

t = 1;
for k = 1:K - 1
    map_elevation = MatVal(Maps{k}, map_id(valid_grid_id, :));
    if isempty(find(map_elevation > line_height, 1))
        blockage_type = 0;
    else
        blockage_type = k / (K - 1);
    end
    t = min(t, 1 - blockage_type);
end

% % Debug area
% if norm(p2(1:2) - [135.1, 782.7]) < 0.5 || norm(p2(1:2) - [135.1, 782]) < 0.5
%     xmin = min(map_id(:, 1)); xmax = max(map_id(:, 1));
%     ymin = min(map_id(:, 2)); ymax = max(map_id(:, 2));
%     figure, 
%     [xmat, ymat] = meshgrid(xmin:xmax, ymin:ymax);
%     surf(xmat.', ymat.', Maps{1}(xmin:xmax, ymin:ymax) + Maps{2}(xmin:xmax, ymin:ymax));
%     view(0, 90);hold on
%     for i = 1:n_grid_points
%         plot3(map_id(i, 1), map_id(i, 2), 50, 'ro');
%     end
%     hold off
% end
% % Debug area end

end 

function V = MatVal(A, index2d)
% Input a set of indices {(x, y)}, output a set of values {A(x,y)}.
    nrow = size(A, 1);
    index1d = (index2d(:, 2)  - 1) * nrow + index2d(:, 1);
    V = A(index1d);
end