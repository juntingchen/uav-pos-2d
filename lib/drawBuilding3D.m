function drawBuilding3D(Blds, U)
% drawBuilding3D(Blds, U)
% 
% INPUT
%   Blds    Nblds * 2 Cell array, each row: a matrix to specify the
%           vertices of a polygon, and a scalar for the building height
%   U       Urban environment configuration

if nargin > 1
    LengthX = U.LengthX;
    LengthY = U.LengthY;
end
    
Nbld = size(Blds, 1);
wallcolor = [0.8 0.8 0.8];
for ib = 1:Nbld
%     Building = [X1 Y1
%                 X2 Y1
%                 X2 Y2
%                 X1 Y2];
    x1 = Blds{ib, 1}(1, 1);
    x2 = Blds{ib, 1}(2, 1);
    y1 = Blds{ib, 1}(1, 2);
    y2 = Blds{ib, 1}(3, 2);
    h = Blds{ib, 2};

    p1 = [x1, y1, 0];
    p2 = [x2, y1, 0];
    p3 = [x2, y2, 0];
    p4 = [x1, y2, 0];
    
    p5 = [x1, y1, h];
    p6 = [x2, y1, h];
    p7 = [x2, y2, h];
    p8 = [x1, y2, h];
    
    alpha_val = 0.5;
    
    % bottom
    h = fill3([p1(1), p2(1), p3(1), p4(1)], ...
          [p1(2), p2(2), p3(2), p4(2)], ...
          [p1(3), p2(3), p3(3), p4(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    hold on 
    alpha(h, alpha_val);
    
    % back
    h = fill3([p3(1), p4(1), p8(1), p7(1)], ...
          [p3(2), p4(2), p8(2), p7(2)], ...
          [p3(3), p4(3), p8(3), p7(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    alpha(h, alpha_val);
    hold on  
    
    % left
    h = fill3([p4(1), p1(1), p5(1), p8(1)], ...
          [p4(2), p1(2), p5(2), p8(2)], ...
          [p4(3), p1(3), p5(3), p8(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    hold on  
    alpha(h, alpha_val);
    
    % top
    h = fill3([p5(1), p6(1), p7(1), p8(1)], ...
          [p5(2), p6(2), p7(2), p8(2)], ...
          [p5(3), p6(3), p7(3), p8(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    hold on  
    alpha(h, alpha_val);
    
    % right
    h = fill3([p2(1), p3(1), p7(1), p6(1)], ...
          [p2(2), p3(2), p7(2), p6(2)], ...
          [p2(3), p3(3), p7(3), p6(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    hold on 
    alpha(h, alpha_val);
    
    % front
    h = fill3([p1(1), p2(1), p6(1), p5(1)], ...
              [p1(2), p2(2), p6(2), p5(2)], ...
              [p1(3), p2(3), p6(3), p5(3)], wallcolor);
    % set(h, 'edgecolor', 'none');
    hold on 
    alpha(h, alpha_val);
    
end
hold off

if nargin > 1
    xlim([0, LengthX]);
    ylim([0, LengthY]);

    xlabel('x');
    ylabel('y');
end

