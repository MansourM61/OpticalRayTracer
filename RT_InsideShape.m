function flag = RT_InsideShape(Point, Shape, Tol)

NoP = size(Shape, 1);  % number of points

NoH = 0;  % number of hits

for Index = 1:(NoP - 1)  % go through the points
    x0 = Shape(Index, 1);  % start x point
    y0 = Shape(Index, 2);  % start y point
    
    x1 = Shape(Index + 1, 1);  % end x point
    y1 = Shape(Index + 1, 2);  % end y point
    
    line = [x0, y0, x1, y1];  % line geometry
    
    [~, ~, ~, ~, flag_int] = RT_Intersection([Point, 0], line, Tol);  % find the intersection point
    
    dirFlag = flag_int(1);  % direction of propagation
    hitFlag = flag_int(2);  % hit/miss
    
    if(dirFlag && hitFlag)  % check if insection is valid
        NoH  = NoH + 1;  % increase the number of impact
    end
    
end

if(mod(NoH, 2) == 0)  % if the source is outside the shape
    flag = false;  % update the return flag
else   % if the source is inside the shape
    flag = true;  % update the return flag
end
