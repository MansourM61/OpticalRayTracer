function ShapePoints = RT_GeometryQuantizer(Geometries, Delta, Coeff)

MAX_MEM = 1000;  % maximum length of temporary memory

dx = Delta(1);  % x resolution
dy = Delta(2);  % y resolution

NoS = size(Geometries, 1);  % number of embedded shapes

Temp = zeros(MAX_MEM, 2); % temperary memory

Mem_Index = 1;  % memory index

for Index = 1:NoS
    t_s = Geometries{Index, 2}(1);  % start point of the shape
    t_e = Geometries{Index, 2}(2);  % end point of the shape
    
    if(t_s < t_e)  % if t parameter is increasing
        dt_0 = +min(dx, dy);  % initial delta t
    else % if t parameter is decreasing
        dt_0 = -min(dx, dy);  % initial delta t
    end
    
    dt = dt_0;  % intialize delta t
    
    t = t_s;  % start from the beginning
    p_s = Geometries{Index, 1}(t);  % start point
    Temp(Mem_Index, :) = p_s;  % insert the first point
    
    t_1 = t;  % previous t parameter
    
    if(t_s < t_e)  % if t parameter is increasing
        Cond = t < t_e;  % loop condition
    else % if t parameter is decreasing
        Cond = t > t_e;  % loop condition
    end
    
    while( Cond )  % while the piece is not finished
        t = t_1 + dt;  % calculate next t
        p_s = Geometries{Index, 1}(t);  % next point
        
        while( (abs(p_s(1) - Temp(Mem_Index, 1)) > dx) || (abs(p_s(2) - Temp(Mem_Index, 2)) > dy) )
            dt = dt*Coeff;  % set a new delta t
            t = t_1 + dt;  % calculate next t
            p_s = Geometries{Index, 1}(t);  % next point
        end
        
        Mem_Index = Mem_Index + 1;  % update memory index
        t_1 = t;  % update previous t parameter
        dt = dt_0;  % intialize delta t
        Temp(Mem_Index, :) = p_s;  % insert the point
        
        if(t_s < t_e)  % if t parameter is increasing
            Cond = t < t_e;  % loop condition
        else % if t parameter is decreasing
            Cond = t > t_e;  % loop condition
        end
        
    end
end

ShapePoints = Temp([1:Mem_Index - 1, 1], :);  % generate quantised shape matrix