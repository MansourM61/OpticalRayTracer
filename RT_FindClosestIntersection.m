function [Index_geo, Index_int] = RT_FindClosestIntersection(Ray, Shapes, TochBound, Tol)

Index_int = [];  % intersection index
Index_geo = [];  % intersection index
d_min = [];  % minimum distance

NoG = length(Shapes);  % number of embedded shapes

for Index_G = 1:NoG  % go though each geometry
    Points = Shapes{Index_G};  % pick the geometry
    
    NoP = size(Points, 1);  % number of points

    for Index = 1:(NoP - 1)  % go through all the points and find the closest intersection
        
        if (Index_G == TochBound(1)) && (Index == TochBound(2))
            continue;
        end

        x0 = Points(Index, 1);  % start x point
        y0 = Points(Index, 2);  % start y point
        
        x1 = Points(Index + 1, 1);  % end x point
        y1 = Points(Index + 1, 2);  % end y point
        
        line = [x0, y0, x1, y1];  % line geometry
        
        [~, ~, ~, dist, flag_int] = RT_Intersection(Ray, line, Tol);  % find the intersection point
        
        dirFlag = flag_int(1);  % direction of propagation
        hitFlag = flag_int(2);  % hit/miss
        
        if(dirFlag && hitFlag)  % check if insection is valid
            if(isempty(d_min))  % if this is the first intersection?
                d_min = dist;  % update the minimum distance
                Index_geo = Index_G;  % update the index of geometry
                Index_int = Index;  % update the index of minimum intersection
            else  % if this is not the first intersection?
                if(dist < d_min)  % if the new intersection is closer?
                    d_min = dist;  % update the minimum distance
                    Index_geo = Index_G;  % update the index of geometry
                    Index_int = Index;  % update the index of minimum intersection
                end
            end
        end
    end
end