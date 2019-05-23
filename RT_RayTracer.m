function [RT_Array] = RT_RayTracer(Source, Geometries, n_background, n_shape, RT_Param, Tol)

DEBUG = false;  % no debugging message

NoL = 100;  % binary-tree ray array size

x_s = Source(1);  % source x point
y_s = Source(2);  % source y point
t_s = Source(3);  % source angle (Deg)
p_s = Source(4);  % source power

norm_len = RT_Param(1);  % normal line length
max_rt_d = RT_Param(2);  % maximum ray tracing depth
valid_ratio = RT_Param(3);  % minimum valid power ratio

BT(NoL) = struct('pos', [0, 0], 'ang', 0, 'power', 0,...
                 't_index', 0, 'r_index', 0, 'p_index', 0,...
                 'b_index', [0, 0], 'normal', [0, 0, 0, 0]);  % binary-tree ray array

for Index = 1:NoL  % go through all the elements in the tree
    BT(Index).pos = [0, 0];  % set all the nodes position to [0, 0]
    BT(Index).ang = 0;  % set all the nodes angle to 0 Deg
    BT(Index).power = 0;  % set all the nodes power to 0
    BT(Index).t_index = -2;  % set all the nodes to not been decided
    BT(Index).r_index = -2;  % set all the nodes to not been decided
    BT(Index).p_index = -2;  % set all the nodes to not been decided
    BT(Index).b_index = [0, 0];  % set all the nodes to not on boundary
    BT(Index).normal = [0, 0, 0, 0];  % set all the nodes to not on boundary
end

TR_depth = 1;  % initial ray tracing depth

Index_c = 1;  % index of current element
Index_e = 2;  % index of next available element
Index_b = [0, 0];  % index of boundary; no boundary is touched

BT(Index_c).pos = [x_s, y_s];  % update current node position
BT(Index_c).ang = t_s;  % update current node angle
BT(Index_c).power = p_s;  % update current node power
BT(Index_c).p_index = 0;  % update current node power

% p_index (= 0 source, <> 0 others)
% r/l_index (= 0 termination, = -1 farfield, = -2 not been decided)

while(true)  % loop throught the traces
    
    Index_t = BT(Index_c).t_index;  % get the transmit index
    Index_r = BT(Index_c).r_index;  % get the receive index
    Index_p = max(0, BT(Index_c).p_index);  % get the parent index
    
    if( (Index_t ~= -2) && (Index_r ~= -2) && (Index_p == 0) )  % if current node is the source and both sides are decided or reflection is decided and transmit is tir
        if(DEBUG == true)
            disp('end of tracing');  % end the loop
        end
        break;  % leave the loop
    end
    
    if (Index_p == 0)  % if current node is root
        if ( BT(Index_c).t_index == -2 ) % if transmit side is not decided
            calcFlag = 'T';  % set the flag to calculate the transmit
        elseif ( BT(Index_c).r_index == -2 )  % if receive side is not decided
            calcFlag = 'R';  % set the flag to calculate the receive
        else  % default state at root
            Index_c = 1;  % set the current element to the root
            Index_b = [0, 0];  % set boundary collision for root
            TR_depth = 1;  % reset the ray tracing depth
        end
        
    elseif (Index_p ~= 0) && (BT(Index_p).t_index ~= 0)  % if current node is not root and there is no tir
        if ( BT(Index_c).t_index == -2 ) % if transmit side is not decided
            calcFlag = 'T';  % set the flag to calculate the transmit
        elseif ( BT(Index_c).r_index == -2 )  % if receive side is not decided
            calcFlag = 'R';  % set the flag to calculate the receive
        else  % if both sides are decided
            Index_c = BT(Index_c).p_index;  % set the current element to the parent
            Index_b = BT(Index_p).b_index;  % update the touched boundary
            TR_depth = TR_depth - 1;  % update the depth of ray tracing
            continue;
        end
        
    elseif (Index_p ~= 0) && (BT(Index_p).t_index == 0)  % if current node is not root and there is tir
        if ( BT(Index_c).r_index == -2 )  % if receive side is not decided
            calcFlag = 'R';  % set the flag to calculate the receive
        else  % if both sides are decided
            Index_c = BT(Index_c).p_index;  % set the current element to the parent
            Index_b = BT(Index_p).b_index;  % update the touched boundary
            TR_depth = TR_depth - 1;  % update the depth of ray tracing
            continue;
        end
        
    else  % to avoid any unknown state
        if(DEBUG == true)  % if debugging is enabled
            print('Unknown state!');  % print the error
        end
    end
    
    ray = [BT(Index_c).pos, BT(Index_c).ang];  % ray geometry
    
    [Index_Geo, Index_int] = RT_FindClosestIntersection(ray, Geometries, Index_b, Tol);  % find closest intersection
    
    if ( isempty(Index_int) )  % if there is no intersection
        BT(Index_c).t_index = -1;  % set all the nodes to farfield
        BT(Index_c).r_index = -1;  % set all the nodes to farfield
        
        if (Index_p ~= 0)  % if the current node is not the root
            Index_c = BT(Index_c).p_index;  % update the parent index to the node above
            TR_depth = TR_depth - 1;  % update the depth of ray tracing
        end
        
        %        Index_p = max(1, BT(Index_c).p_index);  % get the parent index
        Index_b = BT(Index_c).b_index;  % update the touched boundary
        
        continue;  % start over the loop
    end
    
    Index_b = [Index_Geo, Index_int];  % update the touched boundary
    
    x0 = Geometries{Index_Geo}(Index_int, 1);  % intersected boundary x point
    y0 = Geometries{Index_Geo}(Index_int, 2);  % intersected boundary y point
    x1 = Geometries{Index_Geo}(Index_int + 1, 1);  % intersected boundary x point
    y1 = Geometries{Index_Geo}(Index_int + 1, 2);  % intersected boundary y point
    
    src = ray;  % source geometry
    line = [x0, y0, x1, y1];  % boundary geometry
    
    [v_s, v_b, p_i, ~, flags] = RT_Intersection(src, line, Tol);  % find the intersection
    
    insideFlag = flags(3);  % inside/outside flag
    
    if insideFlag  % if point is inside the shape
        n_ray = n_shape(Index_Geo);  % set the ray refractive index to geometry
        n_bnd = n_background;  % set the boundary refractive index to geometry
    else  % if point is outside the shape
        n_ray = n_background;  % set the ray refractive index to geometry
        n_bnd = n_shape(Index_Geo);  % set the boundary refractive index to geometry
    end
    
    source = [src, n_ray];  % ray geometry
    boundary = [line, n_bnd];  % boundary geometry
    
    [v_t, v_r, Rs, Rp, pn, tir] = RT_SnellsLaw(v_s, source, v_b, boundary, p_i, norm_len, Tol);  % calculate transmittance/reflectance
    
    Ts = 1 - Rs;  % transmittance = tranmsissivity = power transmission coefficient, s polarization
    Tp = 1 - Rp;  % transmittance = tranmsissivity = power transmission coefficient, p polarization
    
    R_eff = (Rs + Rp)/2;  % effective reflectance
    T_eff = (Ts + Tp)/2;  % effective transmittance
    
    if ( calcFlag == 'T' )  % if transmit is to be calculated
        
        if tir == true  % if there is tir
            BT(Index_c).t_index = 0;  % set the transmit node to termination
            Index_p = max(1, BT(Index_c).p_index);  % get the parent index
            Index_b = BT(Index_c).b_index;  % update the touched boundary
            
        else
            BT(Index_c).t_index = Index_e;  % set the transmit index to the next available element
            
            Index_p = Index_c;  % set parent index to current element
            Index_c = Index_e;  % set current index to new element
            
            Ang = atan2d(v_t(2), v_t(1));  % transmit ray angle
            p_t = BT(Index_p).power*T_eff;  % transmitted power
            
            BT(Index_c).pos = p_i;  % update current node position with intersection
            BT(Index_c).ang = Ang;  % update current node angle with calculated angle
            BT(Index_c).power = p_t;  % update current node power
            BT(Index_c).p_index = Index_p;  % update current node power
            BT(Index_c).b_index = Index_b;  % update the boundary collision
            BT(Index_c).normal = pn;  % update the normal line
            
            Index_e = Index_e + 1;  % update the next available element index
            
            TR_depth = TR_depth + 1;  % update the depth of ray tracing
            
        end
        
    elseif ( calcFlag == 'R' )  % if reflect is to be calculated
        
        BT(Index_c).r_index = Index_e;  % set the transmit index to the next available element
        
        Index_p = Index_c;  % set parent index to current element
        Index_c = Index_e;  % set current index to new element
        
        Ang = atan2d(v_r(2), v_r(1));  % transmit ray angle
        p_r = BT(Index_p).power*R_eff;  % transmitted power
        
        BT(Index_c).pos = p_i;  % update current node position with intersection
        BT(Index_c).ang = Ang;  % update current node angle with calculated angle
        BT(Index_c).power = p_r;  % update current node power
        BT(Index_c).p_index = Index_p;  % update current node power
        BT(Index_c).b_index = Index_b;  % update the boundary collision
        BT(Index_c).normal = pn;  % update the normal line
        
        Index_e = Index_e + 1;  % update the next available element index
        
        TR_depth = TR_depth + 1;  % update the depth of ray tracing
        
    else  % in case of unknown state
        if(DEBUG == true)  % if debugging is enabled
            disp('undefined state!');  % print the error
        end
    end
    
    if(TR_depth > max_rt_d)  || (BT(Index_c).power/p_s < valid_ratio)  % check for maximum depth of ray tracing or minimum ratio of power
        BT(Index_c).t_index = 0;  % set all the nodes to termination
        BT(Index_c).r_index = 0;  % set all the nodes to termination
        BT(Index_c).p_index = Index_p;  % update current node power
        
        if (Index_p ~= 0)  % if the current node is not the root
            Index_c = BT(Index_c).p_index;  % update the parent index to the node above
            TR_depth = TR_depth - 1;  % update the depth of ray tracing
        end
        
%         Index_p = max(1, BT(Index_c).p_index);  % get the parent index
        Index_b = BT(Index_c).b_index;  % update the touched boundary
    end
    
end

RT_Array = BT(1:(Index_e - 1));  % extract the results