function [Ray, Normal] = RT_RayPlotPreparation(RT_Array, Farfield, Tol)

RT_Array_len = length(RT_Array);  % ray tracing array length

x_f = Farfield(1);  % farfield x centre
y_f = Farfield(2);  % farfield y centre
R_f = Farfield(3);  % farfield radius

R_f_def = 2;  % default farfield radius

X_ray_temp = zeros(1,2*RT_Array_len);  % temporary ray x position
Y_ray_temp = zeros(1,2*RT_Array_len);  % temporary ray y position
U_ray_temp = zeros(1,2*RT_Array_len);  % temporary ray x direction vector
V_ray_temp = zeros(1,2*RT_Array_len);  % temporary ray y direction vector
C_ray_temp = zeros(1,2*RT_Array_len);  % temporary ray color coefficient

Index_ray = 1;  % current index of rays

for Index = 1:RT_Array_len  % go though the nodes to draw arrows
    Index_t = RT_Array(Index).t_index;  % next transmit element
    Index_r = RT_Array(Index).r_index;  % next reflect element
    
    x_0 =  RT_Array(Index).pos(1);  % x position
    y_0 =  RT_Array(Index).pos(2);  % y position
    t_0 =  RT_Array(Index).ang;  % line angle
    
    if( (Index_t == 0) || (Index_t <= -2) ) &&...
            ( (Index_r == 0) || (Index_r <= -2) ) % it is a termination or invalid node or tir
        % do nothing and pass
    elseif( (Index_t == 0) || (Index_t <= -2))  % it is a termination or invalid node
        x_1 =  RT_Array(Index_r).pos(1);  % x end
        y_1 =  RT_Array(Index_r).pos(2);  % y end
        u = (x_1 - x_0);  % x direction vector
        v = (y_1 - y_0);  % y direction vector
        
        X_ray_temp(Index_ray) = x_0;  % update the arrow x
        Y_ray_temp(Index_ray) = y_0;  % update the arrow y
        U_ray_temp(Index_ray) = v;  % update the arrow x direction vector
        V_ray_temp(Index_ray) = u;  % update the arrow y direction vector
        
        C_ray_temp(Index_ray) = RT_Array(Index).power;  % update the power coefficeint
        
        Index_ray = Index_ray + 1;  % update the current index of rays
        
    elseif (Index_t == -1)  % it is a farfield
        
        cx_1 = (x_0 - x_f)*cosd(t_0);  % dummy param x 1
        cy_1 = (y_0 - y_f)*sind(t_0);  % dummy param y 1
        
        cx_2 = (x_0 - x_f)^2;  % dummy param x 2
        cy_2 = (y_0 - y_f)^2;  % dummy param y 2
        
        a_p = 1;  % quadratic equation a parameter
        b_p = 2*(cx_1 + cy_1);  % quadratic equation b parameter
        c_p = cx_2 + cy_2 - R_f^2;  % quadratic equation c parameter
        
        
        R_f_ray = (-b_p + sqrt(b_p^2 - 4*a_p*c_p))/(2*a_p);  % farfield of ray
        
        if (abs(imag(R_f_ray)) > Tol) || (real(R_f_ray) < 0)  % if the calculated radius is invalid
            R_f_ray = R_f_def;  % set the radius to default farfield radius
        end
        
        x_1 = x_0 + R_f_ray*cosd(t_0);  % x end at farfield
        y_1 = y_0 + R_f_ray*sind(t_0);  % y end at farfield

        u = (x_1 - x_0);  % x direction vector
        v = (y_1 - y_0);  % y direction vector
        
        X_ray_temp(Index_ray) = x_0;  % update the arrow x
        Y_ray_temp(Index_ray) = y_0;  % update the arrow y
        U_ray_temp(Index_ray) = v;  % update the arrow x direction vector
        V_ray_temp(Index_ray) = u;  % update the arrow y direction vector
        
        C_ray_temp(Index_ray) = RT_Array(Index).power;  % update the power coefficeint
        
        Index_ray = Index_ray + 1;  % update the current index of rays
        
    else  % it is a ray
        x_1 =  RT_Array(Index_t).pos(1);  % x end
        y_1 =  RT_Array(Index_t).pos(2);  % y end
        u = (x_1 - x_0);  % x direction vector
        v = (y_1 - y_0);  % y direction vector
        
        X_ray_temp(Index_ray) = x_0;  % update the arrow x
        Y_ray_temp(Index_ray) = y_0;  % update the arrow y
        U_ray_temp(Index_ray) = v;  % update the arrow x direction vector
        V_ray_temp(Index_ray) = u;  % update the arrow y direction vector
        
        C_ray_temp(Index_ray) = RT_Array(Index).power;  % update the power coefficeint
        
        Index_ray = Index_ray + 1;  % update the current index of rays
        
    end
end

Index_ray = Index_ray - 1;  % remove the last item

X_ray = X_ray_temp(1:Index_ray);  % cut the extra x locations
Y_ray = Y_ray_temp(1:Index_ray);  % cut the extra y locations
U_ray = U_ray_temp(1:Index_ray);  % cut the extra x direction vectors
V_ray = V_ray_temp(1:Index_ray);  % cut the extra y direction vectors
C_ray = C_ray_temp(1:Index_ray);  % cut the extra color coeefcients

Ray = struct('x', X_ray, 'y', Y_ray, 'u', U_ray, 'v', V_ray, 'c', C_ray);  % create Ray struct

X_norm_temp = zeros(1,2*RT_Array_len);  % temporary normal x position
Y_norm_temp = zeros(1,2*RT_Array_len);  % temporary normal y position
U_norm_temp = zeros(1,2*RT_Array_len);  % temporary normal x direction vector
V_norm_temp = zeros(1,2*RT_Array_len);  % temporary normal y direction vector

Index_normal = 1;  % current index of normals

for Index = 2:RT_Array_len  % go through the normals to draw normal lines
    normal = RT_Array(Index).normal;  % normal line
    
    X_norm_temp(Index_normal) = normal(1);  % nomal x start
    Y_norm_temp(Index_normal) = normal(2);  % nomal y start
    U_norm_temp(Index_normal) = normal(3) - normal(1);  % nomal x direction vector
    V_norm_temp(Index_normal) = normal(4) - normal(2);  % nomal y direction vector
    Index_normal = Index_normal + 1;  % update the current index of normals
    
end

Index_normal = Index_normal - 1;  % remove the last item

X_norm = X_norm_temp(1:Index_normal);  % cut the extra x locations
Y_norm = Y_norm_temp(1:Index_normal);  % cut the extra y locations
U_norm = U_norm_temp(1:Index_normal);  % cut the extra x direction vectors
V_norm = V_norm_temp(1:Index_normal);  % cut the extra y direction vectors

Normal = struct('x', X_norm, 'y', Y_norm, 'u', U_norm, 'v', V_norm);  % create Normal struct