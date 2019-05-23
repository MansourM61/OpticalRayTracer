function Farfield_Rays = RT_EstimateFarfield(RT_Array, Farfield, Tol)

x_f = Farfield(1);  % farfield x centre
y_f = Farfield(2);  % farfield y centre
R_f = Farfield(3);  % farfield radius
Res = Farfield(4);  % farfield resolution
t_ref = Farfield(5);  % farfield reference angle (Deg)

R_f_def = 2;  % default farfield radius

NoR = length(RT_Array);  % number of rays

Farfield_Flag = zeros(1, NoR);  % farfield flag array

for Index = 1:NoR  % go through all the rays
    if (RT_Array(Index).t_index == -1) && (RT_Array(Index).r_index == -1)  % if node is farfield
        Farfield_Flag(Index) = 1;  % set the farfield flag
    end  
end

if (sum(Farfield_Flag) < 1)  % no farfield ray exists
    Farfield_Rays = [];  % set the return matrix to empty
    return;  % return to the caller function
end

Farfield_Array = RT_Array(find(Farfield_Flag == 1));  % create farfield array

NoF = length(Farfield_Array);  % number of farfield array

Angle = linspace(0, 360, Res);  % full space angle

Farfield_index = zeros(1, NoF);  % farfield index array
Farfield_angle = zeros(1, NoF);  % farfield angle array
Farfield_distance = zeros(1, NoF);  % farfield distance array
Farfield_power = zeros(1, NoF);  % farfield distance array

for Index = 1:NoF  % go though the farfield arrayfun

    x_p =  Farfield_Array(Index).pos(1);  % x position
    y_p =  Farfield_Array(Index).pos(2);  % y position
    t_p =  Farfield_Array(Index).ang;  % line angle
         
    cx_1 = (x_p - x_f)*cosd(t_p);  % dummy param x 1
    cy_1 = (y_p - y_f)*sind(t_p);  % dummy param y 1
    
    cx_2 = (x_p - x_f)^2;  % dummy param x 2
    cy_2 = (y_p - y_f)^2;  % dummy param y 2
    
    a_p = 1;  % quadratic equation a parameter
    b_p = 2*(cx_1 + cy_1);  % quadratic equation b parameter
    c_p = cx_2 + cy_2 - R_f^2;  % quadratic equation c parameter    
    
    R_f_ray = (-b_p + sqrt(b_p^2 - 4*a_p*c_p))/(2*a_p);  % farfield of ray
    
    if (abs(imag(R_f_ray)) > Tol) || (real(R_f_ray) < 0)  % if the calculated radius is invalid
        R_f_ray = R_f_def;  % set the radius to default farfield radius
    end
    
    x_r_f = x_p + R_f_ray*cosd(t_p);  % ray x at farfield
    y_r_f = y_p + R_f_ray*sind(t_p);  % ray y at farfield
    
    V_p_x = x_p - x_f;  % centre to point x vector
    V_p_y = y_p - y_f;  % centre to point y vector
    
    V_f_x = x_r_f - x_p;  % point to farfield x vector
    V_f_y = y_r_f - y_p;  % point to farfield y vector
    
    V_r_f_x = V_p_x + V_f_x;  % centre to farfield x vector
    V_r_f_y = V_p_y + V_f_y;  % centre to farfield y vector
    
    t_r_f = atan2d(V_r_f_y, V_r_f_x);  % point to farfield angle
    t_r_f = t_r_f - t_ref;  %  set the angle based on the reference
    if (t_r_f < 0)  % if the angle is from -180 to 0
        t_r_f = t_r_f + 360;  % make the angle between 180 to 360
    end

    [~, Index_min] = min(abs(Angle - t_r_f));  % find the closest angle
    Farfield_index(Index) = Index_min;  % update farfield index array 

    angle_f = Angle(Index_min);  % pick the  closest angle
    Farfield_angle(Index) = angle_f;  % update farfield angle array 
    
    dist_f = norm([V_r_f_x, V_r_f_y]);  % point to farfield distance
    Farfield_distance(Index) = dist_f;  % update farfield distance array 
    
    power_f = Farfield_Array(Index).power;  % farfield power
    Farfield_power(Index) = power_f;  % update farfield distance array 

end

Farfield_Rays = struct('index', Farfield_index, 'angle', Farfield_angle,...
    'distance', Farfield_distance, 'power', Farfield_power);  % create the farfield measurements array
