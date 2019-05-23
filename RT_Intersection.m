function [V_s, V_b, P_i, Dist, Flags] = RT_Intersection(Ray, Line, Tol)

DEBUG = false;  % no debugging message

x_ray = Ray(1);  % incoming ray x position (m)
y_ray = Ray(2);  % incoming ray y position (m)
t_ray = Ray(3);  % propagation angle (deg)

x0 = Line(1);  % line x start (m) - CW direction
y0 = Line(2);  % line y start (m)
x1 = Line(3);  % line x stop (m)
y1 = Line(4);  % line y stop (m)

ms = tand(t_ray); % source ray line slope
cnt_s = y_ray - tand(t_ray)*x_ray;  % source ray line constant
ax_s = cosd(t_ray);  % source x vector
ay_s = sind(t_ray);  % source y vector
V_s = [ax_s, ay_s];  % source vector

mb = (y1 - y0)/(x1 - x0 + eps);  % boundary slope
cnt_b = y1 - mb*x1;  % boundary line constant
ax_b = x1 - x0;  % boundary x vector
ay_b = y1 - y0;  % boundary y vector
V_b = [ax_b, ay_b];  % boundary vector

xi = (cnt_b - cnt_s)/(ms - mb + eps);  % intersection x coordinate
yi1 = ms*xi + cnt_s;  % yi #1 calculation
yi2 = mb*xi + cnt_b;  % yi #2 calculation
if ( ~isinf(yi1) && ~isnan(yi1) )  % check if yi1 is valid
    yi = yi1;  % yi1 is valid
else  % check if yi1 in not valid
    yi = yi2;  % yi2 is assumed to be valid
end

P_i = [xi, yi];  % intersect point

Dist = sqrt( (x_ray - xi)^2 + (y_ray - yi)^2 );  % source to insertion point distance

ax_i = xi - x_ray;  % intersection x vector
ay_i = yi - y_ray;  % intersection y vector
v_i = [ax_i, ay_i];  % intersection vector

cos_t_ray_i = dot(V_s, v_i) / ( norm(V_s) * norm(v_i) + eps);  % cos of angle between the ray and incident vector

if(abs(cos_t_ray_i - 1.0) < Tol)
    dirFlag = true;  % going towards the boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('going towards the boundary!');  % print the message
    end
else
    dirFlag = false;  % coming from the boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('getting away from the boundary!');  % print the message
    end
end

x_i_f = ((max(x0, x1) + Tol) >= xi) && ((min(x0, x1) - Tol) <= xi);   % check if x_i is whithin the boundary range
y_i_f = ((max(y0, y1) + Tol) >= yi) && ((min(y0, y1) - Tol) <= yi);   % check if y_i is whithin the boundary range

if((x_i_f == true) && (y_i_f == true))
    hitFlag = true;  % ray hitting boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('hitting the boundary!');  % print the message
    end
else
    hitFlag = false;  % ray missing boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('missing the boundary!');  % print the message
    end
end

v_d = [V_b(2), -V_b(1)];  % centrifugal vector

cos_d_i = dot(v_d, v_i) / ( norm(v_d) * norm(v_i) + eps);  % cos of angle between the ray and centrifugal vector

if(cos_d_i <= 0)
    sideFlag = true;  % the ray is inside the boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('inside the boundary!');  % print the message
    end
else
    sideFlag = false;  % the ray is outside the boundary
    if(DEBUG == true)  % if debugging is enabled
        disp('outside the boundary!');  % print the message
    end
end

Flags = [dirFlag, hitFlag, sideFlag];  % intersection flags

end