function [V_t, V_r, R_s, R_p, P_n, TIRFlag] = RT_SnellsLaw(V_s, Source, V_b, Boundary, Pos_i, Norm_l, Tol)

DEBUG = false;  % no debugging message

x_s = Source(1);  % source x position (m)
y_s = Source(2);  % source y position (m)
t_s = Source(3);  % propagation angle (deg)
n_s = Source(4);  % source refractive index

x_b_0 = Boundary(1);  % boundary x start (m) - CW direction
y_b_0 = Boundary(2);  % boundary y start (m)
x_b_1 = Boundary(3);  % boundary x stop (m)
y_b_1 = Boundary(4);  % boundary y stop (m)
n_b = Boundary(5); % boundary refractive index

X_i = Pos_i(1);  % intersection x position
Y_i = Pos_i(2);  % intersection x position

ti_ = acosd( dot(V_s, V_b) / ( norm(V_s) * norm(V_b) ) );  % incindent angle (deg)
td = acosd( dot([1, 0], V_b) / ( norm(V_b) ) );  % direction of the boundary

if td < 90  % the choice of (x0, y0) and (x1, y1) is correct
    dFlipBoundary = false;  % boundary vector is straight
    tb = atan2d( +(y_b_1 - y_b_0), +(x_b_1 - x_b_0) );  % boundary angle (deg)
    x_p = x_b_0;  % pivot x coordinate
    y_p = y_b_0;  % pivot x coordinate
    
else  % the choice of (x0, y0) and (x1, y1) is wrong
    dFlipBoundary = true;  % boundary vector is flipped
    tb = atan2d( -(y_b_1 - y_b_0), -(x_b_1 - x_b_0) );  % boundary angle (deg)
    x_p = x_b_1;  % pivot x coordinate
    y_p = y_b_1;  % pivot x coordinate
    
end

if ti_ < 90  % check if the incident angle is in range
    ti = ti_;  % set the angle
else  % check if the incident angle is out of range
    ti = 180 - ti_;  % set the angle
end

r00 = +cosd(-tb); % rotation coefficients
r01 = -sind(-tb); % from http://www.euclideanspace.com/maths/geometry/affine/aroundPoint/matrix2d/
r10 = +sind(-tb);
r11 = +cosd(-tb);

xs_ = r00*x_s + r01*y_s + x_p - r00*x_p - r01*y_p;  % rotated x position of source
ys_ = r10*x_s + r11*y_s + y_p - r10*x_p - r11*y_p;  % rotated y position of source

xi_ = r00*X_i + r01*Y_i + x_p - r00*x_p - r01*y_p;  % rotated x position of incidence
yi_ = r10*X_i + r11*Y_i + y_p - r10*x_p - r11*y_p;  % rotated y position of incidence

ax_i_ = xi_ - xs_;  % intersection x vector - local coordinate
ay_i_ = yi_ - ys_;  % intersection y vector - local coordinate

if ax_i_ < 0  % check for positive x propagation - local coordinate
    ti_x_lc_dir = -1;  % positive propagation
else  % check for negative x propagation
    ti_x_lc_dir = +1;  % negative propagation
end
if ay_i_ < 0  % check for positive y propagation - local coordinate
    ti_y_lc_dir = -1;  % positive propagation
else  % check for negative y propagation
    ti_y_lc_dir = +1;  % negative propagation
end

t1 = 90 - ti;  % angle in medium 1

sin_t2 = (n_s/n_b) * sind(t1);  % sin (angle) in medium 2
t2 = asind(sin_t2);  % transmit angle in medium 2
t2_ = 90 - t2;  % transmit angle in medium 2

if(ti_x_lc_dir >= 0) && (ti_y_lc_dir >= 0) && (dFlipBoundary == false) % moving to +x, +y, Quarter 1, no flipping
    ax_t_ = +cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = +sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = +cosd(ti); % reflection x vector - local coordinate
    ay_r_ = -sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir < 0) && (ti_y_lc_dir >= 0) && (dFlipBoundary == false) % moving to -x, +y, Quarter 2, no flipping
    ax_t_ = -cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = +sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = -cosd(ti); % reflection x vector - local coordinate
    ay_r_ = -sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir < 0) && (ti_y_lc_dir < 0) && (dFlipBoundary == false) % moving to -x, -y, Quarter 3, no flipping
    ax_t_ = -cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = -sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = -cosd(ti); % reflection x vector - local coordinate
    ay_r_ = +sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir >= 0) && (ti_y_lc_dir < 0) && (dFlipBoundary == false) % moving to +x, -y, Quarter 4, no flipping
    ax_t_ = +cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = -sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = +cosd(ti); % reflection x vector - local coordinate
    ay_r_ = +sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir >= 0) && (ti_y_lc_dir >= 0) && (dFlipBoundary == true) % moving to +x, +y, Quarter 1, flipping
    ax_t_ = +cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = +sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = +cosd(ti); % reflection x vector - local coordinate
    ay_r_ = -sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir < 0) && (ti_y_lc_dir >= 0) && (dFlipBoundary == true) % moving to -x, +y, Quarter 2, flipping
    ax_t_ = -cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = +sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = -cosd(ti); % reflection x vector - local coordinate
    ay_r_ = -sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir < 0) && (ti_y_lc_dir < 0) && (dFlipBoundary == true) % moving to -x, -y, Quarter 3, flipping
    ax_t_ = -cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = -sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = -cosd(ti); % reflection x vector - local coordinate
    ay_r_ = +sind(ti); % reflection y vector - local coordinate
    
elseif(ti_x_lc_dir >= 0) && (ti_y_lc_dir < 0) && (dFlipBoundary == true) % moving to +x, -y, Quarter 4, flipping
    ax_t_ = +cosd(t2_); % transmit x vector - local coordinate
    ay_t_ = -sind(t2_); % transmit y vector - local coordinate
    
    ax_r_ = +cosd(ti); % reflection x vector - local coordinate
    ay_r_ = +sind(ti); % reflection y vector - local coordinate
end

if((sin_t2 - 1) < Tol)  % no total internal reflection
    TIRFlag = false;  % no TIR
    if(DEBUG == true)
        disp('no total internal reflection!');  % print the message
    end
else  % total intertnal reflection
    TIRFlag = true;  % TIR
    if(DEBUG == true)
        disp('total internal reflection!');  % print the message
    end
end

xt_ = xi_ + ax_t_;  % transmit x position
yt_ = yi_ + ay_t_;  % transmit y position

xr_ = xi_ + ax_r_;  % reflection x position
yr_ = yi_ + ay_r_;  % reflection y position

r00 = +cosd(tb); % rotation coefficients
r01 = -sind(tb); % from http://www.euclideanspace.com/maths/geometry/affine/aroundPoint/matrix2d/
r10 = +sind(tb);
r11 = +cosd(tb);

xt = r00*xt_ + r01*yt_ + x_p - r00*x_p - r01*y_p;  % rotated x position of transmittance
yt = r10*xt_ + r11*yt_ + y_p - r10*x_p - r11*y_p;  % rotated y position of transmittance

xr = r00*xr_ + r01*yr_ + x_p - r00*x_p - r01*y_p;  % rotated x position of reflectance
yr = r10*xr_ + r11*yr_ + y_p - r10*x_p - r11*y_p;  % rotated y position of reflectance

xn_0 = r00*xi_ + r01*(yi_ - Norm_l) + x_p - r00*x_p - r01*y_p;  % rotated x position of normal start
yn_0 = r10*xi_ + r11*(yi_ - Norm_l) + y_p - r10*x_p - r11*y_p;  % rotated y position of normal start

xn_1 = r00*xi_ + r01*(yi_ + Norm_l) + x_p - r00*x_p - r01*y_p;  % rotated x position of normal start
yn_1 = r10*xi_ + r11*(yi_ + Norm_l) + y_p - r10*x_p - r11*y_p;  % rotated y position of normal start

P_n = [xn_0, yn_0, xn_1, yn_1];  % normal geometry

ax_t = xt - X_i; % transmit x vector - global coordinate
ay_t = yt - Y_i; % transmit y vector - global coordinate

V_t = [ax_t, ay_t];  % transmittance vector

ax_r = xr - X_i; % reflection x vector - global coordinate
ay_r = yr - Y_i; % reflection y vector - global coordinate

V_r = [ax_r, ay_r];  % reflectance vector

R_s = ( abs(n_s*cosd(t1) - n_b*cosd(t2)) / abs(n_s*cosd(t1) + n_b*cosd(t2)) )^2;  % reflectance = reflectivity = power reflection coefficient, s polarization
R_p = ( abs(n_s*cosd(t2) - n_b*cosd(t1)) / abs(n_s*cosd(t2) + n_b*cosd(t1)) )^2;  % reflectance = reflectivity = power reflection coefficient, p polarization
