% Ray Tracer
% Written by Mojtaba Mansour Abadi
% 23/05/2018
%
% This code is tested with MATLAB R2016a and OCTAVE 4.4.1.
% This code performs ray tracing and plots rays as well as farfield
% estimation of the rays.
% The optics are specified by several geometries. The geometries are
% approximated by piecewise lines and the intersection of the rays with
% each line piece is calculated.
% The reflection and refraction angles and directions are calculated using
% Snell's law where as the power splitting ratios are obtained by means of
% Fresnel equations.

%% Cleaning Environment
clc;
clear all;
close all;


%% Parameters
% Change the parameters in this section to define geometries and sources,
% etc.

% Change this value according to the background material in the space
n_b = 1.45;  % background refractive index


% Define all closed geometries in this part.
% Each geometry is composed of several pieces. Each piece is a parametric
% function handle. The geometry pieces must be defined clockwise.
% Geometry 1
geometry_1 = { @(t) [t, +5], [+1, +3];
             @(t) [+3, t], [+5, +1];
             @(t) [t, +1], [+3, +1];
             @(t) [+1, t], [+1, +5];
           };  % geometry 1 definiton
n_g_1 = 5;  % geometry 1 refractive index

% Geometry 2
geometry_2 = { @(t) [t, +5], [-3, -1];
             @(t) [-1, t], [+5, +1];
             @(t) [t, +1], [-1, -3];
             @(t) [-3, t], [+1, +5];
           };  % geometry 2 definiton
n_g_2 = 2;  % geometry 2 refractive index

% Geometry 2
geometry_3 = { @(t) [t, -1], [-3, +3];
             @(t) [+3, t], [-1, -3];
             @(t) [t, -3], [+3, -3];
             @(t) [-3, t], [-3, -1];
           };  % geometry 1 definiton
n_g_3 = 10;  % geometry 1 refractive index


% Define all sources in this section.
% Each source has a cartesdian location as well as propagation angle
% reletive to +x axis. Each source shoots a ray based on the given
% properties carrying a power < 1.0
% Source 1
x_s_1 = 2;  % source 1 x position
y_s_1 = 2;  % source 1 y position
t_s_1 = 50;  % source 1 propagation angle (Deg)
p_s_1 = 1.0;  % source 1 power

% Source 2
x_s_2 = 1;  % source 2 x position
y_s_2 = 0;  % source 2 y position
t_s_2 = 130;  % source 2 propagation angle (Deg)
p_s_2 = 1.0;  % source 2 power

% Source 3
x_s_3 = 2;  % source 3 x position
y_s_3 = 0;  % source 3 y position
t_s_3 = 240;  % source 3 propagation angle (Deg)
p_s_3 = 1.0;  % source 3 power

% Source 4
x_s_4 = ones(1, 5)*-2;  % source 4 x position
y_s_4 = ones(1, 5)*0;  % source 3 y position
t_s_4 = linspace(-30, 30, 5) + 90;  % source 4 propagation angle (Deg)
p_s_4 = ones(1, 5).*cosd(t_s_4 - 90);  % source 4 power

% ray tracing parameters
dx = 0.1;  % x resolution for changing the geometry into line pieces
dy = 0.1;  % y resolution for changing the geometry into line pieces
coeff = 0.9;  % for changing the geometry into line pieces, if smaller line is required, this value specifies the line size reduction
norm_len = 0.2;  % normal line length at the ray-boundary incindence location
TOL = 10*eps;  % numerical tolerance for calculation errors
max_rt_d = 10;  % if the ray is bounced more, the tracing stops
valid_ratio = 0.01;  % if the bounced power to ray power is less, the tracing stops


% fafield parameters
% The parameters in this part define a circle to be the boundary of
% farfield. The radiation pattern is calculated from the rays approaching
% this boundary.
x_f = 0;  % farfield x centre
y_f = 0;  % farfield y centre
R_f = 7;  % farfield radius
Res = 360;  % farfield resolution; the number of angles from 0 to 360 Deg
t_ref = 0;  % farfield reference angle (Deg); reference for farfield angles


%% Initialisation
% In this section the required arrays for ray tracing algorithm is
% generated. Change the corresponding values according to the defined
% geometries or sources.

% Add or remove geometries and corresponding refractive index to be 
% included in the ray tracing.
Geometry = {geometry_1, geometry_2, geometry_3};  % array of all geometries
n_g = [n_g_1, n_g_2, n_g_3];  % geometry refractive index


% Add or remove position, angle and power of sources to be included in the
% ray tracing.
x_s = [x_s_1, x_s_2, x_s_3, x_s_4];  % source x position
y_s = [y_s_1, y_s_2, y_s_3, y_s_4];  % source y position
t_s = [t_s_1, t_s_2, t_s_3, t_s_4];  % source propagation angle (Deg)
p_s = [p_s_1, p_s_2, p_s_3, p_s_4];  % source power

Farfield = [x_f, y_f, R_f, Res, t_ref];  % farfield parameters
delta = [dx, dy];  % resolution in x and y direction
NoG = length(Geometry);  % number of geometries


%% Step 1: Quantise the Geometry

ShapePoints = cell(1, NoG);  % array of shape points

for Index_G = 1:NoG  % go through the geometries
    ShapePoints{Index_G} = RT_GeometryQuantizer(Geometry{Index_G}, delta, coeff);  % quantise the geometry
    NoP = size(ShapePoints, 1);  % number of points
end


%% Step 2: Check Location of the Source

NoS = length(x_s);  % number of sources

n_s = zeros(1, NoS);  % source refractive index array

for Index_S = 1:NoS  % go through the sources

    point = [x_s(Index_S), y_s(Index_S)];  % source coordinate
    
    flag = false;  % initial flag

    for Index_G = 1:NoG  % go through the geometries
        flag = RT_InsideShape(point, ShapePoints{Index_G}, TOL);  % check if the point is inside the shape
        if flag == true
            break;  % stop scanning more geometries
        end
    end

    if(~flag)  % if the source is outside the shape
        n_s(Index_S) = n_b;  % update the source refractive index with background medium
    else   % if the source is inside the shape
        n_s(Index_S) = n_g(Index_G);  % update the source refractive index with geometry medium
    end

end


%% Step 3: Ray Tracing
RT_Array = cell(1, NoS);  % array of rays

for Index_S = 1:NoS  % go through the sources
    point = [x_s(Index_S), y_s(Index_S)];  % source coordinate
    source = [point, t_s(Index_S), p_s(Index_S)];  % source property
    
    rt_param = [norm_len, max_rt_d, valid_ratio];  % ray tracing parameters

    RT_Array{Index_S} = RT_RayTracer(source, ShapePoints, n_b, n_g, rt_param, TOL);  % perform ray tracing for given source and shape
end


%% Step 4: Calculate Farfield
Farfield_Rays = cell(1, NoS);  % array of rays

Angle = linspace(0, 360, Res);  % full space angle
RadPat = zeros(1, Res);  % initial radiation pattern

for Index_S = 1:NoS  % go through the sources
    Farfield_Rays{Index_S} = RT_EstimateFarfield(RT_Array{Index_S}, Farfield, TOL);  % calculate the farfield
    
    if isempty(Farfield_Rays{Index_S})  % if there is no farfield
        continue;  % skip the loop
    end
    
    Ang_f = Farfield_Rays{Index_S}.index;  % index of the farfield angles
    Power_f = Farfield_Rays{Index_S}.power;  % index of the farfield powers
    
    RadPat(Ang_f) = RadPat(Ang_f) + Power_f;  % accumulate the power
    
end


%% Step 5: Ray Tracing Presentation
F_RT_H = figure;  % create a new figure
hold on;  % hold the drawings
box on; % create the box around the figure

RT_RayPlotter(gca, [], ShapePoints, n_g, [], [], 'Geometry', TOL);  % plot the geometries

for Index_S = 1:NoS  % go through the sources
    source = [x_s(Index_S), y_s(Index_S)];  % source coordinate

    [Ray_H, Arrow_H, Normal_H] = RT_RayPlotter(gca, source, [], [],...
    RT_Array{Index_S}, Farfield, 'Ray', TOL);  % plot the ray tracing
end

source = [x_s(Index_S), y_s(Index_S)];  % source coordinate

RT_RayPlotter(gca, [], [], [], [], Farfield, 'Farfield', TOL);  % plot the farfield

xlabel('x');  % x axis label
ylabel('y');  % y axis label
title('Ray Traces');  % set the title of figure
axis equal;  % set the aspect ratio of figure to 1

MakeitPretty(F_RT_H, [10, 9], ['L', 'L'], [12, 0.5, 5, 10], 'Ray_Tracing');  % save ray tracing plot


%% Step 6: Farfield Presentation
F_FF_H = figure;  % create a new figure
hold on;  % hold the drawings
box on; % create the box around the figure

MarkerStyle = {'', 'o', '*', 's', '^', 'h', 'x', '+', 'd', 'v', '<', '>', 'p'};  % define marker styles

MarkerPlot(Angle, RadPat, 'b', '-', MarkerStyle{2}, 10);

xlabel('\phi (Deg)');  % x axis label
ylabel('P');  % y axis label
title('Farfield Radiation Pattern');  % set the title of figure
axis([0, 360, 0, max(RadPat)]);  % set the axis limits
grid on;  % switch on the grids

MakeitPretty(F_FF_H, [10, 9], ['L', 'L'], [12, 1, 5, 10], 'Farfield_RadPat');  % save farfield radiation pattern plot