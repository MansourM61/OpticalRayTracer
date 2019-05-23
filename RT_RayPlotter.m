function [Ray_H, Arrow_H, Normal_H] = RT_RayPlotter(Axe_H, Source, Geometries, n_Geometry, RT_Array, Farfield, drawFlag, Tol)

DEBUG = false;  % no debugging message

ArrowLen = 0.5;  % arrow vector length
HeadSize = 0.5;  % arrow head size
NormalSize = 2;  % normal line size
RayBaseColor = [1, 0, 0];  % ray base color
FarfieldBaseColor = [1, 0, 1];  % farfield base color
GeometryBaseColor = [0, 1, 1];  % Geometry base color


axes(Axe_H);  % se teh ecurrent axes

if strcmp(drawFlag, 'Geometry')  % if function is called for drawing geometries
 
    NoG = length(Geometries);  % number of geometries
    
    n_max = max(n_Geometry);  % maximum refractive index
    n_min = min(n_Geometry);  % minimum refractive index
    
    BaseColor = rgb2hsv(GeometryBaseColor);  % map the RGB color to HSV space
    
    H_value = BaseColor(1);  % hue value of HSV color
    V_Value = BaseColor(3);  % value value of HSV color
    
    Sat_min = 0.25;  % minimum saturation
    Sat_max = 1.00;  % maximum saturation
    
    if abs(n_max - n_min) > Tol  % if materials are not the same
        m_color = (Sat_max - Sat_min)/(n_max - n_min + eps);  % linear gradient of saturation with respect to refractive index
    else  % if materials are the same
        m_color = 0;  % set the saturation gadient to zero
    end
    
    for Index = 1:NoG
        
        S_value = n_Geometry(Index)*m_color + Sat_max - n_max*m_color;  % calculate S value
        
        Color_hsv = [H_value, S_value, V_Value];  % geometry color
        Color = hsv2rgb(Color_hsv);  % map the HSV color to RGB space
        
        fill(Geometries{Index}(:, 1), Geometries{Index}(:, 2), Color);  % draw the shape fill
        plot(Geometries{Index}(:, 1), Geometries{Index}(:, 2), 'k-', 'LineWidth', 2);   % draw the shape boundary

        x_b_0 = min(Geometries{Index}(:, 1));  % x0 bounding box
        x_b_1 = max(Geometries{Index}(:, 1));  % x1 bounding box
        y_b_0 = min(Geometries{Index}(:, 2));  % y0 bounding box
        y_b_1 = max(Geometries{Index}(:, 2));  % y1 bounding box

        plot([x_b_0, x_b_1, x_b_1, x_b_0, x_b_0], [y_b_0, y_b_0, y_b_1, y_b_1, y_b_0], 'b-.');   % draw the shape bounding box
    end
    
elseif strcmp(drawFlag, 'Ray')  % if function is called for drawing rays

    x_s = Source(1);  % source x position
    y_s = Source(2);  % source y position

    plot(x_s, y_s, 'Color', RayBaseColor, 'Marker', 'o', 'LineStyle', '-', 'MarkerFaceColor', 'r', 'MarkerSize', 10);  % draw the source

    [Ray, Normal] = RT_RayPlotPreparation(RT_Array, Farfield(1:3), Tol);  % create the rays and normals arrays

    X_ray = Ray.x;  % ray x locations
    Y_ray = Ray.y;  % ray y locations
    U_ray = Ray.u;  % ray x direction vectors
    V_ray = Ray.v;  % ray y direction vectors
    C_ray = Ray.c;  % ray colors

    Ray_H = zeros(1, length(X_ray));  % rays handle
    Arrow_H = zeros(1, length(X_ray));  % arrows handle

    for Index = 1:length(X_ray)  % go through the nodes
        Color = C_ray(Index)*RayBaseColor;  % update the color
        
        Ray_H(Index) = quiver(X_ray(Index), Y_ray(Index),...
            V_ray(Index), U_ray(Index), 'Color', Color,...
            'AutoScale', 'off', 'ShowArrowHead', 'off');  % plot the rays
        Arrow_H(Index) = quiver(X_ray(Index), Y_ray(Index),...
            V_ray(Index)*ArrowLen, U_ray(Index)*ArrowLen,...
            'Color', Color, 'MaxHeadSize', HeadSize);  % plot the arrows
    end

    X_norm = Normal.x;  % raw normal x locations
    Y_norm = Normal.y;  % raw normal y locations
    U_norm = Normal.u;  % raw normal x direction vectors
    V_norm = Normal.v;  % raw normal y direction vectors

    [X_sorted, I_sorted] = sort(X_norm);  % sort x locations
    Y_sorted = Y_norm(I_sorted);  % sort y locations
    U_sorted = U_norm(I_sorted);  % sort x direction vectors
    V_sorted = V_norm(I_sorted);  % sort y direction vectors

    dx = diff(X_sorted);  % calculate x difference
    dy = diff(Y_sorted);  % calculate y difference
    du = diff(U_sorted);  % calculate u difference
    dv = diff(V_sorted);  % calculate v difference

    I_dup = find((abs(dx) < Tol) & (abs(dy) < Tol) & (abs(du) < Tol) & (abs(dv) < Tol));  % find similar normals

    X_norm(I_dup) = [];  % remove x location duplicates
    Y_norm(I_dup) = [];  % remove y location duplicates
    U_norm(I_dup) = [];  % remove x direction vector duplicates
    V_norm(I_dup) = [];  % remove y direction vector duplicates

    Normal_H = zeros(1, length(X_norm));  % normal handle

    for Index = 1:length(X_norm)  % go through the normals
        Normal_H(Index) = quiver(X_norm(Index), Y_norm(Index),...
            U_norm(Index), V_norm(Index),...
            'k:', 'AutoScale', 'off', 'ShowArrowHead', 'off', 'LineWidth', NormalSize);  % plot the normals
    end

elseif strcmp(drawFlag, 'Farfield')  % if function is called for drawing farfield


    x_f = Farfield(1);  % farfield x centre
    y_f = Farfield(2);  % farfield y centre
    R_f = Farfield(3);  % farfield radius
    Res = Farfield(4);  % farfield resolution
    t_ref = Farfield(5);  % farfield reference angle (Deg)

    t = linspace(0, 360, Res);  % theta angle (Deg)

    x = x_f + R_f*cosd(t);  % x locus of farfield
    y = y_f + R_f*sind(t);  % y locus of farfield

    plot(x, y, 'Color', FarfieldBaseColor, 'LineStyle', '--');  % plot farfield circle

    x_ref = x_f + [0, R_f*cosd(t_ref)];  % reference x locus of farfield
    y_ref = y_f + [0, R_f*sind(t_ref)];  % reference y locus of farfield

    plot(x_ref, y_ref, 'Color', FarfieldBaseColor, 'LineStyle', '--');  % plot reference angle

    plot(x_f, y_f, 'Color', FarfieldBaseColor, 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10);  % plot farfield centre

else  % in case of unknown state
    if(DEBUG == true)  % if debugging is enabled
        print('Unknown state!');  % print the error
    end
end