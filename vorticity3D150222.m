

clear all;
close all;
beep off;
clc;

%PARAMETERS
% filename = ""matlab_pred_linear_ML_2"";
filename = "hpc/ML200_5sl/region1predicted";
% filenames = ["matlab_right_true", "hpc/ML600_5sl/region1predicted", "hpc/LIML600_5sl/region1predicted"];hpc/ML600_5sl/region1predicted
% filename = "hpc/ML600_5sl/region1predicted";
Dimensions = [256 128 160]; % Dimensions of 3d volume
surface_num = 1; %number of iso-surfaces used
surface_min = 1.5; %surface_min * average = value for min voticity surface
surface_max = 2; %surface_max * average = value for max voticity surface. Max used when only one surface shown



% Volume dimensions
m_nQueryPoints = Dimensions;

fprintf('Plotting %ix%ix%i grid points\n', Dimensions(1:3));


uvw3D_field = load(filename,'uvw3D');

field = uvw3D_field.uvw3D(1,:,:,:,:);

% Get velocity gradient
np = py.importlib.import_module('numpy');
[dudx, dudy, dudz] = gradient(field(1,:,:,:,1));
[dvdx, dvdy, dvdz] = gradient(field(1,:,:,:,2));
[dwdx, dwdy, dwdz] = gradient(field(1,:,:,:,3));
% 

m_result9 = [reshape(dudx,1,[]); reshape(dudy,1,[]); reshape(dudz,1,[]); reshape(dvdx,1,[]); reshape(dvdy,1,[]); reshape(dvdz,1,[]); reshape(dwdx,1,[]); reshape(dwdy,1,[]); reshape(dwdz,1,[])];


m_vorticity3 = calculateVorticity(m_result9);    
[w_x, w_y, w_z, ISO] = parseVector(m_vorticity3, m_nQueryPoints);

x = linspace(0, Dimensions(1), Dimensions(1));
y = linspace(0, Dimensions(2), Dimensions(2));                
z = linspace(0, Dimensions(3), Dimensions(3));   
[m_X1, m_X2, m_X3] = meshgrid(x, y, z);

x_figure = figure(1);
m_screenSize = get(0, 'ScreenSize');
if m_screenSize(3) > m_screenSize(4)
    max = m_screenSize(4);
else
    max = m_screenSize(3);
end

set(x_figure, 'Position', [0 0 round(max*1.1) max]);

hold on;
avg = mean(ISO(ISO>0));

ISO = permute(ISO, [2 1 3]); 
x_patch = drawIsoPatch(ISO, m_X1, m_X2, m_X3, surface_min*avg, surface_max*avg, surface_num); %%%plot less for more complex flow     surface_max*avg   0.0572 surface_max*avg
surface_min*avg
surface_max*avg

% a_patch = patch([1 1 256 256], [1 128 128 1], [16 16 16 16], [1 0 0]);
% b_patch = patch([1 1 256 256], [1 128 128 1], [48 48 48 48], [1 0 0]);
% c_patch = patch([1 1 256 256], [1 128 128 1], [80 80 80 80], [1 0 0]);
% d_patch = patch([1 1 256 256], [1 128 128 1], [112 112 112 112], [1 0 0]);
% e_patch = patch([1 1 256 256], [1 128 128 1], [144 144 144 144], [1 0 0]);


xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('z', 'FontSize', 12, 'FontWeight', 'bold');
view(3);
view(-60,45)
%axis vis3d;
camlight;
lighting phong;
axis equal;
axis([0 256 0 128 0 160]);
set(gca, 'FontSize', 11);
xticks([0 32 64 96 128 160 192 224 256])
yticks([0 32 64 96 128])
zticks([0 32 64 96 128 160])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
grid;
% a_patch.FaceVertexAlphaData = 0.2;
% a_patch.FaceAlpha = 'flat' ; 
% x_patch.FaceVertexAlphaData = 0.5;
% x_patch.FaceAlpha = 'flat' ;
% alpha(a_patch,0.2);   %%transparancy
% alpha(b_patch,0.2);   %%transparancy
% alpha(c_patch,0.3);   %%transparancy
% alpha(d_patch,0.2);   %%transparancy
% alpha(e_patch,0.2);   %%transparancy
alpha(x_patch,0.5);   %%transparancy
title(sprintf('%s','Iso-Vorticity Surfaces'), 'FontSize', 13, 'FontWeight', 'bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functions



function [u, v, w, mag] = parseVector(m_results, m_nQueryPoints)

    magv = sqrt(m_results(1,:) .* m_results(1,:) + ...
                m_results(2,:) .* m_results(2,:) + ...
                m_results(3,:) .* m_results(3,:));

    mag = reshape(magv, m_nQueryPoints);
    u = reshape(m_results(1,:), m_nQueryPoints);
    v = reshape(m_results(2,:), m_nQueryPoints);
    w = reshape(m_results(3,:), m_nQueryPoints);

end

function result = calculateVorticity(gradient)
    result(1,:) = gradient(8,:) - gradient(6,:);
    result(2,:) = gradient(3,:) - gradient(7,:);
    result(3,:) = gradient(4,:) - gradient(2,:);
end



function x_patch = drawIsoPatch(scalar, m_X1, m_X2, m_X3, startf, endf, npoints)


    m_ones = ones(size(scalar));

    for i = linspace(startf, endf, npoints)
        % Check if this iso level exists
        if ~isempty(scalar(scalar>i))
            [x_faces, x_verts, x_colors] = isosurface(m_X1, m_X2, m_X3, scalar, i, m_ones*i);

            x_patch = patch('Vertices', x_verts, 'Faces', x_faces, ... 
                'FaceVertexCData', x_colors, ...
                'FaceColor','interp', ... 
                'edgecolor', 'none');
        end

    end

    if npoints > 1
        colorbar;
        colormap 'Jet';
        colorbar('FontSize', 12);
    end

end






%
% Turbmat-Tools - a Matlab library for querying, processing and visualizing
% data from the JHU Turbulence Database
%   
% TurbCache, part of Turbmat-Tools
%

%
% Written by:
% 
% Edo Frederix 
% The Johns Hopkins University / Eindhoven University of Technology 
% Department of Mechanical Engineering 
% edofrederix@jhu.edu, edofrederix@gmail.com
%

%
% This file is part of Turbmat-Tools.
% 
% Turbmat-Tools is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
% 
% Turbmat-Tools is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Turbmat-Tools.  If not, see <http://www.gnu.org/licenses/>.
%


%
% ---- Initiate ----
%

