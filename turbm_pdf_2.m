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

clear all;
close all;
beep off;
clc;

TT = TurbTools(0);


%PARAMETERS
% filename = "matlab_pred_400_lap";
% filename = "matlab_test_2";
%filenames = ["matlab_right_true", "matlab_pred_400", "matlab_pred_linear_ML_2"];
filenames = ["matlab_right_true", "hpc/regiona"];%%%, "hpc/LIML200_5sl/regionbpredicted"];
Dimensions = [256 128 160]; % Dimensions of 3d volume
steps = 60;   % Number of steps in the PDF histogram
i_nondim = 0; % nondimensionalized axis?
i = 1; %%% velocity component

label1 = 'DNS region 1';
label2 = 'DNS region A';
%label3 = 'LI + 3D CNN region b';

% Volume dimensions
m_nQueryPoints = Dimensions;

fprintf('Plotting %ix%ix%i grid points\n', Dimensions(1:3));


m_colors = TT.createColors(3);

%
% ---- Calculate Velocity PDF ----
%

x_figure = TT.startFigure(1);

for i = 1:3
    mins =[];
    maxs =[];
    for ind = 1:length(filenames)
        uvw3D_field = load(filenames(ind),'uvw3D');
        field = uvw3D_field.uvw3D(1,:,:,:,:);
        m_result4.('fieldname'+string(ind)) = [reshape(field(:,:,:,:,1),1,[]); reshape(field(:,:,:,:,2),1,[]); reshape(field(:,:,:,:,3),1,[])];
        mins = [mins; min(m_result4.('fieldname'+string(ind))(i,:))];
        maxs = [maxs; max(m_result4.('fieldname'+string(ind))(i,:))];
    end

    mn = min(mins);
    mx = max(maxs); 
    x = linspace(mn, mx, steps); % output x axis 
    space = ones(1, steps).*((mx-mn)/steps);
    for ind = 1:length(filenames)
        edges = linspace(mn, mx, steps+1);
        PDF = histcounts(m_result4.('fieldname'+string(ind))(i,:), edges,'Normalization','pdf');  % hist function
        % set the surface under the curve equal to 1
    %     surf = sum(PDF.*space);
        y.('fieldname'+string(ind)+string(i)) = PDF;%%%./surf;
    end
end




% u
subplot(3,1,1);
x_bar = zeros(1,3);
x_plot = zeros(1,3);
x_bar(1) = bar(x, y.('fieldname11'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'b', 'Edgecolor', 'b'); hold on;

x_bar(2) = bar(x, y.('fieldname21'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'g', 'Edgecolor', 'g'); hold on;

%x_bar(3) = bar(x, y.('fieldname31'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'r', 'Edgecolor', 'r');

% Style figure
grid;
ylabel('Pdf(v_i)', 'FontSize', 12, 'FontWeight', 'bold');
if i_nondim; xlabel('v_i/{\sigma_{v_i}}', 'FontSize', 12, 'FontWeight', 'bold'); else xlabel('V_i', 'FontSize', 12, 'FontWeight', 'bold'); end
%legend([x_bar, x_plot(1)], 'v_x', 'v_y', 'v_z', 'outline');
legend(label1, label2, 'Location', 'Eastoutside'); %label3,
title(sprintf('PDF of velocity u'), 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');

% v
subplot(3,1,2);
x_bar = zeros(1,3);
x_plot = zeros(1,3);

x_bar(1) = bar(x, y.('fieldname12'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'b', 'Edgecolor', 'b'); hold on;

x_bar(2) = bar(x, y.('fieldname22'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'g', 'Edgecolor', 'g'); hold on;

%x_bar(3) = bar(x, y.('fieldname32'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'r', 'Edgecolor', 'r');

% Style figure
grid;
ylabel('Pdf(v_i)', 'FontSize', 12, 'FontWeight', 'bold');
if i_nondim; xlabel('v_i/{\sigma_{v_i}}', 'FontSize', 12, 'FontWeight', 'bold'); else xlabel('V_i', 'FontSize', 12, 'FontWeight', 'bold'); end
legend(label1, label2, 'Location', 'Eastoutside'); %label3,
title(sprintf('PDF of velocity v'), 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');

% w
subplot(3,1,3);
x_bar = zeros(1,3);
x_plot = zeros(1,3);

x_bar(1) = bar(x, y.('fieldname13'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'b', 'Edgecolor', 'b'); hold on;
x_bar(2) = bar(x, y.('fieldname23'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'g', 'Edgecolor', 'g'); hold on;
%x_bar(3) = bar(x, y.('fieldname33'), 'facealpha', 0.5, 'edgealpha', 0.5, 'Facecolor', 'r', 'Edgecolor', 'r');

% Style figure
grid;
ylabel('Pdf(v_i)', 'FontSize', 12, 'FontWeight', 'bold');
if i_nondim; xlabel('v_i/{\sigma_{v_i}}', 'FontSize', 12, 'FontWeight', 'bold'); else xlabel('V_i', 'FontSize', 12, 'FontWeight', 'bold'); end
legend(label1, label2, 'Location', 'Eastoutside'); %label3,
title(sprintf('PDF of velocity w'), 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
% TT.makeFigureSquare(x_figure);


% % % % % % % % % 


% % Get velocity gradient
% np = py.importlib.import_module('numpy');
% [dudx, dudy, dudz] = gradient(field(1,:,:,:,1));
% [dvdx, dvdy, dvdz] = gradient(field(1,:,:,:,2));
% [dwdx, dwdy, dwdz] = gradient(field(1,:,:,:,3));
% % 
% %%m_result9 = TT.callDatabase('getVelocityGradient', i_points, m_points, f_time, 0);
% m_result9 = [reshape(dudx,1,[]); reshape(dudy,1,[]); reshape(dudz,1,[]); reshape(dvdx,1,[]); reshape(dvdy,1,[]); reshape(dvdz,1,[]); reshape(dwdx,1,[]); reshape(dwdy,1,[]); reshape(dwdz,1,[])];
% 
% subplot(2,2,4);
% steps = 40;
% 
% % transverse and longitudinal gradient components
% m_tvGrad = m_result9([2 3 4 6 7 8], :);
% m_ltGrad = m_result9([1 5 9], :);
% [m_tv, m_tvGradPDF] = TT.calculatePDF(m_tvGrad, steps, 10, i_nondim, 1);
% [m_lt, m_ltGradPDF] = TT.calculatePDF(m_ltGrad, steps, 10, i_nondim, 1);
% 
% plot(m_lt, log10(m_ltGradPDF), 'Color', m_colors(1,:), 'LineWidth', 1.3); hold on;
% plot(m_tv, log10(m_tvGradPDF), 'Color', m_colors(end,:), 'LineWidth', 1.3);
% plot(m_lt, log10(m_ltGradPDF), 'r.', 'LineWidth', 1.3);
% plot(m_tv, log10(m_tvGradPDF), 'r.', 'LineWidth', 1.3);
% 
% % Style figure
% grid;
% ylabel('Pdf(J_{i,j})', 'FontSize', 12, 'FontWeight', 'bold');
% if i_nondim; xlabel('{J_{i,j}}/{\sigma_{J_{i,j}}}', 'FontSize', 12, 'FontWeight', 'bold'); else xlabel('{J_{i,j}}', 'FontSize', 12, 'FontWeight', 'bold'); end
% legend('Longitudinal J_{i,i}', 'Transverse J_{i,j}, i <> j', 'Location', 'NorthWest');
% title('PDF of J_{i,j} = {\delta}v_i/{\delta}x_j', 'FontSize', 12, 'FontWeight', 'bold');
% set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
% TT.makeFigureSquare(x_figure);
