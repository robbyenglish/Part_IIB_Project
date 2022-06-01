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
% Modified by:
%
% Jason Graham
% The Johns Hopkins University
% Department of Mechanical Engineering
% jgraha8@gmail.com
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


% ---- User input ----
%
% Number of lines
i_lines = 1000;  
domain = [256 128 128];
PLOT = 0;

%create lines  --- 128x128x128cube
            
lines = struct();
len = struct();
i_points = 0;

for l = 1:i_lines
    %select random start point
    p = zeros(2,3);
    ind = 1:3;
    a = floor(rand*3)+1; % pick random axis
    %a = 1; %lines in x only    % use this when resolution not constant across axis
    ind = setxor(ind,a); % find the other two axis
    p(1,a) = 1; % zero chosen axis
    p(1,ind(1)) = ceil(rand*domain(ind(1))); % random point on other two axis
    p(1,ind(2)) = ceil(rand*domain(ind(2)));
    key = char(strcat('line', num2str(l)));
    %select end point
    p(2,:) = p(1,:);
    len.(key) = domain(a);   
    p(2,a) = len.(key);
    i_points = i_points+len.(key);
    %store line in struct

    lines.(key).x = linspace(p(1,1), p(2,1), len.(key));
    lines.(key).y = linspace(p(1,2), p(2,2), len.(key));
    lines.(key).z = linspace(p(1,3), p(2,3), len.(key));
    lines.(key).dir = a;
end

m_points = TT.fillLines(lines);

%filename = "matlab_test_3";
%filenames = ["matlab_right_true","matlab_right_pred","matlab_linear_pred","matlab_pred_linear_ML_2"]%%%,"matlab_pred_400"];

% filenames = ["matlab_right_true", "hpc/ML200_5sl/region1predicted", "hpc/ML600_5sl/region1predicted"];
% colours = ["#0000ff","#00ff00","#008844"];
% leg = ["DNS data","2D-3D CNN 200 snapshots","2D-3D CNN 600 snapshots"];
% 
% filenames = ["matlab_right_true", "hpc/LIML200_5sl/region1predicted", "hpc/LIML600_5sl/region1predicted"];
% colours = ["#0000ff","#ff0000","#ffff00"];
% % leg = ["DNS data","LI + 3D CNN 200 snapshots","LI + 3D CNN 600 snapshots"];
% 
% filenames = ["matlab_right_true", "hpc/ML200_5sl/region1predicted", "hpc/LIML200_5sl/region1predicted"];
% colours = ["#0000ff","#00ff00","#ff0000"];
% leg = ["DNS data","2D-3D CNN 200 snapshots","LI + 3D CNN 200 snapshots"];

% filenames = ["matlab_right_true", "hpc/LIML200_5sl/region1predicted", "hpc/LIML200_9sl/region1predicted"];
% colours = ["#0000ff","#ff0000","#ffff00"];
% leg = ["DNS data","LI + 3D CNN 5 slices","LI + 3D CNN 9 slices"];

filenames = ["matlab_right_true", "hpc/LIML200_5sl/region1predicted", "hpc/LIML200_9sl/region1predicted", "hpc/LIML200_19sl/region1predicted", "hpc/LIML200_39sl/region1predicted"];
filenames = ["matlab_right_true", "hpc/regionb", "hpc/ML200_5sl/regionbpredicted", "hpc/regionx", "hpc/regionz"];
colours1 = [0 0 1];
colours2 = [1 0 0];
colours3 = [0 1 0];
colours4 = [1 0.66 0];
colours5 = [1 1 0];
%leg = ["DNS data","LI + 3D CNN region 1 (trained)","LI + 3D CNN region x","LI + 3D CNN region y","LI + 3D CNN region z"];
leg = ["DNS data","LI + 3D CNN 5 slices","LI + 3D CNN 9 slices","LI + 3D CNN 19 slices","LI + 3D CNN 39 slices"];
leg = ["Region 1 DNS","Region B DNS","Region B CNN output","regionx","regionz"];


for ind = 1:length(filenames)


    uvw3D_field = load(filenames(ind),'uvw3D');

    m_result3 = zeros(3,i_points);

    for i = 1:i_points
      m_result3(:,i) = uvw3D_field.uvw3D(1,m_points(1,i),m_points(2,i),16+m_points(3,i),:);
    end

    %
    % ---- Calculate spectra ----
    %


    keys = fieldnames(lines);
    index = 1;
    vInlineStruct = struct();
    for i = 1:numel(keys)
        key = char(keys(i));
        dir = lines.(key).dir;
        inc = numel(lines.(key).x);
%         s_inlineVel.(key) = vecnorm([m_result3(1, index:(index+inc-1));m_result3(2, index:(index+inc-1));m_result3(3, index:(index+inc-1))],2,1);
%         s_inlineVel.(key) = [m_result3(1, index:(index+inc-1));m_result3(2, index:(index+inc-1));m_result3(3, index:(index+inc-1))];
        s_inlineVel.(key) = m_result3(dir, index:(index+inc-1));
        index = index+inc;
    end
    
%     s_inlineVel = TT.parseLines(m_result3, lines);
    s_inlineVel = TT.calculateZeroMean(s_inlineVel);

    keys = fieldnames(s_inlineVel);
    allx = zeros(1, i_points);
    i_start = 0;
    i_next = 0;
    bins = 1024;
    dft.('fieldname'+string(ind)) = zeros(1,bins);
    pwr.('fieldname'+string(ind)) = zeros(1,bins);
    PXX.('fieldname'+string(ind)) = zeros();

    for i = 1:numel(keys)
        key = char(keys(i));
        x = s_inlineVel.(key);          %signal
        m = length(x);                  %window length
        n = pow2(nextpow2(m));          %transform length
        y = fft(x,bins);                   %DFT
        power = y.*conj(y)/m;           %power of the DFT
        y_pxx = pwelch(x,[],0,1024);


        %collect
        dft.('fieldname'+string(ind)) = dft.('fieldname'+string(ind)) + abs(y);
        pwr.('fieldname'+string(ind)) = pwr.('fieldname'+string(ind)) + abs(power);
        PXX.('fieldname'+string(ind)) = PXX.('fieldname'+string(ind)) + abs(y_pxx);

        i_start = i_next + 1;
        i_next = i_next + m;
        allx(1, i_start:i_next) = x;


    end


    % 
    Pxx.('fieldname'+string(ind)) = pwelch(allx,1024,0,1024);
    
    %average
    dft.('fieldname'+string(ind)) = 2*dft.('fieldname'+string(ind))/numel(keys);

    pwr.('fieldname'+string(ind)) = 2*pwr.('fieldname'+string(ind))/numel(keys);

    PXX.('fieldname'+string(ind)) = PXX.('fieldname'+string(ind))/numel(keys);

    

    output = TT.calculateStatProperties(s_inlineVel);


    fprintf('Mean inline velocity: %1.4f. Variance: %1.4f. Mean squared: %1.4f. Standard deviation: %1.4f\n', output);

    if PLOT
        %
        % ---- Plot all signals ----
        %

        x_figure = TT.startFigure(1);
        subplot(2,10,[2 9]);

        keys = fieldnames(s_inlineVel);

        m_colors = TT.createColors(i_lines);
        for i = 1:numel(keys)
            key = char(keys(i));
            x = linspace(0,2*pi,len.(key));
            plot(x, s_inlineVel.(key), 'LineWidth', 1.3, 'Color', m_colors(i,:));
            hold on;
        end

        % Style figure
        grid;
        title(sprintf('Inline velocity for %i random lines', i_lines), 'FontSize', 12, 'FontWeight', 'bold');
        TT.setFigureAttributes('1d', {'x', 'v'});
        % Adjust axis from default
        xlim([0.0 2*pi]);

        %
        % ---- Plot random lines in space ----
        %

        subplot(2,10,[17 20]);
        keys = fieldnames(lines);
        for i = 1:numel(keys);
            key = char(keys(i));
            plot3(lines.(key).x, lines.(key).y, lines.(key).z, 'Color', m_colors(i,:), 'LineWidth', 1.3);
            hold on;
        end

        title(sprintf('Randomly chosen lines in domain'), 'FontSize', 12, 'FontWeight', 'bold');
        TT.setFigureAttributes('3d', {'x', 'y', 'z'});

    end




    
end 


  
% velocity measurements at specific point

kEta2 = (2:100).* TT.KOLMOGOROV_LENGTH;
EIS = 18/55 * 1.6 * kEta2.^(-5/3);





fi5=figure('Color','w','Position',[969 102 811 357]);

load('spec1024.mat')

k__ = (0:1023);
% [kEta_Pxx, E_Pxx] = scaleEnergySpectrum(k__(1:512), Pxx.(filenames(1))(2:513));
% loglog(kEta_Pxx, E_Pxx/170.1); hold on;



[kEta_pwr, E_pwr] = scaleEnergySpectrum(k__(1:512), pwr.('fieldname1')(1:512));
loglog(kEta_pwr, E_pwr/1024, 'color', colours1, 'LineWidth', 1.4); hold on;

% [kEta_Pxx, E_Pxx] = scaleEnergySpectrum(k__(1:512), Pxx.('fieldname1')(1:512));
% loglog(kEta_Pxx, E_Pxx/170.1, 'color', colours2, 'LineWidth', 1.4); hold on;
% 
[kEta_pwr, E_pwr] = scaleEnergySpectrum(k__(1:512), pwr.('fieldname2')(1:512));
loglog(kEta_pwr, E_pwr/1024, 'color', colours2, 'LineWidth', 1.4); hold on;

% [kEta_Pxx, E_Pxx] = scaleEnergySpectrum(k__(1:512), Pxx.('fieldname2')(1:512));
% loglog(kEta_Pxx, E_Pxx/170.1, 'color', colours4, 'LineWidth', 1.4); hold on;

[kEta_pwr, E_pwr] = scaleEnergySpectrum(k__(1:512), pwr.('fieldname3')(1:512));
loglog(kEta_pwr, E_pwr/1024, 'color', colours3, 'LineWidth', 1.4); hold on;

% % % [kEta_Pxx, E_Pxx] = scaleEnergySpectrum(k__(1:512), Pxx.(filenames(4))(2:513));
% % % % loglog(kEta_Pxx, E_Pxx/170.1); hold on;
% % 
% [kEta_pwr, E_pwr] = scaleEnergySpectrum(k__(1:512), pwr.('fieldname4')(1:512));
% loglog(kEta_pwr, E_pwr/1024, 'color', colours4, 'LineWidth', 1.4); hold on;
% % % 
% % % % % [kEta_Pxx, E_Pxx] = scaleEnergySpectrum(k__(1:512), Pxx.(filenames(4))(2:513));
% % % % % % loglog(kEta_Pxx, E_Pxx/170.1); hold on;
% % % % 
% [kEta_pwr, E_pwr] = scaleEnergySpectrum(k__(1:512), pwr.('fieldname5')(1:512));
% loglog(kEta_pwr, E_pwr/1024, 'color', colours5, 'LineWidth', 1.4); hold on;


%loglog(kEta_1024, E_1024/1024); hold on;
plot(kEta2, 0.1*EIS, 'color', '#000000', 'LineWidth', 1.6);


% Style figure
grid;
%text(kEta2(20)*2, EIS(20)*2, 'E_{11} = (18/55) * 1.6 * (k\eta)^{-5/3}');
title('Power spectrum', 'FontSize', 12, 'FontWeight', 'bold');
%legend('true field pwelch','true field fft power','predicted field pwelch 5 slices','predicted field fft power 5 slices','predicted field pwelch 7 slices','predicted field fft power 7 slices','-5/3')
%%legend('DNS data','predicted field ML','predicted field linear interpolation','predicted field 400 snapshots','predicted field linear interpolation with ML 2','-5/3','location','bestoutside')
%legend(leg(1),leg(2),'-5/3','location','bestoutside')
legend(leg(1),leg(2),leg(3),leg(4),leg(5),'-5/3','location','bestoutside')
ylabel('E_{11}/(\epsilon\nu^5)^{1/4}', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('k\eta', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');

% TT.makeFigureSquare(x_figure);



%%%%%%%%%%%%%%%%
% This function scales the provided power spectrum with the
% Kolmogorov length scale, dissipation rate and viscosity
function [kEta, E] = scaleEnergySpectrum(k, pwr)
    KOLMOGOROV_LENGTH = 0.00287;
    VISCOSITY = 0.000185;
    DISSIPATION_RATE = 0.0928;

    kEta = k.* KOLMOGOROV_LENGTH;
    E = pwr./(DISSIPATION_RATE * VISCOSITY^5)^(1/4);

end
