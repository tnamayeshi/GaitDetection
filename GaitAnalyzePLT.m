function GaitAnalyzePLT(reducedTrcFile, grfFile, mass, plotOptions, isRightLeg, subplt)

if(nargin < 1)
    [file, path] = uigetfile('GaitDataTrc.csv');
    if isequal(file,0)
       disp('User selected Cancel');
       return;
    else
       reducedTrcFile = fullfile(path,file);
    end
end

if(nargin < 2)
    [file, path] = uigetfile('GRFData.csv');
    if isequal(file,0)
       disp('User selected Cancel');
       return;
    else
       grfFile = fullfile(path,file);
    end
end



if nargin < 4
   plotOptions = 6; 
end

if nargin < 6
   subplt = 0; 
end

plotExtra = bitand(plotOptions, 1) > 0;
showPlot  = bitand(plotOptions, 2) > 0;
savePlot  = ~showPlot || bitand(plotOptions, 4) > 0;

fid = fopen(reducedTrcFile);
header = fgetl(fid);
if nargin < 5
   isRightLeg = contains(path, '\RD\', 'IgnoreCase', true) || contains(path, '\RA\', 'IgnoreCase', true);
end

fclose(fid);

trcData = readmatrix(reducedTrcFile, 'HeaderLines', 1);
grfData = readmatrix(grfFile, 'HeaderLines', 1);

trc_fs = 1 / (trcData(2, 1) - trcData(1, 1));
grf_fs = 1 / (grfData(2, 1) - grfData(1, 1));

if grf_fs > trc_fs
    fs = trc_fs;
else
    fs = grf_fs;
end

if abs(fs - trc_fs) > 1
    [U, D] = rat(fs / trc_fs);
    trcData = resample(trcData, U, D);
    trcData = trcData((U+D):end-(U+D), :);
end

if abs(fs - grf_fs) > 1
    [U, D] = rat(fs / grf_fs);
    grfData = resample(grfData, U, D);
    grfData = grfData((U+D):end-(U+D), :);
end

clear trc_fs grf_fs;

FyL = grfData(:, 2);
FyR = grfData(:, 3);

if(nargin < 3)
	avgLen = floor(length(FyL) / 3);
    mass = (mean(FyL(10:avgLen)) + mean(FyR(10:avgLen))) / 9.81;
end

%%
if showPlot
    if subplt < 2
        fig_h = figure;
    end
    
    if subplt > 0
        savePlot = 0;
        subplot(2, 2, subplt)
    end
else
    fig_h = figure('visible', 'off');
end

ftime_ms = grfData(:, 1) * 1000;
hold on;

if plotExtra
    plot(ftime_ms, FyR/mass, 'LineWidth', 1.7, 'Color', [0.07,0.62,1.00], 'DisplayName', 'GRF R');
    plot(ftime_ms, FyL/mass, 'LineWidth', 1.7, 'Color', [0.07,0.62,1.00], 'DisplayName', 'GRF L');
elseif isRightLeg
    plot(ftime_ms, FyR/mass, 'LineWidth', 1.7, 'Color', [0.07,0.62,1.00], 'DisplayName', 'Vertical GRF');
else
    plot(ftime_ms, FyL/mass, 'LineWidth', 1.7, 'Color', [0.07,0.62,1.00], 'DisplayName', 'Vertical GRF');
end

xlabel('time (ms)');
ylabel('Force (N/kg)');

pos = trcData(:, 3) / 1000;
vel_r = (pos(2:end) - pos(1:end-1)) * fs;
pos = trcData(:, 5) / 1000;
vel_l = (pos(2:end) - pos(1:end-1)) * fs;


vtime_ms = trcData(1:end-1, 1) * 1000.0;
delta_sample =  floor(((ftime_ms(1)-vtime_ms(1)) / 1000) * fs);

if isRightLeg
    vel = vel_r;
else
    vel = vel_l;
end
yyaxis right;
ylabel('Velocity (m/s)');
plot(vtime_ms, vel, 'LineWidth', 1.4, 'LineStyle', '-.', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');

[hill_strikes_r, toe_offs_r, hs_vr, to_vr] = extractLegGaitEvents(FyR, vel_r, delta_sample, mass);
[hill_strikes_l, toe_offs_l, hs_vl, to_vl] = extractLegGaitEvents(FyL, vel_l, delta_sample, mass);


%% Perturbation sample find



if isRightLeg
    hill_strikes_x = hill_strikes_r;
    %toe_offs_x = toe_offs_r;
    vel = vel_r;
    hs_v = hs_vr;
    to_v = to_vr;
else
    hill_strikes_x = hill_strikes_l;
    %toe_offs_x = toe_offs_l;
    vel = vel_l;
    hs_v = hs_vl;
    to_v = to_vl;
end

max_dv = -1;
max_dv_i = -1;
%max_v_t = [];
%max_v_v = [];




for i = 1:length(hill_strikes_x)
    
    vti = hs_v(i):to_v(i);
    subvel = vel(vti);
    vt = vtime_ms(vti);
    if isempty(vt)
        continue;
    end
    avg = mean(subvel);
    if i == 1
        plot(vt, subvel,'LineWidth',2, 'LineStyle',':', 'Color',[1.00,0.41,0.16],'DisplayName','Leg Speed');
    else
        plot(vt, subvel,'LineWidth',2, 'LineStyle',':', 'Color',[1.00,0.41,0.16],'HandleVisibility', 'off');
    end
    
    
    dv = abs(avg + 1);
    if dv > max_dv
        max_dv = dv;        
        max_dv_i = i;
        %max_v_t = vt;
        %max_v_v = subvel;
    end
end

yyaxis left;
perturb_sample = hill_strikes_x(max_dv_i);
%perturb_sample_end = toe_offs_x(max_dv_i);
%% Plot data

hill_strikes_r_ms = ftime_ms(1) + (hill_strikes_r * 1000.0 / fs);
toe_offs_r_ms = ftime_ms(1) + (toe_offs_r * 1000.0 / fs);
hill_strikes_l_ms = ftime_ms(1) + (hill_strikes_l * 1000.0 / fs);
toe_offs_l_ms = ftime_ms(1) + (toe_offs_l * 1000.0 / fs);
perturb_time_ms = ftime_ms(1) + (perturb_sample * 1000.0 / fs);
%perturb_sample_end_ms = ftime_ms(1) + (perturb_sample_end * 1000.0 / fs);

if plotExtra
    plot(trcData(:, 1) * 1000.0, trcData(:, 3));
end

if plotExtra || isRightLeg
    %plot(hill_strikes_r_ms, FyR(hill_strikes_r)/mass, 'o');
    %plot(toe_offs_r_ms, FyR(toe_offs_r)/mass, 'o');
end

if plotExtra || ~isRightLeg
    %plot(hill_strikes_l_ms, FyL(hill_strikes_l)/mass, 'o');
    %plot(toe_offs_l_ms, FyL(toe_offs_l)/mass, 'o');
end

xline(perturb_time_ms, '--b', {''}, 'Color',[0.717647058823529 0.274509803921569 1],...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'DisplayName','Perturbation',...
    'Label',{''}, 'Alpha', 1);
%xline(perturb_sample_end_ms, '--r');
%perturb_rng_t = [max_v_t(1)-0.1 max_v_t(end)+0.1];
%perturb_rng_v = [min(max_v_v)-0.1 max(max_v_v)+0.1];

%patch([perturb_rng_t(1) perturb_rng_t(1) perturb_rng_t(2) perturb_rng_t(2)], [perturb_rng_v(1) perturb_rng_v(2) perturb_rng_v(2) perturb_rng_v(1)], [0.8 0.8 0.8]);
%alpha(0.1)
if subplt == 1
legend();
end
%% Save Data

data = {};
data.version = 3;
data.mass = mass;
data.hill_strikes_r_ms = hill_strikes_r_ms;
data.hill_strikes_l_ms = hill_strikes_l_ms;
data.toe_offs_r_ms = toe_offs_r_ms;
data.toe_offs_l_ms = toe_offs_l_ms;
data.perturb_time_ms = perturb_time_ms;
[filepath, ~, ~] = fileparts(reducedTrcFile);

save(fullfile(filepath, 'neededDataV3.mat'), 'data');

yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
yyaxis left; yliml = get(gca,'Ylim');
if yliml(2)*ratio<yliml(1)
    set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
else
    set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
end

tt = filepath(19:end);
tt = strrep(tt, '\', '-');
title(tt, 'Interpreter', 'none');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

if savePlot
    print(fig_h, fullfile(filepath, 'Perturbation.png'), '-dpng', '-r300');
    set(fig_h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    savefig(fig_h, fullfile('./figs/', [tt '.fig']), 'compact');
end

if ~showPlot
   close(fig_h);
end

function [hs, to, hs_v, to_v] = extractLegGaitEvents(grf, velo, delta_s, mass)

	thrshld = max(mass / 8, 10);
    vel = velo * 100;
    dv_th = 0.5;
	
	a = grf > thrshld;
	b = a(2:end) - a(1:end-1);

	hs = find(b > 0);
	to = find(b < 0) + 1;


	if length(hs) < 2 || length(to) < 2
		hs = [];
		to = [];
	   return; 
	end
	
	
	% Remove unpaired hill strikes or toe offs

	if hs(1) > to(1)
		to = to(2:end);
	end

	if hs(end) > to(end)
		hs = hs(1:end-1);
    end
        
	i = 1;
	while i <= length(hs)
        
        e = hs(i);
        if i == 1
            s = 1;
        else
            s = to(i - 1);
        end
        
        grf_dv = grf((s + 1):e) - grf(s:(e - 1));
        x = find(grf_dv <= dv_th, 1, 'last');
        
        if isempty(x)
           hs(i) = []; 
           to(i) = [];
           continue;
        end
        
        e = e - (length(grf_dv) - x);
        hs(i) = e - 1;
        
        s = to(i);
        if i == length(hs)
            e = length(grf);
        else
            e = hs(i + 1);
        end
        
        grf_dv = grf((s + 1):e) - grf(s:(e - 1));
        x = find(grf_dv >= -dv_th, 1);
        
        if isempty(x)
           hs(i) = []; 
           to(i) = [];
           continue;
        end
        
        s = s + x - 2;
        to(i) = s;
		i = i + 1;
    end
    
    
    hs = hs(:)';
    to = to(:)';
    
    
    % remove short/invalid steps
    
    i = 1;
	while i <= length(hs)
       dt = hs(i):to(i);
       if length(dt) < 10
           hs(i) = [];
           to(i) = [];
           continue;
       end
       
       stp = grf(dt);
       
       if (max(stp) - min(stp)) < mass / 8
           hs(i) = [];
           to(i) = [];
           continue;
       end
       
       i = i + 1;
    end
    
    
    i = 1;
    hs_v = zeros(1, length(hs));
    to_v = zeros(1, length(hs));
    
	while i <= length(hs)
       vt_i_e = delta_s + to(i);
       if vt_i_e > length(vel)
           vt_i_e = length(vel);
       end
       vt_i = (delta_s + hs(i)):vt_i_e;
       
       svel = vel(vt_i);
       ds = find(svel < -1, 1);
       de = length(svel) - find(svel < -1, 1, 'last');
       
       if isempty(ds)
           ds = 1;
           de = 1;
       end
       
       hs_v(i) = delta_s + (hs(i) + ds);
       to_v(i) = delta_s + (to(i) - de);
       
       if (ds + de) == 2
           i = i + 1;
           continue;
       end
       
       if to_v(i) > length(vel)
           to_v(i) = length(vel);
       end
       
       vt_i2 = hs_v(i):to_v(i);
       svel2 = vel(vt_i2);
       x = find(svel2 > 10, 1);
       
       if ~isempty(x)
           x = x + ds;
           de2 = find(svel(1:x) < -1, 1, 'last');
           ds2 = x + find(svel(x:end) < -1, 1);
           to = [to(1:i-1) hs(i)+de2 to(i:end)];
           hs = [hs(1:i) hs(i)+ds2 hs(i+1:end)];
           svel = vel(hs_v(i)+(1:x));
           de2 = find(svel < -1, 1, 'last');
           to_v(i) = hs_v(i) + de2;
       end
       
       i = i + 1;
    end
    
end
end