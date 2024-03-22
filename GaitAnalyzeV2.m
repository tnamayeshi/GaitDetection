function GaitAnalyzeV2(trialDir, mass, plotOptions, isRightLeg)

if(nargin < 1)
    trialDir = uigetdir();
    if isequal(trialDir,0)
       error('User selected Cancel');
    end
end

if nargin < 3
   plotOptions = 14; % 14 for vel 
end

plotExtra = bitand(plotOptions, 1) > 0;
showPlot  = bitand(plotOptions, 2) > 0;
savePlot  = ~showPlot || bitand(plotOptions, 4) > 0;
plotVel = bitand(plotOptions, 8) > 0;

if nargin < 4
   isRightLeg = contains(trialDir, '\RD\', 'IgnoreCase', true) || contains(path, '\RA\', 'IgnoreCase', true);
end

trcData = load(fullfile(trialDir, 'GaitDataTrcV2.csv.mat'));
grfData = load(fullfile(trialDir, 'GRFData.csv.mat'));

trc_fs = 1 / (trcData.Time(2) - trcData.Time(1));
grf_fs = 1 / (grfData.time(2) - grfData.time(1));

if grf_fs > trc_fs
    fs = trc_fs;
else
    fs = grf_fs;
end

if abs(fs - trc_fs) > 1
    [U, D] = rat(fs / trc_fs);
    trcData = resampleStruct(trcData, U, D);
end

if abs(fs - grf_fs) > 1
    [U, D] = rat(fs / grf_fs);
    grfData = resampleStruct(grfData, U, D);
end

clear trc_fs grf_fs;

FyL = grfData.ground_force_1_vy;
FyR = grfData.ground_force_2_vy;

if(nargin < 2)
	avgLen = floor(length(FyL) / 3);
    mass = (mean(FyL(10:avgLen)) + mean(FyR(10:avgLen))) / 9.81;
end

%%
if showPlot
    fig_h = figure;
else
    fig_h = figure('visible', 'off');
end

ftime_ms = grfData.time * 1000;
hold on;

if plotExtra
    plot(ftime_ms, FyR, 'LineWidth', 1.7, 'Color', [0.32,0.68,0.92], 'DisplayName', 'GRF R');
    plot(ftime_ms, FyL, 'LineWidth', 1.7, 'Color', [0.07,0.62,1.00], 'DisplayName', 'GRF L');
elseif isRightLeg
    plot(ftime_ms, FyR, 'LineWidth', 1.7, 'Color', [0.32,0.68,0.92], 'DisplayName', 'Vertical GRF');
else
    plot(ftime_ms, FyL, 'LineWidth', 1.7, 'Color', [0.32,0.68,0.92], 'DisplayName', 'Vertical GRF');
end

xlabel('time (ms)');
ylabel('Force (N)');

pos = trcData.RTOE_X / 1000;
vel_r = (pos(2:end) - pos(1:end-1)) * fs;
pos = trcData.LTOE_X / 1000;
vel_l = (pos(2:end) - pos(1:end-1)) * fs;


vtime_ms = trcData.Time(1:end-1) * 1000.0;
delta_sample =  floor(((ftime_ms(1)-vtime_ms(1)) / 1000) * fs);

if isRightLeg
    vel = vel_r;
else
    vel = vel_l;
end

if plotVel
	yyaxis right;
	ylabel('Velocity (m/s)');
	plot(vtime_ms, vel, 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
end

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
	
	if plotVel
		
		if i == 1
			plot(vt, subvel,'LineWidth',1.7, 'Marker', 'none', 'LineStyle','-.', 'Color',[0.96,0.06,0.60],'DisplayName','Leg Speed');
		else
			plot(vt, subvel,'LineWidth',1.7, 'Marker', 'none', 'LineStyle','-.', 'Color',[0.96,0.06,0.60],'HandleVisibility', 'off');
		end
		
	end
	
    % plot([vt(1),vt(end)], [avg avg], 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'LineStyle', '-', 'Marker', 'none');
    
    dv = abs(avg + 1);
    if dv > max_dv
        max_dv = dv;        
        max_dv_i = i;
        %max_v_t = vt;
        %max_v_v = subvel;
    end
end

if plotVel
    yyaxis left;
end

if max_dv_i > 0
    perturb_sample = hill_strikes_x(max_dv_i);
else
    perturb_sample = -1;
	warning('No perturb detected!');
end
%perturb_sample_end = toe_offs_x(max_dv_i);
%% Plot data

hill_strikes_r_ms = ftime_ms(1) + (hill_strikes_r * 1000.0 / fs);
toe_offs_r_ms = ftime_ms(1) + (toe_offs_r * 1000.0 / fs);
hill_strikes_l_ms = ftime_ms(1) + (hill_strikes_l * 1000.0 / fs);
toe_offs_l_ms = ftime_ms(1) + (toe_offs_l * 1000.0 / fs);
perturb_time_ms = ftime_ms(1) + (perturb_sample * 1000.0 / fs);
%perturb_sample_end_ms = ftime_ms(1) + (perturb_sample_end * 1000.0 / fs);

if plotExtra || isRightLeg
    plot(hill_strikes_r_ms, FyR(hill_strikes_r), 'ro','DisplayName', 'Hill Strike');
    plot(toe_offs_r_ms, FyR(toe_offs_r), 'go','DisplayName', 'Toe Off');
end

if plotExtra || ~isRightLeg
    plot(hill_strikes_l_ms, FyL(hill_strikes_l), 'ro','DisplayName', 'Hill Strike');
    plot(toe_offs_l_ms, FyL(toe_offs_l), 'go','DisplayName', 'Toe Off');
end

if perturb_sample > 0
    xline(perturb_time_ms, '--b', {''}, 'Color',[0.717647058823529 0.274509803921569 1],...
	    'LineStyle','-.',...
	    'LineWidth',1.5,...
	    'DisplayName','Perturbation',...
	    'Label',{''}, 'Alpha', 1);
end

plot(trcData.Time * 1000, trcData.RHEE_X, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Right Heel X');
plot(trcData.Time * 1000, trcData.LHEE_X, 'Color', [0.39, 0.83, 0.07], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Left Heel X');
plot(trcData.Time * 1000, abs(trcData.RHEE_X-trcData.LHEE_X), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '-.', 'DisplayName', 'Heels X Distance');
plot(trcData.Time * 1000, abs(trcData.RHEE_Z-trcData.LHEE_Z), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Heels Z Distance');

%xline(perturb_sample_end_ms, '--r');
%perturb_rng_t = [max_v_t(1)-0.1 max_v_t(end)+0.1];
%perturb_rng_v = [min(max_v_v)-0.1 max(max_v_v)+0.1];

%patch([perturb_rng_t(1) perturb_rng_t(1) perturb_rng_t(2) perturb_rng_t(2)], [perturb_rng_v(1) perturb_rng_v(2) perturb_rng_v(2) perturb_rng_v(1)], [0.8 0.8 0.8]);
%alpha(0.1)

legend();
%% Save Data



data = {};
data.version = 4;
data.mass = mass;
data.perturb_leg = 'L';

if isRightLeg
	data.perturb_leg = 'R';
end

data.hill_strikes_r_ms = hill_strikes_r_ms;
data.hill_strikes_l_ms = hill_strikes_l_ms;
data.toe_offs_r_ms = toe_offs_r_ms;
data.toe_offs_l_ms = toe_offs_l_ms;
data.perturb_time_ms = perturb_time_ms;

save(fullfile(trialDir, 'neededDataV4.mat'), 'data');

if plotVel
	yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
	yyaxis left; yliml = get(gca,'Ylim');

	if yliml(2)*ratio<yliml(1)
		set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
	else
		set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
	end
	
	ax = gca;
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'k';
end 

tt = trialDir(19:end);
tt = strrep(tt, '\', '-');
title(tt, 'Interpreter', 'none');


if savePlot
    print(fig_h, fullfile(trialDir, 'Perturbation.png'), '-dpng', '-r300');
    set(fig_h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    savefig(fig_h, fullfile(trialDir, 'Perturbation.fig'), 'compact');
end

if ~showPlot
   close(fig_h);
end


function sdatar = resampleStruct(sdata, U, D)
	fns = fieldnames(sdata);
	for i = 1:length(fns)
		d = sdata.(string(fns(i)));
		d = resample(d, U, D);
		d = d((U+D):end-(U+D));
		sdata.(string(fns(i))) = d;
	end
	sdatar = sdata;
end
end