function GaitAnalyzeV3(trialDir, plotOptions)

if(nargin < 1)
    trialDir = uigetdir();
    if isequal(trialDir,0)
       error('User selected Cancel');
    end
end

if nargin < 2
   plotOptions = 14; % 14 for vel 
end

rawData = load(fullfile(trialDir, 'rawdata.mat'));

plotExtra = bitand(plotOptions, 1) > 0;
showPlot  = bitand(plotOptions, 2) > 0;
savePlot  = ~showPlot || bitand(plotOptions, 4) > 0;
plotVel = bitand(plotOptions, 8) > 0;

isRightLeg = rawData.perturbedLeg == 'R';
pertOnsetIsHs = rawData.perturbationOnset == 'HS';
mass = rawData.weight;
rawData.weight = mass * 10;
rawData.dataType = 'gait';
rawData.mass = mass;

trc_fs = 1 / (rawData.trc.time(2) - rawData.trc.time(1));
mot_fs = 1 / (rawData.mot.time(2) - rawData.mot.time(1));

fs = 100;

if abs(fs - trc_fs) > 1
    [U, D] = rat(fs / trc_fs);
    rawData.trc = resampleStruct(rawData.trc, U, D);
end

if abs(fs - mot_fs) > 1
    [U, D] = rat(fs / mot_fs);
    rawData.mot = resampleStruct(rawData.mot, U, D);
end

clear trc_fs mot_fs;

FyL = rawData.mot.ground_force_1_vy;
FyR = rawData.mot.ground_force_2_vy;


%%
if showPlot
    fig_h = figure;
else
    fig_h = figure('visible', 'off');
end

ftime_ms = rawData.mot.time * 1000;
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

pos = rawData.trc.RTOE_X / 1000;
vel_r = (pos(2:end) - pos(1:end-1)) * fs;
pos = rawData.trc.LTOE_X / 1000;
vel_l = (pos(2:end) - pos(1:end-1)) * fs;


vtime_ms = rawData.trc.time(1:end-1) * 1000.0;
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

[hill_strikes_r, toe_offs_r, hs_vr, to_vr] = extractLegGaitEvents(FyR, vel_r, delta_sample, mass, trialDir, 'R Leg');
[hill_strikes_l, toe_offs_l, hs_vl, to_vl] = extractLegGaitEvents(FyL, vel_l, delta_sample, mass, trialDir, 'L Leg');


%% Perturbation sample find



if isRightLeg
    hill_strikes_x = hill_strikes_r;
    toe_offs_x = toe_offs_r;
    vel = vel_r;
    hs_v = hs_vr;
    to_v = to_vr;
else
    hill_strikes_x = hill_strikes_l;
    toe_offs_x = toe_offs_l;
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
    if pertOnsetIsHs
        perturb_sample = hill_strikes_x(max_dv_i);
    else
        perturb_sample = hill_strikes_x(max_dv_i);
    end
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

plot(rawData.trc.time * 1000, rawData.trc.RHEE_X, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Right Heel X');
plot(rawData.trc.time * 1000, rawData.trc.LHEE_X, 'Color', [0.39, 0.83, 0.07], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Left Heel X');
plot(rawData.trc.time * 1000, abs(rawData.trc.RHEE_X-rawData.trc.LHEE_X), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '-.', 'DisplayName', 'Heels X Distance');
plot(rawData.trc.time * 1000, abs(rawData.trc.RHEE_Z-rawData.trc.LHEE_Z), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Heels Z Distance');

%xline(perturb_sample_end_ms, '--r');
%perturb_rng_t = [max_v_t(1)-0.1 max_v_t(end)+0.1];
%perturb_rng_v = [min(max_v_v)-0.1 max(max_v_v)+0.1];

%patch([perturb_rng_t(1) perturb_rng_t(1) perturb_rng_t(2) perturb_rng_t(2)], [perturb_rng_v(1) perturb_rng_v(2) perturb_rng_v(2) perturb_rng_v(1)], [0.8 0.8 0.8]);
%alpha(0.1)

legend();
%% Save Data

perturb_s_index_r = -1;
perturb_s_index_l = -1;

steps_r_i = [];
steps_l_i = [];

if perturb_time_ms > 0
    
   perturb_s_index_r = find(toe_offs_r_ms > perturb_time_ms, 1);
   perturb_s_index_l = find(toe_offs_l_ms > perturb_time_ms, 1);
   
   steps_r_i = zeros(1, length(hill_strikes_r_ms));
   steps_l_i = zeros(1, length(hill_strikes_l_ms));
end

steps_r = zeros(2, length(hill_strikes_r_ms));



for ix = 1:length(hill_strikes_r_ms)
    steps_r(1, ix) = hill_strikes_r_ms(ix);
    steps_r(2, ix) = toe_offs_r_ms(ix);
    if perturb_time_ms > 0
        steps_r_i(ix) = ix - perturb_s_index_r;
    end
end

steps_l = zeros(2, length(hill_strikes_l_ms));
for ix = 1:length(hill_strikes_l_ms)
    steps_l(1, ix) = hill_strikes_l_ms(ix);
    steps_l(2, ix) = toe_offs_l_ms(ix);
    if perturb_time_ms > 0
        steps_l_i(ix) = ix - perturb_s_index_l;
    end
end


rawData.hillStrikesR_ms = hill_strikes_r_ms(:);
rawData.hillStrikesL_ms = hill_strikes_l_ms(:);
rawData.toeOffsR_ms = toe_offs_r_ms(:);
rawData.toeOffsL_ms = toe_offs_l_ms(:);
rawData.perturbationTime_ms = perturb_time_ms;
rawData.perturbationStepIndexR = perturb_s_index_r;
rawData.perturbationStepIndexL = perturb_s_index_l;
rawData.stepsR_ms = steps_r;
rawData.stepsL_ms = steps_l;
rawData.stepsRelativeIndexesR = steps_r_i(:);
rawData.stepsRelativeIndexesL = steps_l_i(:);

save(fullfile(trialDir, 'gaitdata.mat'), '-struct', 'rawData');

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
    fig_h.PaperPositionMode = 'auto';
    axis tight;
    set(fig_h,'units','normalized','outerposition',[0 0 1 1])
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