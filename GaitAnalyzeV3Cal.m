function GaitAnalyzeV3Cal(calDir, plotOptions)

if(nargin < 1)
    calDir = uigetdir();
    if isequal(calDir,0)
       error('User selected Cancel');
    end
end

if nargin < 2
   plotOptions = 15; % 14 for vel 
end

rawData = load(fullfile(calDir, 'caldata.mat'));

plotExtra = bitand(plotOptions, 1) > 0;
showPlot  = bitand(plotOptions, 2) > 0;
savePlot  = ~showPlot || bitand(plotOptions, 4) > 0;
plotVel = bitand(plotOptions, 8) > 0;

mass = rawData.mass;
rawData.weight = mass * 10;
rawData.dataType = 'gait';

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
    plot(ftime_ms, FyL, 'LineWidth', 1.7, 'Color', [0.77,0.62,0.99], 'DisplayName', 'GRF L');
end

xlabel('time (ms)');
ylabel('Force (N)');

pos = rawData.trc.RTOE_X / 1000;
vel_r = (pos(2:end) - pos(1:end-1)) * fs;
pos = rawData.trc.LTOE_X / 1000;
vel_l = (pos(2:end) - pos(1:end-1)) * fs;


vtime_ms = rawData.trc.time(1:end-1) * 1000.0;
delta_sample =  floor(((ftime_ms(1)-vtime_ms(1)) / 1000) * fs);


if plotVel
	yyaxis right;
	ylabel('Velocity (m/s)');
	plot(vtime_ms, vel_r, 'LineWidth', 1, 'LineStyle', '--', 'Color', [0.75 0.85 0.75], 'HandleVisibility', 'off');
    plot(vtime_ms, vel_l, 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.25 0.75 0.75], 'HandleVisibility', 'off');
    yyaxis left;
end



[hill_strikes_r, toe_offs_r, hs_vr, to_vr] = extractLegGaitEvents(FyR, vel_r, delta_sample, mass, calDir, 'R Leg');
[hill_strikes_l, toe_offs_l, hs_vl, to_vl] = extractLegGaitEvents(FyL, vel_l, delta_sample, mass, calDir, 'L Leg');


hill_strikes_r_ms = ftime_ms(1) + (hill_strikes_r * 1000.0 / fs);
toe_offs_r_ms = ftime_ms(1) + (toe_offs_r * 1000.0 / fs);
hill_strikes_l_ms = ftime_ms(1) + (hill_strikes_l * 1000.0 / fs);
toe_offs_l_ms = ftime_ms(1) + (toe_offs_l * 1000.0 / fs);

%perturb_sample_end_ms = ftime_ms(1) + (perturb_sample_end * 1000.0 / fs);

if plotExtra
    plot(hill_strikes_r_ms, FyR(hill_strikes_r), 'ro','DisplayName', 'Hill Strike R');
    plot(toe_offs_r_ms, FyR(toe_offs_r), 'go','DisplayName', 'Toe Off R');
end

if plotExtra
    plot(hill_strikes_l_ms, FyL(hill_strikes_l), 'bo','DisplayName', 'Hill Strike L');
    plot(toe_offs_l_ms, FyL(toe_offs_l), 'co','DisplayName', 'Toe Off L');
end


plot(rawData.trc.time * 1000, rawData.trc.RHEE_X, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Right Heel X');
plot(rawData.trc.time * 1000, rawData.trc.LHEE_X, 'Color', [0.39, 0.83, 0.07], 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'Left Heel X');
plot(rawData.trc.time * 1000, abs(rawData.trc.RHEE_X-rawData.trc.LHEE_X), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '-.', 'DisplayName', 'Heels X Distance');
plot(rawData.trc.time * 1000, abs(rawData.trc.RHEE_Z-rawData.trc.LHEE_Z), 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Heels Z Distance');


legend();
%% Save Data

steps_r_i = zeros(1, length(hill_strikes_r_ms));
steps_l_i = zeros(1, length(hill_strikes_l_ms));

steps_r = zeros(2, length(hill_strikes_r_ms));



for ix = 1:length(hill_strikes_r_ms)
    steps_r(1, ix) = hill_strikes_r_ms(ix);
    steps_r(2, ix) = toe_offs_r_ms(ix);
    steps_r_i(ix) = ix;
end

steps_l = zeros(2, length(hill_strikes_l_ms));
for ix = 1:length(hill_strikes_l_ms)
    steps_l(1, ix) = hill_strikes_l_ms(ix);
    steps_l(2, ix) = toe_offs_l_ms(ix);
    steps_l_i(ix) = ix;
end


rawData.hillStrikesR_ms = hill_strikes_r_ms(:);
rawData.hillStrikesL_ms = hill_strikes_l_ms(:);
rawData.toeOffsR_ms = toe_offs_r_ms(:);
rawData.toeOffsL_ms = toe_offs_l_ms(:);
rawData.stepsR_ms = steps_r;
rawData.stepsL_ms = steps_l;
rawData.stepsRelativeIndexesR = steps_r_i(:);
rawData.stepsRelativeIndexesL = steps_l_i(:);

save(fullfile(calDir, 'gaitdata.mat'), '-struct', 'rawData');

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

tt = calDir(19:end);
tt = strrep(tt, '\', '-');
title(tt, 'Interpreter', 'none');


if savePlot
    fig_h.PaperPositionMode = 'auto';
    axis tight;
    set(fig_h,'units','normalized','outerposition',[0 0 1 1])
    print(fig_h, fullfile(calDir, 'Perturbation.png'), '-dpng', '-r300');
    set(fig_h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
    savefig(fig_h, fullfile(calDir, 'Perturbation.fig'), 'compact');
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