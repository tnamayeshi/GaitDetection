function batchProcessAvg(pd, dataDirs, isRightLeg)

if pd(end) ~= '/' && pd(end) ~= '\'
    pd = strcat(pd, '\');
end
showFigure = true;
doSubplot = true;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir([pd 'Average']);

%dataDirs = {'1', '2', '3'};
fields = {'ground_force_1_vy', 'ground_force_2_vy'};
steps = [-1 0 1];

%% export data from trials

fieldAndOptions = {};
fieldAndOptions.dosum = 0;
fieldAndOptions.scale = 1;
if isRightLeg
    fieldAndOptions.leftLeg = 0;
else
    fieldAndOptions.leftLeg = 1;
end
fieldAndOptions.showAlt = 0;
fieldAndOptions.export = 1;
fieldAndOptions.ho = 0;
fieldAndOptions.continus = 1;
fieldAndOptions.prcnt = 1;
fieldAndOptions.sep = 0;
fieldAndOptions.lowPassFreq = 0;
fieldAndOptions.steps = steps;
fieldAndOptions.selected_fields = fields;

for i = 1:length(dataDirs)
    dd = dataDirs{i};
    PlotFile(load([pd dd '/neededDataV3.mat']), [pd dd '/MOT.mat'], fieldAndOptions, [pd dd '/Exports/'], 1);
end

if fieldAndOptions.dosum > 0
    sum_field = fields{1};
    for i = 2:length(fields)
        sum_field = [sum_field '+' fields{i}];
    end
    fields = {sum_field};
end


%%
for fi = 1:length(fields)
    if showFigure
        fig_h = figure('name', fields{fi});
    else
        fig_h = figure('name', fields{fi}, 'visible', 'off');
    end
    if doSubplot
        subplot(2,2,1);
    end
    title(fields{fi}, 'Interpreter', 'none');
    legend();
    hold on;
    for si = 1:length(steps)
        if doSubplot
            subplot(2,2,1);
            hold on;
        else
            set(0, 'CurrentFigure', fig_h)
        end
        data = cell(1, length(dataDirs));
        data_prv = cell(1, length(dataDirs));
        data_nxt = cell(1, length(dataDirs));
        data_t = cell(1, length(dataDirs));
        max_samples = 0;
        for di = 1:length(dataDirs)
            dd = dataDirs{di};
            datax_file = [pd dd '/Exports/' fields{fi} '/' num2str(steps(si)) '.mat'];
            d = load(datax_file);
            data{di} = d.data;
            data_t{di} = d.data_tm;
            data_prv{di} = d.prv_data;
            data_nxt{di} = d.nxt_data;
            if length(d.data) >= max_samples
                max_samples = length(d.data);
            end
        end
        
        data_mat = zeros(length(dataDirs), max_samples);
        data_t_mat = zeros(length(dataDirs), max_samples);
        
        for di = 1:length(dataDirs)
            if length(data{di}) ~= max_samples
                [U, D] = rat((max_samples)/length(data{di}));
                tmp = [data_prv{di} data{di} data_nxt{di}];
                tmp = resample(tmp, U, D);
                cnt_extra = length(tmp) - max_samples;
                cnt_prv = round(cnt_extra * length(data_prv{di}) * 1.0 / (length(data_prv{di})+length(data_nxt{di})));
                cnt_nxt = cnt_extra - cnt_prv;
                data_mat(di, :) = tmp((1 + cnt_prv):(end - cnt_nxt));
                t = data_t{di};
                data_t_mat(di, :) = linspace(t(1), t(end), max_samples);
            else
                data_mat(di, :) = data{di};
                data_t_mat(di, :) = data_t{di};
            end
        end
        
        data_max = max(data_mat);
        data_min = min(data_mat);
        data_avg = sum(data_mat) / length(dataDirs);
        
        
        
        %figure('name', fields{fi});
        %hold on;
        
        %        plot(data_t_mat(1,:), data_max, 'LineStyle', '--', 'Color', 0.4 * [1 1 1], 'HandleVisibility', 'off');
        %        plot(data_t_mat(1,:), data_min, 'LineStyle', '--', 'Color', 0.4 * [1 1 1], 'HandleVisibility', 'off');
        %        title(['Step ' num2str(steps(si))]);
        %
        x = [data_t_mat(1,:), fliplr(data_t_mat(1,:)) data_t_mat(1,1)];
        inBetween = [data_max, fliplr(data_min) data_max(1)];
        % fill(x, inBetween, 0.95 * [1 1 1], 'LineStyle', 'none');
        
        plot(data_t_mat(1,:), data_avg, 'DisplayName', ['Step ' num2str(steps(si))], 'LineStyle', '-');
        
        fname = [fields{fi} ' Step ' num2str(steps(si))];
        save([pd 'Average\' fname], 'data_avg');
        
        if doSubplot
            subplot(2,2,1 + si);
        else
            if showFigure
                fig2_h = figure('name', fields{fi});
            else
                fig2_h = figure('name', fields{fi}, 'visible', 'off');
            end
        end
        
        title([fields{fi} ' Step ' num2str(steps(si))], 'Interpreter', 'none');
        legend('Location', 'southeast');
        hold on;
        plot(data_t_mat(1,:), data_avg, 'DisplayName', 'Average', 'LineStyle', '-', 'LineWidth', 2);
        for i = 1:length(dataDirs)
            sd = round(std(data_mat(i,:)-data_avg, 0, 2), 1);
            plot(data_t_mat(1,:), data_mat(i,:), 'DisplayName', ['Trial #' num2str(i) ' SD=' num2str(sd)], 'LineStyle', '-', 'LineWidth', 1);
        end
        
        if ~doSubplot
            fig2_h.PaperPositionMode = 'auto';
            axis tight;
            print(fig2_h, fullfile([pd 'Average\'], [fields{fi} ' S' num2str(steps(si)) '.png']), '-dpng');
            print(fig2_h, fullfile([pd 'Average\'], [fields{fi} ' S' num2str(steps(si)) '.pdf']), '-dpdf');
            set(fig2_h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            savefig(fig2_h, fullfile([pd 'Average\'], [fields{fi} ' S' num2str(steps(si)) '.fig']), 'compact');
        end
        
    end
    fig_h.PaperPositionMode = 'auto';
    axis tight;
    print(fig_h, fullfile([pd 'Average\'], [fields{fi} '.png']), '-dpng');
    print(fig_h, fullfile([pd 'Average\'], [fields{fi} '.pdf']), '-dpdf');
    set(fig_h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(fig_h, fullfile([pd 'Average\'], [fields{fi} '.fig']), 'compact');
end

