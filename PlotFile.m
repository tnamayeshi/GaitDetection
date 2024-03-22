function PlotFile(neededData, filename, fieldAndOptions, exportDir, noPlot)

pf_version = 0.93;
data_version = 3;

%%
disp(['Plot File version ' num2str(pf_version)])

path = '';

if nargin < 2
    [file, path] = uigetfile({'*.mat;*.sto;*.mot'});
    if file == 0
        error('Operation cancel by user!');
    end
    filename = fullfile(path, file);
else
    file = filename;
end


if nargin < 1
	nd_filename = fullfile(path, 'neededDataV3.mat');
	if ~isfile(nd_filename)
		nd_filename = fullfile(path, '../neededDataV3.mat');
	end
    neededData = load(nd_filename);
end

if nargin < 4
    exportDir = './Exports/';
end

if nargin < 5
    noPlot = 0;
end

if isfield(neededData, 'data')
    neededData = neededData.data;
end

if ~isfield(neededData, 'version')
    error('Invalid Needed Data File !');
elseif neededData.version ~= data_version
    error(['Invalid Needed Data Version ! Expected ' num2str(data_version) ' but got ' num2str(neededData.version)]);
end


if isequal(file, 0)
    warning('User selected Cancel');
    return;
else
    [data, fields] = load_data_file(filename);
    time_inx = -1;
    if strcmpi(fields{1}, 'time')
        time_inx = 1;
    elseif strcmpi(fields{2}, 'time')
        time_inx = 2;
    else
        error('Invalid file selected!');
    end
    
    datac = {};
    i_cnt = 0;
    [z, data_size] = size(data);
    
    if z > data_size
        data = data';
        data_size = z;
    end
    
    for i=time_inx:data_size
        i_cnt = i_cnt + 1;
        datac{i_cnt} = data{:, i};
    end
    data = datac;
    data_time_ms = data{1} * 1000.0;
    diff_time = data_time_ms(2:end) - data_time_ms(1:end-1);
    diff_time = round(diff_time(2:end) - diff_time(1:end-1), 5);
    if ~isempty(find(diff_time, 1))
        error('Variable sample rate not supported!');
    end
    data_fs = 1000.0 / (data_time_ms(2)-data_time_ms(1));

    fields = fields((time_inx+1):end);

    if nargin < 3
        [selected_fields, dosum, scale, leftLeg, showAlt, export, lpfFreq] = selectFileds(fields);
    else
        dosum   = fieldAndOptions.dosum;
        scale   = fieldAndOptions.scale;
        leftLeg = fieldAndOptions.leftLeg;
        showAlt = fieldAndOptions.showAlt;
        export  = fieldAndOptions.export;
		lpfFreq = fieldAndOptions.lowPassFreq;
        sf_names = fieldAndOptions.selected_fields;
        selected_fields = [];
        for i = 1:length(fields)
            if any(strcmpi(sf_names, fields{i}))
                selected_fields = [selected_fields i];
            end
        end
    end

    if isempty(selected_fields)
        warning('User selected no field!');
        return;
    end
end


if leftLeg > 0
    disp('Selected left leg for hill strikes and toe offs reference');
    altLegName = 'Right';
    hill_strikes_ms = neededData.hill_strikes_l_ms;
    toe_offs_ms = neededData.toe_offs_l_ms;
    hill_strikes_a_ms = neededData.hill_strikes_r_ms;
    toe_offs_a_ms = neededData.toe_offs_r_ms;
else
    disp('Selected right leg for hill strikes and toe offs reference');
    altLegName = 'Left';
    hill_strikes_ms = neededData.hill_strikes_r_ms;
    toe_offs_ms = neededData.toe_offs_r_ms;
    hill_strikes_a_ms = neededData.hill_strikes_l_ms;
    toe_offs_a_ms = neededData.toe_offs_l_ms;
end


fs = data_fs;

if abs(data_fs - fs) > 0.001
    disp('Resample data');
    [U, D] = rat(fs/data_fs);
    otime = data{1};
    for i = 2:length(data)
        data{i} = resample(data{i}, U, D);
    end
    data{1} = linspace(otime(1), otime(end), length(data{2}));
    if iscolumn(data{2})
        v = data{1};
        data{1} = v(:);
    end
    data_time_ms = data{1} * 1000.0;
    data_fs = 1000.0 / (data_time_ms(2)-data_time_ms(1));
end

first_hs = find(hill_strikes_ms >= data_time_ms(1), 1);
last_to  = find(toe_offs_ms < data_time_ms(end), 1, 'last');

first_hs_a = find(hill_strikes_a_ms >= data_time_ms(1), 1);
last_to_a  = find(toe_offs_a_ms < data_time_ms(end), 1, 'last');

if isempty(first_hs) || isempty(last_to)
    warning('Invalid data time range!');
    return;
end


hill_strikes_ms = hill_strikes_ms(first_hs:last_to);
toe_offs_ms = toe_offs_ms(first_hs:last_to);

hill_strikes_a_ms = hill_strikes_a_ms(first_hs_a:last_to_a);
toe_offs_a_ms = toe_offs_a_ms(first_hs_a:last_to_a);

hill_strikes = floor((hill_strikes_ms - data_time_ms(1)) * fs / 1000.0);
toe_offs = floor((toe_offs_ms - data_time_ms(1)) * fs / 1000.0);

hill_strikes_a = floor((hill_strikes_a_ms - data_time_ms(1)) * fs / 1000.0);
toe_offs_a = floor((toe_offs_a_ms - data_time_ms(1)) * fs / 1000.0);

perturb_sample = floor((neededData.perturb_time_ms - data_time_ms(1)) * fs / 1000.0);
perturb_step_i = find(toe_offs < perturb_sample, 1, 'last') + 1;

if isempty(perturb_step_i)
    if perturb_sample < min(toe_offs)
        perturb_step_i = 1;
    else
        error('No perturbation step found!');
    end
end

if nargin < 3
    [selected_steps, ho, continus, prcnt, sep] = selectSteps(length(hill_strikes), perturb_step_i);
else
    ho       = fieldAndOptions.ho;
    continus = fieldAndOptions.continus;
    prcnt    = fieldAndOptions.prcnt;
    sep      = fieldAndOptions.sep;
    selected_steps = fieldAndOptions.steps + perturb_step_i;
end

if isempty(selected_steps)
    warning('User selected no step!');
    return;
end

if export > 0
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(exportDir);
end


disp([num2str(length(selected_fields)) ' field(s) selected for plotting!']);


if lpfFreq > 0
	[f_num, f_den] = butter(2, lpfFreq / (data_fs / 2));
end

if dosum > 0
    for xfi = 1:length(selected_fields)
        fi = selected_fields(xfi);
		adata = data{1+fi};
		
		if lpfFreq > 0 
			[adata, f_delay] = filter(f_num, f_den, adata);
			f_delay = length(f_delay);
			adata = adata((1 + f_delay):end);
		end
		
        if xfi == 1
            fdata = adata;
            sumtitle = fields{fi};
        else
            fdata = fdata + adata;
            sumtitle = [sumtitle '+' fields{fi}];
        end
    end
    selected_fields = selected_fields(1);
    data{1+selected_fields(1)} = fdata;
end

disp([num2str(length(selected_steps)) ' step(s) selected for plotting!']);


for xfi = 1:length(selected_fields)
    fi = selected_fields(xfi);

    if sep < 0.1 && noPlot < 0.5
        if dosum > 0
            figure('Name', sumtitle);
            title(sumtitle, 'Interpreter', 'none');
        else
            figure('Name', fields{fi});
            title(fields{fi}, 'Interpreter', 'none');
        end

        legend();
    end

    fdata = scale .* data{1+fi};
	if dosum < 0.5 && lpfFreq > 0 
		[fdata, f_delay] = filter(f_num, f_den, fdata);
		f_delay = length(f_delay);
		fdata = fdata((1 + f_delay):end);
	end
	
    for xsi = 1:length(selected_steps)
        si = selected_steps(xsi);

        exp_file = {};

        if dosum > 0
            exp_file.field = sumtitle;
        else
            exp_file.field = fields{fi};
        end


        if sep > 0 && noPlot < 0.5
            if dosum > 0
                figure('Name', sumtitle);
            else
                figure('Name', fields{fi});
            end
            title(['Step : ' num2str(si-perturb_step_i)]);
            legend();
        end

        if noPlot < 0.5
            hold on;
        end

        ssi = hill_strikes(si);
        sei = toe_offs(si);

        if si < length(hill_strikes)
            nssi = hill_strikes(si + 1);
        else
            nssi = length(fdata);
        end

        if si > 1
            psei = toe_offs(si - 1);
        else
            psei = 1;
        end

        if showAlt > 0
            hs_a = hill_strikes_a(hill_strikes_a > ssi & hill_strikes_a < sei );
            to_a = toe_offs_a( toe_offs_a > ssi & toe_offs_a < sei );
        else
            hs_a =[];
            to_a =[];
        end

        if prcnt > 0
            t = linspace(0, 100, length(ssi:sei));
			pt = linspace(-length(psei:ssi) * 100.0 / length(ssi:sei), 0, length(psei:ssi));
			nt = linspace(100, 100 + (length(sei:nssi) * 100.0 / length(ssi:sei)), length(sei:nssi));
            hs_a = (hs_a - ssi) .* 100 / length(ssi:sei);
            to_a = (to_a - ssi) .* 100 / length(ssi:sei);
            x_label = 'Percent';
        elseif ho < 0.1
            t = (ssi:sei) * 1000.0 / fs;
            pt = (psei:ssi) * 1000.0 / fs;
            nt = (sei:nssi) * 1000.0 / fs;
            hs_a = hs_a * 1000.0 / fs;
            to_a = to_a * 1000.0 / fs;
            x_label = 'Time (ms)';
        else
            t = linspace(0, length(ssi:sei) * 1000.0 / fs, length(ssi:sei));
            pt = linspace(-length(psei:ssi) * 1000.0 / fs, 0, length(psei:ssi));
            nt = linspace(length(ssi:sei) * 1000.0 / fs, length(ssi:nssi) * 1000.0 / fs, length(sei:nssi));
            hs_a = (hs_a-sei) * 1000.0 / fs;
            to_a = (to_a-sei) * 1000.0 / fs;
            x_label = 'Time (ms)';
        end

        exp_file.prcnt = prcnt;
        exp_file.ho = ho;
        exp_file.sep = sep;
        exp_file.continus = continus;
        exp_file.step = si-perturb_step_i;

        if showAlt > 0 && noPlot < 0.5
            ax = gca;
            ccoi = ax.ColorOrderIndex;
            rng = [min(fdata(ssi:sei)) max(fdata(ssi:sei))];
            df = rng(2) - rng(1);
            rng = - df * 0.1 + (rng .* 1.2);
            for i = 1:length(hs_a)
                if i == 1 && (xsi == 1 || sep > 0)
					xline(hs_a, '--b', 'DisplayName', ['HS ' altLegName]);
                else
					xline(hs_a, '--b', 'HandleVisibility','off');
                end

            end

            for i = 1:length(to_a)
                if i == 1 && (xsi == 1 || sep > 0)
					xline(to_a, '--r', 'DisplayName', ['TO ' altLegName]);
                else
					xline(to_a, '--r', 'HandleVisibility','off');
                end
            end
            ax.ColorOrderIndex = ccoi;
        end

        exp_file.alt_hs = hs_a;
        exp_file.alt_to = to_a;

        if noPlot < 0.5
            plot(t, fdata(ssi:sei), 'LineWidth', 2, 'DisplayName', ['Step : ' num2str(si-perturb_step_i)]);
        end

        exp_file.data_t = linspace(0, length(ssi:sei) * 1000.0 / fs, length(ssi:sei));
        exp_file.data_tm = t;
        exp_file.data = fdata(ssi:sei)';


        if continus > 0 && noPlot < 0.5
            ax = gca;
            ccoi = ax.ColorOrderIndex;
            plot(pt(:), fdata(psei:ssi), 'LineWidth', 2, 'LineStyle', '-.', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            plot(nt(:), fdata(sei:nssi), 'LineWidth', 2, 'LineStyle', '-.', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            ax.ColorOrderIndex = ccoi;
        end

        exp_file.prv_t = linspace(-length(psei:ssi) * 1000.0 / fs, 0, length(psei:ssi));
		exp_file.prv_tm = pt;
		exp_file.prv_data = fdata(psei:ssi)';

		exp_file.nxt_t = linspace(length(ssi:sei) * 1000.0 / fs, length(ssi:nssi) * 1000.0 / fs, length(sei:nssi));
		exp_file.nxt_tm = nt;
		exp_file.nxt_data = fdata(sei:nssi)';

        if export > 0
            mkdir([exportDir exp_file.field]);
            save([exportDir exp_file.field '/' num2str(exp_file.step)], '-struct', 'exp_file');
        end
    end
end



%% Helper functions

    function [selected_index, ho, cont, prcnt, sep] = selectSteps(fc, perturb_step_i)
        cc = 4;
        ho = 0;
        cont = 0;
        prcnt = 0;
        sep = 0;
        rc = floor((fc + cc - 1.0) / cc);
        figsize = [200, 200, cc * 100, 240 + rc * 40];
        scrsize = get(0,'ScreenSize');
        figsize(1:2) = (scrsize(3:4)-figsize(3:4))/2;
        h.f = figure('numbertitle', 'off', 'Name', 'Select Steps', 'units','pixels','position', figsize, 'toolbar','none','menu','none');
        h.l = uicontrol('style', 'text','units','pixels', 'position',[10, 190 + rc * 40, 120, 30],'string','Select steps to plot:');
        for r = 0:(rc-1)
            for c = 0:(cc-1)
                i = c + r * cc + 1;
                if i >= fc
                    break;
                end

                if i == perturb_step_i
                    txt = 'P';
                    h.c(i) = uicontrol('style','checkbox','units','pixels','string', txt, 'ForegroundColor', 'r', ...
                        'position', [20 + c * 90, 120 + (rc+1) * 40 - r * 40, 80, 35]);
                else
                    txt = num2str(abs(i - perturb_step_i));
                    h.c(i) = uicontrol('style','checkbox','units','pixels','string', txt, ...
                        'position', [20 + c * 90, 120 + (rc+1) * 40 - r * 40, 80, 35]);
                end

            end
        end

        h.p = uicontrol('style','pushbutton','units','pixels', 'position',[40,30,120,30],'string', 'OK');
        h.ho = uicontrol('style','checkbox','units','pixels','string', 'Hold On', 'position', [40, 120, 200, 35]);
        h.cont = uicontrol('style','checkbox','units','pixels','string', 'Continus', 'position', [140, 120, 200, 35], 'Value', 1);
        h.prcnt = uicontrol('style','checkbox','units','pixels','string', 'Percent', 'position', [40, 80, 200, 35]);
        h.sep = uicontrol('style','checkbox','units','pixels','string', 'Seperate', 'position', [140, 80, 200, 35]);
        selected_index = [];

        set(h.p, 'callback', @(src, event) p_call(src, event, h));
        uiwait();

        function p_call(~, ~, h)
            vals = get(h.c,'Value');
            selected_index = find([vals{:}]);
            if isempty(selected_index)
                return;
            end
            ho = get(h.ho, 'Value');
            cont = get(h.cont, 'Value');
            prcnt = get(h.prcnt, 'Value');
            sep = get(h.sep, 'Value');
            close(h.f);
        end
    end

    function [selected_index, dosum, scale, leftLeg, showAlt, export, lpfFreq] = selectFileds(fields)
        fc = length(fields);
        selected_index = [];
        dosum = 0;
        export = 0;
        lpfFreq = -1;
        scale = 1;
        leftLeg = 0;
        showAlt = 0;
        if fc > 100
            cc = 7;
        elseif fc > 50
            cc = 6;
        elseif fc > 20
            cc = 5;
        else
            cc = 4;
        end

        rc = floor((fc + cc - 1.0) / cc);
        figsize = [50, 50, cc * 200, 120 + rc * 28];
        scrsize = get(0,'ScreenSize');
        figsize(1:2) = (scrsize(3:4)-figsize(3:4))/2;
        h.f = figure('numbertitle', 'off', 'Name', 'Select Fields', 'units','pixels','position', figsize, 'toolbar','none','menu','none');
        h.l = uicontrol('style', 'text','units','pixels', 'position',[10, 70 + rc * 28, 120, 30],'string','Select fields to plot:');
        for r = 0:(rc-1)
            for c = 0:(cc-1)
                i = c + r * cc + 1;
                if i >= fc
                    break;
                end
                h.c(i) = uicontrol('style','checkbox','units','pixels','string', fields(i), ...
                    'position', [10 + c * 200, (rc + 1.8 - r) * 28, 180, 25]);
            end
        end

        h.sum = uicontrol('style','checkbox','units','pixels','string', 'Sum', 'position', [180, 40, 60, 35]);
        h.lbl = uicontrol('style','text','units','pixels','string', 'Scale :', 'position', [225, 28, 50, 35]);
        h.scale = uicontrol('style','edit','units','pixels','string', '1', 'position', [270, 42, 100, 25]);
        h.filtEnable = uicontrol('style','checkbox','units','pixels','string', 'Filter data with LPF @ freq :', 'position', [180, 12, 180, 25]);
        h.lpfFreq = uicontrol('style','edit','units','pixels','string', '6', 'position', [340, 12, 100, 25]);
        h.leftLeg = uicontrol('style','checkbox','units','pixels','string', 'Left as Reference', 'position', [385, 38, 120, 35]);
        h.showAlt = uicontrol('style','checkbox','units','pixels', 'Value', 1,'string', 'Draw Alt HS-TO', 'position', [505, 38, 100, 35]);
        h.export = uicontrol('style','checkbox','units','pixels','string', 'Export File', 'position', [615, 38, 100, 35]);
        h.p = uicontrol('style','pushbutton','units','pixels', 'position',[40,40,120,30],'string','OK');

        set(h.p, 'callback', @(src, event) p_call(src, event, h));
        uiwait();
        function p_call(~, ~, h)
            vals = get(h.c,'Value');
            selected_index = find([vals{:}]);
            if isempty(selected_index)
                return;
            end
            dosum = get(h.sum, 'Value');
            filtEnable = get(h.filtEnable, 'Value');
            hs = get(h.scale, 'String');
            leftLeg = get(h.leftLeg, 'Value');
            showAlt = get(h.showAlt, 'Value');
            export = get(h.export, 'Value');
            scale = str2double(hs);
            if filtEnable > 0
                hs = get(h.lpfFreq, 'String');
                lpfFreq = str2double(hs);
            end
            close(h.f);
        end
    end

    function [data, fields] = load_data_file(filename)
        [~,fn,ext] = fileparts(filename);
        if strcmpi(ext, '.mat')
            disp('Load MATLAB data file ...');
            if endsWith(fn, '.csv','IgnoreCase',true)
                data = load(filename);
                fields = fieldnames(data);
                data = struct2array(data);
            else
                data = struct2array(load(filename));
                fields = fieldnames(data);
            end
        elseif strcmpi(ext, '.sto')
            disp('Load STO file ...');
            data = load_sto_file(filename);
            fields = fieldnames(data);
            data = struct2cell(data);
        elseif strcmpi(ext, '.mot')
            disp('Load MOT file ...');
            data = load_sto_file(filename);
            fields = fieldnames(data);
            data = struct2cell(data);
        end
        disp('Data loaded !');
    end

    function out = load_sto_file(filename)
        [file_data,s_data]= readtext(filename, '\t', '', '', 'empty2NaN');

        % search the numerical data (in s_data.numberMask) to find when the block
        % of data starts

        a = find(abs(diff(sum(s_data.numberMask,2)))>0);
        [m,n] = size(file_data);

        % create an array with all of the data
        num_dat = [file_data{a(end)+1:end,1:sum(s_data.numberMask(a(end)+1,:),2)}];

        % reshape to put back into columns and rows
        data = reshape(num_dat,m-a(end),sum(s_data.numberMask(a(end)+1,:),2));

        % now find the column headings if there are any
        if sum(s_data.stringMask(a(end),:)) == sum(s_data.numberMask(a(end)+1,:))
            data_label = file_data(a(end),:);
            b = a(end)-1;
        else b = a(end);
        end

        % go through the data labels and find any that are duplicates (this occurs
        % in the ground reaction force data where each forceplate has the same
        % column headings) and add a number to distinguish the duplicates.

        for i = 1:length(data_label)
            tf = strcmp(data_label(i),data_label);
            c = find(tf>0);
            if length(c) > 1
                for j = 1:length(c)
                    data_label(c(j)) = cellstr([data_label{c(j)} num2str(j)]);
                end
            end
        end

        % now create the output structure with the field names from the data labels
        % and the corresponding data from the columns of the data array
        for i = 1:length(data_label)
            f_name = data_label{i};
            % find any spaces and replace with underscore
            e = findstr(f_name, ' ');
            if ~isempty(e)
                f_name(e) = '_';
            end
            e = findstr(f_name, '.');
            if ~isempty(e)
                f_name(e) = '_';
            end
            if ~isempty(str2num(f_name(1)))
                f_name = ['N' f_name];
            end
            out.(f_name) = data(:,i);
        end
    end

    function [data, result]= readtext(fname, delimiter, comment, quotes, options)

        % Read (or set to default) the input arguments.
        if((nargin < 1) || ~ischar(fname) || isempty(fname))		% Is there a file name?
            error('First argument must be a file name!');
        end
        if(nargin < 2), delimiter=	',';				end			% Default delimiter value.
        if(nargin < 3), comment=	'';					end			% Default comment value.
        if(nargin < 4), quotes=		'';					end			% Default quotes value.
        if(nargin < 5), options=	[];					end			% Default options value.

        options=		lower(options);
        op_waitbar=		~isempty(strfind(options, 'usewaitbar'));	% Do waitbar calls.
        op_numeric=		~isempty(strfind(options, 'numeric'));		% Read as numerical.
        op_textual=		~isempty(strfind(options, 'textual')) && ~op_numeric;	% Read as textual.
        op_empty=		[];											% Ignore empties, ...
        if(~isempty(strfind(options, 'empty2zero')))
            op_empty=		0;										% ... or replace by zero ...
        elseif(op_numeric || ~isempty(strfind(options, 'empty2nan')))
            op_empty=		NaN;									% ... or replace by NaN.
        end
        if(op_textual), op_empty= num2str(op_empty);	end			% Textual 'empty'.
        if(~ischar(comment) || (length(comment) > 1))
            error('Argument ''comment'' must be a string of maximum one character.');
        end
        if(~ischar(quotes) || (length(quotes) > 2))
            error('Argument ''quotes'' must be a string of maximum two characters.');
        end

        % Set the default return values.
        result.min=		Inf;
        result.max=		0;
        result.quote=	0;

        % Read the file.
        [fid, errmess]=	fopen(fname, 'r');							% Open the file.
        if(fid < 0), error(['Trying to open ' fname ': ' errmess]); end
        text=			fread(fid, 'uchar=>char')';					% Read the file.
        fclose(fid);												% Close the file.

        if(op_waitbar)
            th= waitbar(0, '(readtext) Initialising...');			% Show waitbar.
            thch=			findall(th, '-property', 'Interpreter');
            set(thch, 'Interpreter', 'none');						% No (La)TeX) formatting.
        end

        % Clean up the text.
        eol=			char(10);
        text=			strrep(text, [char(13) char(10)], eol);		% Replace Windows-style eol.
        text=			strrep(text, char(13), eol);				% Replace MacClassic-style eol.
        if(~isempty(comment))										% Remove comments.
            text=	regexprep(text, ['^' comment '[^\n]*\n'], '');	% Remove commented lines.
            text=	regexprep(text, [comment '[^\n]*'], '');		% Remove commented line endings.
        end
        if(text(end) ~= eol), text= [text eol];				end		% End string with eol, if none.

        % Find column and row dividers.
        delimiter=		strrep(delimiter, '\t', char( 9));			% Convert to one char, quicker?
        delimiter=		strrep(delimiter, '\n', char(10));
        delimiter=		strrep(delimiter, '\r', char(10));
        delimiter=		strrep(delimiter, '\f', char(12));
        if(1 == length(delimiter))									% Find column dividers quickly.
            delimS=		find((text == delimiter) | (text == eol));
            delimE=		delimS;
        elseif(isempty(regexp(delimiter, '[\+\*\?\|\[^$<>]', 'once')))	% Find them rather quickly.
            delimS=		strfind(text, delimiter);
            eols=		find(text == eol);
            delimE=		union(eols, delimS + length(delimiter) - 1);
            delimS=		union(eols, delimS);
        else														% Find them with regexp.
            [delimS, delimE]=	regexp(text, [delimiter '|' eol]);
        end
        divRow=			[0, find(text == eol), length(text)];		% Find row dividers+last.

        % Keep quoted text together.
        if(~isempty(quotes))										% Should we look for quotes?
            if((length(quotes) == 1) || (quotes(1) == quotes(2)))	% Opening char == ending.
                exclE=			find(text == quotes(1));
                exclS=			exclE(1:2:end);
                exclE=			exclE(2:2:end);
            else													% Opening char ~= closing.
                exclS=			find(text == quotes(1));
                exclE=			find(text == quotes(2));
            end
            if((length(exclS) ~= length(exclE)) || (sum(exclS > exclE) > 0))
                if(op_waitbar), close(th); 	end						% Close waitbar or it'll linger.
                error('Opening and closing quotes don''t match in file %s.', fname);
            end
            if(~isempty(exclS))										% We do have quoted text.
                if(op_waitbar), waitbar(0, th, '(readtext) Doing quotes...'); end	% Inform user.
                r=		1;
                rEnd=	length(exclS);
                n=		1;
                nEnd=	length(delimS);
                result.quote=	rEnd;
                while((n < nEnd) && (r < rEnd)) % "Remove" delimiters and newlines within quyotes.
                    while((r <= rEnd) && (delimS(n) > exclE(r))), r= r+1;	end
                    while((n <= nEnd) && (delimS(n) < exclS(r))), n= n+1;	end
                    while((n <= nEnd) && (delimS(n) >= exclS(r)) && (delimS(n) <= exclE(r)))
                        delimS(n)=	0;
                        n=			n+1;
                    end
                    if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd); end	% Update waitbar.
                end
                if(op_waitbar), waitbar(1);	end;
                delimE=	delimE(delimS > 0);
                delimS=	delimS(delimS > 0);
            end
        end
        delimS=		delimS-1;										% Last char before delimiter.
        delimE=		[1 delimE(1:end-1)+1];							% First char after delimiter.

        % Do the stuff: convert text to cell (and maybe numeric) array.
        if(op_waitbar), waitbar(0, th, sprintf('(readtext) Reading ''%s''...', fname));	end
        r=				1;
        c=				1;											% Presize data to optimise speed.
        data=			cell(length(divRow), ceil(length(delimS)/(length(divRow)-1)));
        nums=			zeros(size(data));							% Presize nums to optimise speed.
        nEnd=			length(delimS);								% Prepare for a waitbar.
        for n=1:nEnd
            temp=			text(delimE(n):delimS(n));
            data{r, c}= 	temp;									% Textual item.
            if(~op_textual), nums(r, c)= str2double(temp);	end		% Quicker(!) AND better waitbar.
            if(text(delimS(n)+1) == eol)							% Next row.
                result.min=		min(result.min, c);					% Find shortest row.
                result.max=		max(result.max, c);					% Find longest row.
                r=				r+1;
                c=				0;
            end
            c=				c+1;
            if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd);	end	% Update waitbar.
        end
        if(op_waitbar), waitbar(1);	end

        % Clean up the conversion and do the result statistics.
        if(op_waitbar), waitbar(0, th, '(readtext) Cleaning up...');	end	% Inform user.
        data=				data(1:(r-1), 1:result.max);			% In case we started off to big.
        if(~op_textual), nums= nums(1:(r-1), 1:result.max);	end		% In case we started off to big.
        if(exist('strtrim', 'builtin')), data= strtrim(data);		% Not in Matlab 6.5...
        else							 data= deblank(data);
        end
        while(all(cellfun('isempty', data(end, :))))				% Remove empty last lines.
            data=	data(1:end-1, :);
            nums=	nums(1:end-1, :);
            r=		r-1;
        end
        while(all(cellfun('isempty', data(:, end))))				% Remove empty last columns.
            data=	data(:, 1:end-1);
            nums=	nums(:, 1:end-1);
            c=		c-1;
        end
        result.rows=		r-1;
        empties=			cellfun('isempty', data);				% Find empty items.
        result.emptyMask=	empties;
        if(op_textual)
            result.numberMask=	repmat(false, size(data));			% No numbers, all strings.
            result.stringMask=	~empties;							% No numbers, all strings.
            data(empties)=		{op_empty};							% Set correct empty value.
        else
            result.numberMask=	~(isnan(nums) & ~strcmp(data, 'NaN'));	% What converted well.
            if(op_numeric)
                nums(empties)=		op_empty;						% Set correct empty value.
                data=				nums;							% Return the numeric array.
                result.stringMask=	~(empties | result.numberMask);	% Didn't convert well: so strs.
            else
                data(result.numberMask)= num2cell(nums(result.numberMask));	% Copy back numerics.
                data(empties)=		{op_empty};						% Set correct empty value.
                result.stringMask=	cellfun('isclass', data, 'char');	% Well, the strings.
            end
        end
        result.empty=		sum(result.emptyMask(:));				% Count empties.
        result.numberMask=	result.numberMask & ~result.emptyMask;	% Empties don't count.
        result.number=		sum(result.numberMask(:));				% Count numbers.
        result.stringMask=	result.stringMask & ~result.emptyMask;	% Empties don't count.
        result.string=		sum(result.stringMask(:));				% Count strings.

        if(op_waitbar), close(th);	end								% Removing the waitbar.
    end

end
