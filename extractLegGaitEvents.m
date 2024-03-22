function [hs, to, hs_v, to_v] = extractLegGaitEvents(grf, legVel, delta_s, mass, trial_dir, hint_msg)

do_log = 1;

if nargin < 5
    trial_dir = 'Unknown Trial';
    do_log = 0;
end

if nargin < 6
    hint_msg = '';
    do_log = 0;
end

trial_hint = [trial_dir ' : ' hint_msg];

msg_logs = {};

hsThreshold   = 8; % Threshold for minimum grf to assume start of hill strike (N)
stepThreshold = mass * 9.81 / 5; % Threshold for minimum grf to assume it's an step (N), default: 20% of weight
minStepLenghtSamples = 10; % Minimum step length (samples)
speedThreshold = 0.01;
stepDetSpeedThreshold = 0.5;

hs_v = [];
to_v = [];

a = grf > hsThreshold;
b = a(2:end) - a(1:end-1);

hs = find(b > 0);
to = find(b < 0);

if length(hs) < 2 || length(to) < 2
    hs = [];
    to = [];
    warning('No Gait Even detected!');
   return; 
end

% Remove unpaired hill strikes or toe offs

if hs(1) > to(1)
	to = to(2:end);
end

if hs(end) > to(end)
	hs = hs(1:end-1);
end

inx = 1;

hs = hs(:)';
to = to(:)';

while inx <= length(hs)
	invalid_step = 0;
	if to(inx) - hs(inx) < minStepLenghtSamples
		invalid_step = 1;
	else
		step_i = hs(inx):to(inx);
		step = grf(step_i);
		if max(step) < stepThreshold
			invalid_step = 1;
		else
			p_to = 1;
            if inx > 1
                p_to = to(inx - 1);
            end
            
            n_hs = length(grf) - 1;
			if inx < length(hs)
                n_hs = hs(inx + 1);
            end
			
            [fhs, fto] = fineTuneStep(grf, p_to, hs(inx), to(inx), n_hs);
            hs(inx) = fhs;
            to(inx) = fto;
            
            
            [fhs_v, fto_v] = fineTuneStepSpeed(legVel, hs(inx) + delta_s, to(inx) + delta_s, speedThreshold);
            
            step_v_i = fhs_v:fto_v;
            step_v = legVel(step_v_i);
            
            sx = find(step_v > stepDetSpeedThreshold, 1);
            if ~isempty(sx)
                is_valid_new_step = 1; % force to 0 so never insert step based on speed!

                n_to = fhs_v + sx - delta_s;

                ex = find(step_v > 0.01, 1, 'last');
                if ~isempty(ex)
                    n_hs = fhs_v + ex - delta_s;
                else
                    disp(['Invalid Leg Speed! Trial "' trial_hint '"']);
                    msg_logs{end + 1} = ['Invalid Leg Speed! Trial "' trial_hint '"'];
                end

                if((n_hs - n_to) < 5)
                   is_valid_new_step = 0;
                   disp(['Check this #1 Trial "' trial_hint '"']);
                   msg_logs{end + 1} = ['Check this #1 Trial "' trial_hint '"'];
                else 
                    n_step_i = n_to:n_hs;
                    n_step = grf(n_step_i);
                    if (max(n_step) - min(n_step)) < 10
                        is_valid_new_step = 0;
                        disp(['Check this #2 Trial "' trial_hint '"']);
                        msg_logs{end + 1} = ['Check this #2 Trial "' trial_hint '"'];
                    end
                end
                
                if(is_valid_new_step)
                    disp(['Check this !!! Trial "' trial_hint '"']);
                    msg_logs{end + 1} = ['Check this !!! Trial "' trial_hint '"'];
                    hs = [hs(1:inx) n_hs hs((inx+1):end)];
                    to = [to(1:(inx-1)) n_to to((inx):end)];
                end

            end
            
		end
	end
	
	
	if invalid_step > 0.1
		hs(inx) = [];
		to(inx) = [];
	else
		inx = inx + 1;
	end
end

hs_v = hs;
to_v = to;

inx = 1;
while inx <= length(hs)
    [fhs, fto] = fineTuneStepSpeed(legVel, hs(inx) + delta_s, to(inx) + delta_s, speedThreshold);
    hs_v(inx) = fhs;
    to_v(inx) = fto;
    inx = inx + 1;
end


if do_log && ~isempty(msg_logs)
    fid = fopen(fullfile(trial_dir, 'GaitLog.txt'), 'a+');
    fprintf(fid, 'Log file @ %s for %s : %s\n', datestr(now,'HH:MM:SS.FFF'), trial_dir, hint_msg);
    for i = 1:length(msg_logs)
       fprintf(fid, '%s\n', msg_logs{i});
    end
    fclose(fid);
end


function [fhs, fto] = fineTuneStep(theGRF, pto, ohs, oto, nhs)
    dv_th = 0.5;
    
	grf_dv = theGRF((pto+1):ohs) - theGRF(pto:(ohs-1));
    
    x = find(grf_dv <= dv_th, 1, 'last');
    
    if isempty(x)
       fhs = ohs;
    else
       fhs = ohs - (length(grf_dv) - x) - 1;
    end
    
    grf_dv = theGRF((oto+1):nhs) - theGRF(oto:(nhs-1));
    x = find(grf_dv >= -dv_th, 1);
    
    if isempty(x)
       fto = oto;
    else
       fto = oto + x - 2;
    end
end

function [fhs, fto] = fineTuneStepSpeed(theV, ohs, oto, dv_th)
    
	v_dv = theV(ohs:oto);
    
    x = find(v_dv <= -dv_th, 1);
    
    if isempty(x)
       fhs = ohs;
       fto = oto;
       warning('Invalid Leg Speed !!');
    else
       fhs = ohs + x;
       x = find(v_dv <= -dv_th, 1, 'last');
       fto = oto - (length(v_dv) - x);
    end
    
end
end














