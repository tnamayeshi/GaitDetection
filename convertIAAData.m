function convertIAAData(iaa_dir)
    
	data = struct();
	
    ver = 1;
    
    data.version = ver;
    
    files_struct = dir(fullfile(iaa_dir, '*.sto'));
	
	files = cell(1, length(files_struct));
    pnames = cell(1, length(files_struct));
	
	for i = 1:length(files_struct)
		files{i} = files_struct(i).name;
		pnames{i} = strrep(strrep(lower(files{i}), 'gaitmodel2392-scaled_inducedaccelerations_', ''), '.sto', '');
	end
    
    for i = 1:length(files)
		if strcmpi(pnames{i}, 'induced_constraint_reactions')
			% continue;
		end
       fn = fullfile(iaa_dir, files{i});
       if isfile(fn)
          tmp = convertMot2Struct(fn);
          data.(pnames{i}) = tmp;
       end
    end
    
    data.dataType = 'iaa';
    
    save(fullfile(iaa_dir, 'iaadata.mat'), '-struct', 'data');
end