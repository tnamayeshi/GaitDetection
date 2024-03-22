function s = convertTrc2Struct(trcFile)
endHeaderLineNumber = 6;
fid = fopen(trcFile);
pnames = '';
for i = 1:4
    pnames = fgetl(fid);
end
fclose(fid);

pnames = replace(strtrim(pnames), char(9), ',');
pnames = replace(pnames, '#', '');
pnames = replace(pnames, 'Time', 'time');
pnames = replace(pnames, 'Frame', 'frame');
while contains(pnames, ',,')
    pnames = replace(pnames, ',,', ',');
end

pnames = split(pnames, ',');
pnames_arr = pnames(1:2);

for i = 3:length(pnames)
    inx = 2 + 3 * (i - 3); 
    pnames_arr(inx + 1) = strcat(pnames(i), '_X');
    pnames_arr(inx + 2) = strcat(pnames(i), '_Y');
    pnames_arr(inx + 3) = strcat(pnames(i), '_Z');
end

data = readtable(trcFile, 'FileType', 'text', 'HeaderLines', endHeaderLineNumber);
data = data(:,1:length(pnames_arr));
data.Properties.VariableNames = pnames_arr;
s = table2struct(data, 'ToScalar', true);
end