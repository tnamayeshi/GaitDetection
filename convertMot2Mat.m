function convertMot2Mat(motFile, matFile)
endHeaderLineNumber = 0;
fid = fopen(motFile);
for i = 1:50
    line = fgetl(fid);
    if contains(line, 'endheader', 'IgnoreCase', true)
        endHeaderLineNumber = i;
        break;
    end
end
fclose(fid);
if endHeaderLineNumber < 1
   warning(['Invalid mot file : ' motFile]);
   return;
end
data = readtable(motFile, 'FileType', 'text', 'HeaderLines', endHeaderLineNumber);
save(matFile, 'data');
end