function saveToJSON(FEMData, jsonFile);

jsonStr = jsonencode(FEMData);

jsonStr = strrep(jsonStr, ',"', sprintf(',\r\t"'));
jsonStr = strrep(jsonStr, '":', sprintf('": '));

jsonStr = strrep(jsonStr, '[[', sprintf('[\r\t\t['));
jsonStr = strrep(jsonStr, ']]', sprintf(']\r\t]'));
jsonStr = strrep(jsonStr, '],[', sprintf('],\r\t\t['));

jsonStr = strrep(jsonStr, '{', sprintf('{\r\t'));
jsonStr = strrep(jsonStr, '}', sprintf('\r}'));

fid = fopen(jsonFile, 'w');
if fid == -1 
    error('Cannot create JSON file'); 
end

fwrite(fid, jsonStr, 'char');
fclose(fid);