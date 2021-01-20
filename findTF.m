function [tfNum,deplMatchIdx] = findTF(fnFileForTF, harpDataSummary)

deplMatch = [];
tfNum = [];
deplMatchIdx = [];
iFile = 1;
while isempty(deplMatch) && iFile<=size(harpDataSummary,1)
    dBaseName = strrep(harpDataSummary.Data_ID{iFile},'-','');
    deplMatch = strfind(lower(fnFileForTF),lower(dBaseName));
    iFile = iFile+1;
end
if isempty(deplMatch)
    warning('No deployment matching %s in HARP database',fnFileForTF)
    
else
    deplMatchIdx = iFile-1;
    tfNum = str2double(harpDataSummary.PreAmp{deplMatchIdx});
end