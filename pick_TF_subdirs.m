function [TFPath, TFName] = pick_TF_subdirs(tfNum,TFbaseDir)

% This is a shamefully hard-coded script for picking a TF out of the
% MBARC_TFs folder.


expectedDirList = {'100-399';'400-499';'500-599';'600-699';'700-799';
                '800-899';'900-999'};
tfRanges = [100,399;400,499;500,599;600,699;700,799;
                800,899;900,999];
                    
dirIdx = find(tfNum>tfRanges(:,1)&tfNum<tfRanges(:,2));
if isdir(fullfile(TFbaseDir,expectedDirList{dirIdx}))
    
    TFPath = fullfile(TFbaseDir,expectedDirList{dirIdx},num2str(tfNum));
    TFName =  dir(fullfile(TFPath,'*.tf'));
else
   disp('Not typical TF folder structure, searching base directory for TF.')
   TFPath = fullfile(TFbaseDir);
   TFName =  dir(fullfile(TFPath,[num2str(tfNum),'*.tf']));

end