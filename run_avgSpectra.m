% run_avgSpectra
% Base function for creating LSGs

baseDir  = 'M:\Shared drives\MBARC_All\LTSAs\GofMX\MC';
outDir = 'F:\LSGs\GOM\MC\dailyAves\';

if ~isdir(outDir)
    mkdir(outDir)
end
% Set parameters for avgSpectra_noTF, 
% I don't recommend modifying these because I haven't tested any other settings
p.NA = 5;     % number of time slices (spectral averages) to read per raw file
p.tres = 1;   % time bin resolution 0 = month, 1 = days, 2 = hours
p.navepd = 5760; 
p.pflag = 1;
p.rm_fifo = 0;
p.dctype = 0;
p.av = [100 100000 10 120];	% plot axis vector

dirList = dir(baseDir);
% dirList = dir([baseDir,'\G*']);%<- Use this form to with first letter of your LTSA folders
% to avoid random non-target folders in your base directory on GDrive.

% could implement a smarter solution.
nDirs = length(dirList);
idxF = 1;
for iDir = 1:nDirs % parfor works here, unless it hits something unusual. for is safer, but slower
    if dirList(iDir).bytes == 0
        % nothing in this folder
        continue
    end
    LTSAdir = fullfile(dirList(iDir).folder,dirList(iDir).name);
    outFile = avgSpectra_noTF(LTSAdir, outDir, p);
    fprintf('Done with folder %s\n',LTSAdir)
end

favDepl = 2;%DC % favorite deployment, just counts in order of deployment start time from the 
% set of deployments in outDir. Doesn't do any matching of deployment numbers or anything

TFbaseDir = 'C:\Users\HARP\Documents\TFs';% ideally points to MBARC_TFs on google drive if you have filestream set up
tfType = 'MBARC_TFs'; % 'Wind' option is not yet implemented. MBARC_TFs option will automatically select the 
% correct TF, but if there are two versions with different dates, it will
% pick the first one. Need to make it smart enough to pick based on date in
% future.

%% Calculate mean power (typical LSG)
LSGType = 1; %1 = mean, 2 = min (attempt at noise floor estimate),
% 3 = top 10th percentile (frequency response estimate)
plotName = 'Mean power';
saveName = 'meanPower';
plotDailyLSG(outDir,plotName,saveName,LSGType,favDepl,tfType,TFbaseDir)

%% Calculate min power (lowest 5th percentile)
LSGType = 2; 
plotName = 'Min power';
saveName = 'minPower';
plotDailyLSG(outDir,plotName,saveName,LSGType,favDepl,tfType,TFbaseDir)

%% Calculate max power (top 10th percentile)
LSGType = 3; 
plotName = 'Top 10th percentile';
saveName = 'top10thPerc';
plotDailyLSG(outDir,plotName,saveName,LSGType,favDepl,tfType,TFbaseDir)



