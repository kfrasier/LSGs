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

favDepl = [2,3,4,5,6];%DC % favorite deployment, just counts in order of deployment start time from the 
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


%% Example of saving transfer function file *.tf
load('F:\LSGs\GOM\MC\dailyAves\SingleLSG_Corrected_meanPower.mat')
myDepl = 9; % pick the the deployment
thisCorrection = corrFac(myDepl,:); 
thisCorrection(600:end) = 0;% set correction to zero above the point where you stop beleiving it;
thisNewTF = smooth(tfSet(:,myDepl)-thisCorrection',50);% apply correction to original TF, and smooth a little.
myTF = [freq',thisNewTF];

% make a folder to store it if needed
tfOutDir = fullfile(outDir, 'tf files\');
if ~isdir(tfOutDir)
    mkdir(tfOutDir)
end
% write TF file
myTFName = fullfile(tfOutDir,sprintf('%0.0f_%s_adjusted_invSensit.tf',tf{myDepl},sn{myDepl}));
fod = fopen(myTFName,'w');
fprintf(fod,'%6.0f   %4.2f\n',myTF');
fclose(fod);
% Test that the TF file works
loadTF(myTFName)
global PARAMS
figure(11);clf
plot(PARAMS.tf.freq,PARAMS.tf.uppc)
hold on
plot(freq,tfSet(:,myDepl))
legend({'adjusted TF','original TF'})
grid on
title(PARAMS.tf.filename)
