baseDir  = 'F:\Data\LTSAs\GOM';
outDir = 'F:\Data\LTSAs\GOM\dailyAves_noTF';

if ~isdir(outDir)
    mkdir(outDir)
end
% initial changable parameters:
% rm_fifo = 0;    % remove FIFO via interpolation on spectra. 0=no, 1=yes
% fsflag = 1;     % sample rate flag for FIFO removal 1=80kHz, 0=all other
p.NA = 10;     % number of time slices (spectral averages) to read per raw file
p.tres = 1;   % time bin resolution 0 = month, 1 = days, 2 = hours
p.navepd = 5760;
p.sflag = 1;
p.pflag = 1;
p.rm_fifo = 0;
p.dctype = 0;
p.av = [100 100000 10 120];	% plot axis vector

dirList = dir([baseDir,'\G*']);
nDirs = length(dirList);
outFileList = [];
idxF = 1;
for iDir = 1:nDirs
    LTSAdir = fullfile(dirList(iDir).folder,dirList(iDir).name);
    outFile = avgSpectra_noTF(LTSAdir, outDir, p);
    if ~isempty(outFile)
        outFileList{idxF,1} = outFile;
        idxF = idxF+1;
    end
    fprintf('Done with folder %s\n',LTSAdir)
end

cAdjust = -70;
plotDailyLSG(outDir,'Mean power, no TF','meanPower_noTF',1,cAdjust)

plotDailyLSG(outDir,'Min power, no TF','minPwr_noTF',2,cAdjust)
plotDailyLSG(outDir,'75th percentile, no TF','75thPerc_noTF',3,cAdjust)