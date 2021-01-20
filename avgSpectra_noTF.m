function outFile = avgSpectra_noTF(LTSAdir, outDir, p)

% Calculate average LTSA values WITHOUT TRANSFER FUNCTION
% given an input folder containing .ltsa files
% Save output somewhere
% Ugly hack of SMW's code, should pull stuff into functions etc.
% only works for 200khz data, full band LTSAs

global PARAMS
outFile = [];

% Figure out pieces of deployment name with a few possible matching options
[~,fullProjectName] = fileparts(LTSAdir);
p.fullName = fullProjectName;
% determine name parts
if ~isempty(regexp(fullProjectName,'(\w*)_(\w*)_(\d*)','match'))
    splitTemp = split(fullProjectName,'_');
    p.projectStr = splitTemp{1};
    p.siteStr = splitTemp{2};
    p.deplStr = splitTemp{3};
elseif ~isempty(regexp(fullProjectName,'(\w*)_(\w*)(\d*)','match'))
    splitTemp = split(fullProjectName,'_');
    p.projectStr = splitTemp{1};
    p.siteStr = splitTemp{2}(1:end-2);
    p.deplStr = splitTemp{2}(end-1:end);
elseif ~isempty(regexp(fullProjectName,'(\w*)(\d*)','match'))
    p.projectStr = fullProjectName(1:end-3);
    p.deplStr = fullProjectName(end-2:end-1);
    p.siteStr = fullProjectName(end);
else
    error('LTSA folder name does not match a known name format eg. GOM_MC_01, or GofMX_MC01, or SOCAL31M')
end
    
% initial changable parameters:
% rm_fifo = 0;    % remove FIFO via interpolation on spectra. 0=no, 1=yes
% fsflag = 1;     % sample rate flag for FIFO removal 1=80kHz, 0=all other
% p.NA = 5;     % number of time slices (spectral averages) to read per raw file
% p.tres = 1;   % time bin resolution 0 = month, 1 = days, 2 = hours
% 
% p.navepd = 5760;
% p.pflag = 1; 
% p.rm_fifo = 0;
% p.dctype = 0;
mnum2secs = 24*60*60;
% p.av = [100 100000 10 120];	% plot axis vector  

sprintf('Processing folder %s\n',LTSAdir)
% get LTSA file names
fList = dir(fullfile(LTSAdir,'*.ltsa'));

% sort if multiple files selected
fn_files = {fList.name}';
if size(fn_files,1)>1
    fn_files = sort(fn_files);
end
fn_pathname = LTSAdir;

[~,dName] = fileparts(LTSAdir);
outName = [dName '_DailyAves.mat'];
outFile = fullfile(outDir,outName);

disp('Calculating Averages for:')
if iscell(fn_files) % then it's a cell full of filenames
    fn = cell(length(fn_files),1);
    for k = 1:length(fn_files)
        fn{k} = fullfile(fn_pathname, fn_files{k});
        disp(fn{k});
    end
else % it's just one file, but put into cell
    fn = cell(1,1);
    fn{1} = fullfile(fn_pathname,fn_files);
    disp(fn);
end


tic
nltsas = length(fn);    % number of LTSA files

% read ltsa headers and sum the total number of raw files for
% pre-allocating vectors/matrices
nrftot = 0;
% fill up header matrix H
H = zeros(nrftot,3);    % 3 columns: filenumber, datenumber, byteloc in filenumber
cnt1 = 1;
cnt2 = 0;
badLTSA = zeros(nltsas,1);

for k = 1:nltsas    % loop over files
    PARAMS.ltsa = [];   % clear
    PARAMS.ltsa.ftype = 1;
    [PARAMS.ltsa.inpath,infile,ext] = fileparts(fn{k});
    PARAMS.ltsa.infile = [infile,ext];
    try
        read_ltsahead % better would be to just read nrftot instead of whole header
    catch
        continue
    end
    nrf = PARAMS.ltsa.nrftot;
    nrftot = nrftot + nrf;
    nf = PARAMS.ltsa.nf;
    if k == 1 % read and set some useful parameters that should be the same across all ltsas
        
        nave = PARAMS.ltsa.nave(1);
        freq = PARAMS.ltsa.freq;
    end

    fs0 = PARAMS.ltsa.fs;
    if fs0 ~=200000
        disp('not 200kHz ltsa, skipping');
        badLTSA(k) = 1;
        continue
    end
    knave = PARAMS.ltsa.nave;
    cnt2 = cnt2 + nrf;
    H(cnt1:cnt2,1) = k.*ones(nrf,1);
    H(cnt1:cnt2,2) = PARAMS.ltsa.dnumStart;
    H(cnt1:cnt2,3) = PARAMS.ltsa.byteloc;
    cnt1 = cnt2 + 1;
end
fn(logical(badLTSA)) = [];
% % for removing 80kHz FIFO
% if max(knave) == 38
%     fsflag = 1;
% else
%     fsflag = 0;
% end
disp('Done collecting LTSA header info')
toc
dvec = datevec(H(:,2));

if p.tres == 0
    mnum = dvec(:,1).*12 + dvec(:,2);   % month number where 1 = Jan 2000
elseif p.tres == 1
    mnum = floor(datenum(dvec));   % day number where 1 = Jan 2000
elseif p.tres == 2
    mnum = floor(datenum(dvec) * 24);   % hour
else
    disp(['Error: unknown time resolution = ',num2str(tres)])
end
%
dur = unique(mnum);   % unique averaging time bins
nm = length(dur); % number of ave time bins

% pwrA = zeros(nf,cnt2);
cnt1 = 1;
cnt2 = 0;
ptime = zeros(nm,1);      % start time of bin average
nmave = zeros(nm,2);    % number of averages possible and used(ie not filtered out) for each bin
mpwr = ones(nf,nm);     % mean power over time period
minPwr = ones(nf,nm); 
perc5 = ones(nf,nm); 
perc75= ones(nf,nm);
percTop10 = ones(nf,nm);

for m = 1:nm    % loop over time average bins
    if p.tres == 0
        disp(['Month = ',num2str(dur(m))])
    elseif p.tres == 1
        disp(['Day = ',num2str(dur(m))])
    elseif p.tres == 2
        disp(['Hour = ',num2str(dur(m))])
    end
    
    I = [];
    I = find(mnum == dur(m));
    nrfM = length(I);
    pwrM = [];
    pwrM = nan(nf,p.NA*nrfM);
    fnum = [];
    fnum = unique(H(I,1));
    nfiles = length(fnum);
    NBO = 0;    % number of averages (taves) read for mean spectra
    for f = 1:nfiles    % loop over files with same month (probably only 2 max)
        % open ltsa file
        % fid = fopen([PARAMS.ltsa.inpath,fn{fnum(f)}],'r');
        fid = fopen(fn{fnum(f)},'r');
        % samples to skip over in ltsa file
        %J = [];
        J = find(H(I,1) == fnum(f));
        nrfRead = length(J);
        skip = H(I(J(1)),3);    % get first byteloc of file for that month
        fseek(fid,skip,-1);    % skip over header + other data
        
        ptime(m) = H(I(J(1)),2);     % save 1st time of this time unit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % number of time slices to use for calcs
        % need to set 'xxxx*int8' value
        prcsn = [num2str(p.NA*nf),'*int8'];
        
        % allocate memory
        %pwrF = [];
        NB = p.NA * nrfRead;
        %pwrF = zeros(nf,NB);
        
        % number of time slices to skip on each read
        SA = nave - p.NA;
        skip = nf*(SA);
        % initial time slice to start at
        if fs0 == 2000 || fs0 == 200000 || fs0 == 10000
            %         IA = (15+5)*nf;
            IA = (nave + p.NA)*nf;  % this should work better for 320kHz vs commented 200kHz hardwired above
        elseif fs0 == 3200 || fs0 == 320000
            IA = (nave + p.NA-1)*nf;
        else
            disp('Error: unknown sample rate')
            disp(['fs0 = ',num2str(fs0)])
        end
        fseek(fid,IA,-0);
        
        % read data into File power
        [pwrF,count] = fread(fid,[nf,NB],prcsn,skip);
        if count ~= nf*NB
            disp('error = did not read enough data')
            %             NBO = NB;
            NB = floor(count/nf);
            %             disp(['NBO = ',num2str(NBO)])
            disp(['NB = ',num2str(NB)])
        end
        
        % fill up Month power
        if nfiles == 1
            pwrM = pwrF;
            NBO = NB;
        else
            if f == 1
                pwrM(1:nf,1:NB) = pwrF;
                NBO = NB;
            else
                pwrM(1:nf,NBO+1:NBO+NB) = pwrF;
                NBO = NBO+NB;
            end
        end
        disp(['Number of 5s time bins in this Average Time Bin = ',num2str(NBO)])
        fclose(fid);
    end  % end for f    
    nmave(m,1) = NBO;
    cnt2 = cnt2 + size(pwrM,2);
    
    cnt1 = cnt2 + 1;
    mpwr(1:nf,m) = nanmean(pwrM,2);
    prcPwr = prctile(pwrM',[1,5,75]);
    minPwr(1:nf,m) = nanmin(pwrM,[],2);
    perc5(1:nf,m) = prcPwr(2,:)';
    perc75(1:nf,m) =  prcPwr(3,:)';
    percTop10(1:nf,m) = -prctile(-pwrM',10);

end

ptime = floor(ptime);

% try to remove bad data points with some heuristics
outlierIdx0 = mode(minPwr)==0;
nmave(outlierIdx0,1) = 0;
if ~isempty(outlierIdx0)
    1;
end
minPwrTemp = minPwr(:,~outlierIdx0);
outlierIdx = find(mean(minPwr)<(mean(mean(minPwrTemp))-2*(std(mean(minPwrTemp))))|...
    mean(minPwr)>(mean(mean(minPwrTemp))+3*(std(mean(minPwrTemp)))));

if ~isempty(outlierIdx)
    1;
end
nmave(outlierIdx,1) = 0;

outlierIdx = find(mode(minPwr)==0);
% filter out partial days (often the first and last day with support ship
% sounds
K = find(nmave(:,1) > 0.9*p.navepd & nmave(:,1) <= p.navepd);

lk = length(K);
pmp = mpwr(:,K);
%YY = year(ptime(K));
%MM = month(ptime(K));
pmptime = ptime(K);
pmperc75 = perc75(:,K);
pmperc5 = perc5(:,K);
pmpMinPwr = minPwr(:,K);
pmpercTop10 = percTop10(:,K);

% save results
save(outFile,'lk','freq','nmave','p','pmpMinPwr','pmperc75','pmperc5','pmpercTop10',...
    'pmptime','pmp')

t = toc;
disp(' ')
disp(['Time Elapsed: Spectra from LTSA ', num2str(t),' secs'])

