function plotDailyLSG(inDir,titleStr,saveName,typeFlag,favDepl,tfType,TFbaseDir)
%
% Single site long spectrogram use the same T and F vector.
% the power matrix for each site is pre-allocated, then filled for days
% when data exists
%
% 170707 smw
%
global PARAMS
%ns = 1; % number of sites
poff = datenum([2000 0 0 0 0 0]);   % needed to show year as 2010 not 10
spflag = 1; % save plot flag yes=1, no=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get user input
fs = 200000;
%nch = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codePath = mfilename('fullpath');
[codeDir,~] = fileparts(codePath);
harpDataSummaryCSV = fullfile(codeDir,'HARPdataSummary_20200122.csv');
harpDataSummary = readtable(harpDataSummaryCSV);

% set parameters based on sample rate
if fs == 2000 || fs == 3200
    npF = 1001; % number of frequency bins
    Fp = 0:1:1000;
    la = 10;        % plot min freq
    lb = 1000;       % plot max freq
    % color min and max in spectral power
    cmn = 39;
    cmx = 110;
    %cmin = 1; cmax = 110;   %[dB]
    fstr = 'Hz';
    ptype = 1;  % high freq linear =0, low freq log =1
elseif fs == 10000
    npF = 501; % number of frequency bins
    %Fp = 0:0.01:5;
    %la = 0;        % plot min freq
    %lb = 5;       % plot max freq
    % color min and max in spectral power
    cmn = 50;
    cmx = 80;
    % cmin = 1; cmax = 80;   %[dB]
    % fstr = 'kHz';
    ptype = 0;  % high freq linear =0, low freq log =1
    if ptype == 1
        Fp = 0:10:5000;
        la = 10;        % plot min freq
        lb = 5000;       % plot max freq
        fstr = 'Hz';
    else
        Fp = 0:0.01:5;
        la = 0;        % plot min freq
        lb = 5;       % plot max freq
        fstr = 'kHz';
        cmx = 70;   %[dB]
    end
elseif fs == 200000
    npF = 1001; % number of frequency bins
    Fp = 0:0.1:100;
    la = 0;        % plot min freq
    lb = 100;       % plot max freq
    % color min and max in spectral power
    cmn = 30;
    cmx = 60;
    
    fstr = 'kHz';
    ptype = 0;  % high freq linear =0, low freq log =1
elseif fs == 320000
    npF = 1601; % number of frequency bins
    Fp = 0:0.1:160;
    la = 0;        % plot min freq
    lb = 160;       % plot max freq
    % color min and max in spectral power
    cmn = 30;
    cmx = 60;
    %cmin = 1; cmax = 70;   %[dB]
    fstr = 'kHz';
    ptype = 0;  % high freq linear =0, low freq log =1
else
    disp(['Error, Unknown Sample Rate : ',num2str(fs)])
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% choose directory with *DailyAvesB.mat files named in increasing temporal
% order
%

fld = dir(fullfile(inDir,'*DailyAves.mat'));
x = struct2cell(fld);
fn_files = x(1,:);

% determine file order from ptime variable
startTime = [];
nFiles = size(fn_files,2);
for iFile = 1:nFiles
    load(fullfile(inDir,fld(iFile).name),'pmptime')
    startTime(iFile,1) = pmptime(1);
end

[~,fOrder] = sort(startTime);
fn_filesSorted = {fn_files{fOrder}};
% read in spectra Aves for each site, already combined over deployments
%
% [fn_files, fn_pathname] = uigetfile('*.mat','Pick Site*_SpectralAves.mat Files','MultiSelect','on');
disp('Spectrum Levels for : ')
if iscell(fn_filesSorted) % then it's a cell full of filenames
    fn = cell(length(fn_filesSorted),1);
    for k = 1:length(fn_filesSorted)
        fn{k} = fullfile(inDir, fn_filesSorted{k});
        disp(fn{k});
    end
else % it's just one file, but put into cell
    fn = cell(1,1);
    fn{1} = fullfile(inDir,fn_filesSorted);
    disp(fn);
end
nS = length(fn);    % number of *_DailyAves.mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% number of days to plot
load(fn{1})
ta = floor(pmptime(1)) + poff - 1;
load(fn{nS})
tb = ceil(pmptime(end)) + poff;
Tp = ta:1:tb;
dD = ta - poff;

nD = tb - ta + 1;

FM = zeros(npF,nD);  % pre-allocate power matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put site values in one matrix
%
sn{nS}=[];
tf_file = [];
meanPmp = [];
noiseFloor = [];
signalEstimate = [];
tfSet = [];
for k = 1:nS
    load(fn{k})   % 'T','P','freq','MP','sname'
    
    % get tf
    if strcmp(tfType,'MBARC_TFs')
        [tfNum,~] = findTF(p.fullName, harpDataSummary);
        if ~isempty(tfNum)
            [TFPath, TFName] = pick_TF_subdirs(tfNum,TFbaseDir);
        else
            TFName = [];
        end
    elseif strcmp(tfType,'Wind')
        error('Wind tf not implemented yet')  
    else
        error('Unknown TF option, please specify "MBARC_TFs" or "Wind"')  
    end
    % correct TF for more than one measurement at a specific frequency
    if isempty(TFName)
        warning('No TF file found for %d', tfNum)
        tfNum = []; % reset to empty because some files may not have this variable
        Ptf = zeros(size(pmp(:,1)));
        tf{k} = [];

    else
        
        loadTF(fullfile(TFName(1).folder,TFName(1).name))
        [C,ia,ic] = unique(PARAMS.tf.freq);
    
        tf{k} = tfNum;
        if length(ia) ~= length(ic)
            disp(['Error: TF file ',tf,' is not monotonically increasing'])
        end
        tf_freq = PARAMS.tf.freq(ia);
        tf_uppc = PARAMS.tf.uppc(ia);
        % Transfer function correction vector
        Ptf = interp1(tf_freq,tf_uppc,freq,'linear','extrap')';
    end
    tfSet(:,k) = Ptf;
    T = pmptime;
    if typeFlag == 1 % plot average
        meanPmp(k,:) = smooth(mean(pmp+Ptf,2),10);
        P = pmp+Ptf;
    elseif typeFlag == 2 % plot noise floor (5th percentile)
        meanPmp(k,:) = smooth(mean(pmpMinPwr+Ptf,2));
        P = pmpMinPwr+Ptf;
    elseif typeFlag == 3 % plot response (75th percentile)
        meanPmp(k,:) = smooth(mean(pmpercTop10+Ptf,2),10);
        P = pmpercTop10+Ptf;
    end
    pmpMinPwrSmooth = [];
    pmpercTop10Smooth = [];
    for Ir = 1:size(pmpMinPwr,2)
        pmpMinPwrSmooth(:,Ir) = smooth(pmpMinPwr(:,Ir)+Ptf,20);
        pmpercTop10Smooth(:,Ir) = smooth(pmpercTop10(:,Ir)+Ptf,20);
    end
    noiseFloor(k,:)= smooth(prctile(pmpMinPwrSmooth,5,2),20);
    signalEstimate(k,:) = smooth(mean(pmpercTop10Smooth,2),20);
    sname = [p.siteStr,'_',p.deplStr];
    %sname = [p.siteStr(2),p.deplStr];
    
    sn{k} = sname;
    deplStart{k} = T(1);
    Ti{k} = T - dD;
    FM(:,Ti{k}) = P(1:npF,:);
    
    
end
legend(gca,sn,'location','eastoutside')
grid on
% make xTick for months
d=datevec(Tp);
%takes every month just once
[a,idxm]=unique(d(:,1:2),'rows');
%takes every year just once
[a,idxy]=unique(d(:,1),'rows');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make figure
%%
hF=figure(750);clf

set(hF,'Units','normalized','Position',[.15,.2,.75,.7])
% pp = [0.25 0.25 8.0 10.5];
pp = [0.25 0.25 10.5 8.0];
set(hF,'PaperPosition',pp)
FS = 24;
FS2 = FS;
%%%%%%%%%%%%%%%%%%%%
% plot spectrogram
%
h = surf(Tp,Fp,FM,'Linestyle','none');
sgAx = get(h,'Parent');
view(sgAx,0,90)
if ptype
    set(sgAx,'yscale','log')
else
    set(sgAx,'yscale','linear')
end
set(sgAx,'TickDir','out')

axis([ta tb la lb])

set(gca,'FontSize',20)
% set(gca,'XtickLabel',[])
set(gca,'XtickLabel',Tp(idxm))
% set(gca,'xtick',Tp(idxy))
%set(gca,'xtick',Tp(idxm))
datetick('x','mmm-yy','keeplimits','keepticks')

% datetick

% % do not put 1st year label unless first month is Jan
% if d(1,2) ~= 1
%     XtL = cellstr(get(gca,'XtickLabel'));
%     if Tp(end)-Tp(1)<365
%         XtL{1} = num2str(year(Tp(1)));
%     else
%         XtL{1} = blanks(1);
%     end
%     set(gca,'XtickLabel',XtL);
% end

% xLabel = cellstr(get(gca,'XTickLabel'));
% janIdx = find(strcmp(xLabel,'Jan'));
% if ~isempty(janIdx) % if deployment is less than a year, specify year on "Jan" label, if exists.
%     xLabel{janIdx} = [ num2str(year(Tp(janIdx))), ' ',xLabel{janIdx}];
%     set(gca,'XTickLabel',xLabel)
%     xtickangle(45)
% end
ax = get(gca);
% minor ticks
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = Tp(idxm);

ap = get(gca,'Position');
set(gca,'Position',[ap(1) ap(2)+0.1 ap(3) ap(4)*.80])

ylabel(['Frequency [',fstr,']'],'FontSize',FS)

% plot frame/box
txtz = 200;
% dtl = datenum([0 0 1 0 0 0]);
dtl = datenum([0 0 0 1 0 0]);
dfa = -1e-10; dfb = 1e-10;
lx = [ta+dtl ta+dtl tb-2*dtl tb-2*dtl ta+dtl];
ly = [la+dfa lb-dfb lb-dfb la+dfa la+dfa];
lz = txtz.*ones(1,5);
hold on
plot3(lx,ly,lz,'-k','LineWidth',2)
hold off
caxis([cmn cmx]) % set color limits
colormap(jet)
z_max = max(max(get(h,'Zdata')));

hold on
% add deployment name and tf name
for iDep = 1:length(deplStart)
    text(deplStart{iDep}+poff,lb*1.06+1500*ptype,0,strrep(sn{iDep},'_','\_'))
    if ~isempty(tf{iDep})
        tfShort = num2str(tf{iDep});
        text(deplStart{iDep}+poff,lb*1.03+500*ptype,0,tfShort)
        line([deplStart{iDep}+poff,deplStart{iDep}+poff],...
            [1e-2,lb], [z_max,z_max],'color','k','LineWidth',2,'linestyle','--')
    end
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change color map to have more change over smaller range from center
% cmap = jet(cmx-cmn+1);
% nc = 100;   % number of colors
% a = 25;  % number to set miniumu color
% b = 75;
% rcc = nc/(b-a); % rate of color change ie slope
% cmap = jet(nc);
% mxc = cmap(nc,:);   % max color
% mnc = cmap(1,:);    % minimum color
% cmap(1,:) = [1 1 1];
% cmap(2:a-1,:) = ones(a-2,1)*mnc;
% cmap(a:b,:) = jet(nc -(b-a)+1);
% cmap(b+1:nc,:) = ones(nc-b,1)*mxc;
% % cmap = jet(c*(cmx-cmn+1) + b);
% % cmap = [ones(1,3); cmap];
% colormap(cmap)

% colorbar
cb = colorbar('Location','SouthOutside','FontSize',FS);%'Position',[0.12 0.125 0.8 0.0125]);
xlabel(cb,'Spectrum Level [dB re 1\muPa^2/Hz]')
%    'HorizontalAlignment','center','Units','normalized','FontSize',FS)
set(cb,'XLim',[cmn+1 cmx])

% add plus sign to last label on colorbar
if 1
    xtl = get(cb,'XTickLabel');
    xtl = [char(xtl) blanks(length(xtl))'];
    xtl(end,end) = '+';
    set(cb,'XTickLabel',xtl)
end

% use path as title - could be better...
warning('off')
hT = title(titleStr,'FontSize',FS2);
set(hT,'Position',get(hT,'Position')+[0,lb*.06+3000*ptype,0])
%%
% Save plot to file
if spflag
    opath = inDir;
    ofile = ['SingleLSG_',saveName,'.jpg'];
    print('-f750','-djpeg','-r600',fullfile(opath,ofile))
    ofile2 = ['SingleLSG_',saveName,'.tif'];
    print('-f750','-dtiff','-r600',fullfile(opath,ofile2))
    ofile3 = ['SingleLSG_',saveName,'.fig'];
    saveas(gca,fullfile(opath,ofile3))
end


%%
%favDepl = 10;
if 1
    hF=figure(751);clf
    FM2 = zeros(npF,nD);  % pre-allocate power matrix
    corrFac = signalEstimate-signalEstimate(favDepl,:);
    %corrFac = sMinusN-sMinusN(favDepl,:);

    
    for iK = 1:length(Ti)
        FM2(:,Ti{iK}) = FM(:,Ti{iK})-corrFac(iK,:)';
        %figure(347);clf
       
        %plot(freq/1000,[noiseFloor(iK,:);signalEstimate(iK,:)])
        1;
    end
    correctedTFs = tfSet-corrFac';
    figure(348);clf
    plot(freq/1000,-corrFac);grid on;legend(gca,strrep(sn,'_','\_'),'location','eastoutside')
    hold on
    plot(freq/1000,Ptf,'k')
    ylim([-20,20]);xlabel('Frequency (kHz)')
    title(sprintf('Estimated correction based on deployment %s',strrep(sn{favDepl},'_','\_')))

    
    figure(349);clf
    plot(freq/1000,correctedTFs);grid on;legend(gca,strrep(sn,'_','\_'),'location','eastoutside')
    hold on
    plot(freq/1000,Ptf,'k')
    ylim([40,100]);xlabel('Frequency (kHz)')
    title(sprintf('Estimated TFs based on deployment %s',strrep(sn{favDepl},'_','\_')))
    figure(hF)
    set(hF,'Units','normalized','Position',[.15,.2,.75,.7])
    pp = [0.25 0.25 10.5 8.0];
    set(hF,'PaperPosition',pp)
    FS = 24;
    FS2 = FS;
    %%%%%%%%%%%%%%%%%%%%
    % plot spectrogram
    h = surf(Tp,Fp,FM2,'Linestyle','none');
    sgAx = get(h,'Parent');
    view(sgAx,0,90)
    if ptype
        set(sgAx,'yscale','log')
    else
        set(sgAx,'yscale','linear')
    end
    set(sgAx,'TickDir','out')
    axis([ta tb la lb])
    set(gca,'FontSize',20)
    set(gca,'XtickLabel',Tp(idxm))
    datetick('x','mmm-yy','keeplimits','keepticks')
    ax = get(gca);
    % minor ticks
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = Tp(idxm);
    ap = get(gca,'Position');
    set(gca,'Position',[ap(1) ap(2)+0.1 ap(3) ap(4)*.80])
    ylabel(['Frequency [',fstr,']'],'FontSize',FS)
    
    % plot frame/box
    txtz = 200;
    % dtl = datenum([0 0 1 0 0 0]);
    dtl = datenum([0 0 0 1 0 0]);
    dfa = -1e-10; dfb = 1e-10;
    lx = [ta+dtl ta+dtl tb-2*dtl tb-2*dtl ta+dtl];
    ly = [la+dfa lb-dfb lb-dfb la+dfa la+dfa];
    lz = txtz.*ones(1,5);
    hold on
    plot3(lx,ly,lz,'-k','LineWidth',2)
    hold off
    caxis([cmn cmx]) % set color limits
    colormap(jet)
    z_max = max(max(get(h,'Zdata')));
    hold on
    % add deployment name and tf name
    for iDep = 1:length(deplStart)
        text(deplStart{iDep}+poff,lb*1.06+1500*ptype,0,strrep(sn{iDep},'_','\_'))
        if ~isempty(tf{iDep})
            tfShort = num2str(tf{iDep});
            text(deplStart{iDep}+poff,lb*1.03+500*ptype,0,tfShort)
            line([deplStart{iDep}+poff,deplStart{iDep}+poff],...
                [1e-2,lb], [z_max,z_max],'color','k','LineWidth',2,'linestyle','--')
        end
    end
    hold off
    % colorbar
    cb = colorbar('Location','SouthOutside','FontSize',FS);%'Position',[0.12 0.125 0.8 0.0125]);
    xlabel(cb,'Spectrum Level [dB re 1\muPa^2/Hz]')
    %    'HorizontalAlignment','center','Units','normalized','FontSize',FS)
    set(cb,'XLim',[cmn+1 cmx])
    
    % add plus sign to last label on colorbar
    xtl = get(cb,'XTickLabel');
    xtl = [char(xtl) blanks(length(xtl))'];
    xtl(end,end) = '+';
    set(cb,'XTickLabel',xtl)
    
    % use path as title - could be better...
    warning('off')
    hT = title(sprintf('%s: "Corrected" LSG based on %s',titleStr,strrep(sn{favDepl},'_','\_')),'FontSize',FS2);
    set(hT,'Position',get(hT,'Position')+[0,lb*.06+3000*ptype,0])
    
    opath = inDir;
    ofile = ['SingleLSG_Corrected',saveName,'.jpg'];
    print('-f750','-djpeg','-r600',fullfile(opath,ofile))
    ofile2 = ['SingleLSG_Corrected',saveName,'.tif'];
    print('-f750','-dtiff','-r600',fullfile(opath,ofile2))
    ofile3 = ['SingleLSG_Corrected',saveName,'.fig'];
    saveas(gca,fullfile(opath,ofile3))
end