function plotDailyAveSpectra_fun(infile,ptime,mpwrtf,freq,nmave,...
        navepd,B,sflag,pflag,rm_fifo,dctype,av,tf_file)
%
% based on plotDailyAveSpectra_180119.m
%
% reads in paramter file including *DailyAves.mat file to filter and adjust
% outputs new *DailyAvesB.mat and plots
%
% original *DailyAves.mat is generated in LTSAdailySpectra.m
%

% clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load plot parameters
% [ipfile, ippath ] = uigetfile('*.txt','Choose Plot Parameter File');
% fid = fopen( fullfile(ippath,ipfile));
% if fid == -1
%     fprintf( 'Matlab couldn''t open %s', fileName );
%     return;
% end
% %loop through param file and evaluate the expression
% while ~feof( fid )
%     evalc( fgetl( fid ) );
% end
% fclose( fid );

% load daily spectral aves
% 'ptime','mpwr','mpwrTF','freq','nmave','tf_file'
% [ifile, ipath ] = uigetfile('*.mat','Choose Daily Average File to Plot');
% infile = fullfile(ipath,ifile);
% tf_file = []; % in case file doesn't have this, make empty one.
% load(infile)
% 
[ipath,inname,ext] = fileparts(infile);
name = strrep(inname,'_','\_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptime = floor(ptime);

% filter out partial days (often the first and last day with support ship
% sounds
K = find(nmave(:,1) > 0.9*navepd & nmave(:,1) <= navepd);

% remove bad days explicitly - B is often found with oneAtaTime.m
if ~isempty(B)
    [~,M] = setdiff(ptime(K),B);
    lk = length(M);
    pmp = mpwrtf(:,K(M));   % plot mean power
    YY = year(ptime(K(M)));
    MM = month(ptime(K(M)));
    pmptime = ptime(K(M));
else
    lk = length(K);
    pmp = mpwrtf(:,K);
    YY = year(ptime(K));
    MM = month(ptime(K));
    pmptime = ptime(K);
end

% apply correction for decimation low pass filter effects
switch dctype
    case 0
        disp('no correction for decimation applied')
    case 1
        load('dftf_R2013b.mat')
    case 2
        load('dftf_R2016b.mat')
    case 3  % two different decimation correction factors
        load('dftf_R2013b.mat')
        dftfa = dftf;
        load('dftf_R2016b.mat')
        dftfb = dftf;
        dctype = 0;
        for m = 1:length(K)
            if ptime(K(m)) <= dcd || ptime(K(m)) >= dcd2
                pmp(:,m) = mpwrtf(:,K(m))+ dftfa';
            else
                pmp(:,m) = mpwrtf(:,K(m))+ dftfb';
            end
        end
    case 4  % SDT HP R2016b
        load('dftf_SDT_HP_R2016b.mat')
    otherwise
        disp('Error: decimation correction flag not set properly')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove FIFO pulses if needed
if rm_fifo
    pf0 = 25;    % first FIFO pulse frequency
    npf = (1000 / pf0) - 1;  % number of fifo pulse frequencies
    for k = 1:npf
        kpf = k*pf0;
        pf = kpf:1:kpf+2; % pulse frequency indices (3 per pulse)
        x = [kpf-1 kpf+3];
        v = pmp(x,:);
        pmp(pf,:) = interp1(x,v,pf);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean differenced spectra
mmf = mean(pmp,2);   % total average of all months for these ltsas
dmf = pmp - mmf * ones(1,lk); % differences from total average

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc monthly averages

YM = 12.*YY + MM;   % year-month
nm = YM(end)-YM(1)+1;
ma = zeros(length(freq),nm);
yr = zeros(nm,1);
mn = zeros(nm,1);
me = zeros(nm,1);   % monthly effort
lstr{nm} = [];
cnt = 0;
for k = YM(1):1:YM(end)
    cnt = cnt+1;
    if dctype > 0
        ma(:,cnt) = mean(pmp(:,YM==k),2) + dftf'; % apply decimation/lpf correction factor
    else
        ma(:,cnt) = mean(pmp(:,YM==k),2); % do not apply decimation correction
    end
    yr(cnt) = floor((k-1)/12);
    mn(cnt) = k - 12.*yr(cnt);
    me(cnt) = (size(pmp(:,YM==k),2)) / eomday(yr(cnt)+2000,mn(cnt)); % effort
    if me(cnt) > 0.9
        lstr{cnt} = datestr([yr(cnt)+2000 mn(cnt) 1 0 0 0],'mmm yyyy');
    else
        lstr{cnt} = [datestr([yr(cnt)+2000 mn(cnt) 1 0 0 0],'mmm yyyy'),'*'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
LW = 3;
FS = 12;
FS2 = 16;

% for aflag = 2 case
fi = 5e3;
fe = 40e3;

% make color matrix for 12 months
mc = zeros(12,3);

mc(1,:) = [0 1 1];          % jan cyan
mc(2,:) = [0 0.75 1];       % feb sky blue
mc(3,:) = [0 0 1];          % mar blue
mc(4,:) = [.54 .17 .89];    % apr blue violet
mc(5,:) = [0.75 0 0.75];    % may purple
mc(6,:) = [1 0 1 ];         % jun magenta
mc(7,:) = [1 0 0];          % jul red
mc(8,:) = [1 0.5 0];        % aug orange
mc(9,:) = [1 1 0];         % sep yellow
mc(10,:) = [0.6 0.8 0.2];   % oct yellow-green
mc(11,:) = [0 0.75 0];      % nov green
mc(12,:) = [0 1 0];          % dec lime green

RGB = [];
for k = 1:cnt
    RGB = [RGB; mc(mn(k),:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% monthly mean
h = figure(199);
h.Position = [50 500 800 600];
set(gcf,'DefaultAxesColorOrder',RGB)
semilogx(freq,ma,'LineWidth',LW)
grid on
axis(av)
legend(lstr,'FontSize',FS)
xlabel('Frequency [Hz]','FontSize',FS2)
ylabel('Spectrum Level [dB re 1\muPa^2/Hz]','FontSize',FS2)
if 0
    title(strrep(name,'_','\_'),'FontSize',FS2)
else
    ht = text(100,0.95*(av(4)-av(3))+av(3),name,'FontSize',FS2,'interpreter','none');
    ht.HorizontalAlignment = 'center';
    ht.EdgeColor = 'k';
end
% text(100,85,strrep(dname,'_','\_'),'FontSize',FS2)
set(gca,'FontSize',FS2)

if pflag
    opfile = fullfile(ipath,[inname,'_MonthlySpectra.jpg']);
    print('-f199','-djpeg','-r300', opfile)
    opfile2 = fullfile(ipath,[inname,'_MonthlySpectra.tif']);
    print('-f199','-dtiff','-r300', opfile2)
end

if sflag
    frs = freq;
    ma = ma(1:1001,:); freq = freq(1:1001); %#ok<NASGU>
    ofile = fullfile(ipath,[inname,'_MonthlySpectra.mat']);
    save(ofile,'ma','freq','yr','mn','lstr','me')
    ofile2 = fullfile(ipath,[inname,'B.mat']);
    save(ofile2,'pmptime','pmp','freq','name','navepd','tf_file')
    freq = frs;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% daily mean
h=figure(200);
h.Position = [780 610 560 420];
semilogx(freq,pmp)
hold on
semilogx(freq,mmf,'k','LineWidth',LW)
grid on
axis(av)
% legend(LSTR)
title([name,' Spectra'],'FontSize',FS)
xlabel('Frequency [Hz]','FontSize',FS)
ylabel('Spectral Level [dB re 1\muPa^2/Hz]','FontSize',FS)
set(gca,'FontSize',FS)
hold off

if pflag
    opfile = fullfile(ipath,[inname,'_Spectra.jpg']);
    print('-f200','-djpeg','-r300', opfile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total mean differenced spectra
h=figure(201);
h.Position = [780 100 560 420];
semilogx(freq,dmf)
grid on
v = axis;
axis([av(1) av(2) v(3) v(4)])
title([name,' Differenced Spectra'],'FontSize',FS)
xlabel('Frequency [Hz]','FontSize',FS)
ylabel('Differenced Spectral Level [dB re 1\muPa^2/Hz]','FontSize',FS)
set(gca,'FontSize',FS)

if pflag
    opfile = fullfile(ipath,[inname,'_DiffSpectra.jpg']);
    print('-f201','-djpeg','-r300', opfile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percentiles
P = prctile(pmp,[1 10 50 90 99],2);
h=figure(202);
h.Position = [1330 600 560 420];
semilogx(freq,P)
grid on
v = axis;
axis([av(1) av(2) v(3) v(4)])
axis(av)
title([name,' Percentile Spectra'],'FontSize',FS)
xlabel('Frequency [Hz]','FontSize',FS)
ylabel('Spectral Level [dB re 1\muPa^2/Hz]','FontSize',FS)
set(gca,'FontSize',FS)

if pflag
    opfile = fullfile(ipath,[inname,'_PercentileSpectra.jpg']);
    print('-f202','-djpeg','-r300', opfile)
    
end

