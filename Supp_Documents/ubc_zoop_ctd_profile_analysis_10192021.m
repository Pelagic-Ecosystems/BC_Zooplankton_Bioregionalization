% Creating the envidata table from CTD and bottle profiles
% Calculating MLD and N2 from CTD profile data
% May 17, 2021

% Use this instead of the ubc_zoop_netmetadata_10182021.m code.
% After this code, merge the satellite data.

% Updated October 19, 2021 to extract the CTD and bottle data to be used in
% the zooplankton bioregionalization analysis.

%% Load CTD data but isolate only keys that match the Net keys
load('netctd_matching_03302021.mat','netmetable','netmetablehead') % only need the net-ctd matching section, disregard other measured data
[ctdkeys, ia, ~] = unique(cellfun(@num2str,netmetable(:,30),'uni',0));
ctdkeys = ctdkeys(2:end);
ia = ia(2:end);
ctdyears = cell2mat(netmetable(ia,32));
netctdkeys = netmetable(ia,1);
nettowdepth = cell2mat(netmetable(ia,13));

% Save table of net keys matching with ctd keys and bottle keys
fol = 'C:\Sync\Zooplankton-data-files\';
keylist = netmetable(:,[1,30,41]);

T1 = array2table(keylist,'VariableNames',...
    {'NetKey','CTDKey','BottleKey'});
writetable(T1,strcat(fol,'net_ctd_bot_keys_10192021.xlsx'));

%% Extracting CTD data
% Modified from ubc_zoop_netmetadata.m

% values are T90 - ITS-90 Potential temperature, and SP is practical
% salinity (PSS-78) ... this need to calculate CT and SA

% Aggregate all CTD tows to use first
PT = nan(length(ia), 8000);
SP = nan(length(ia), 8000);
depth = nan(length(ia), 8000);
pres = nan(length(ia), 8000);
oxyg = nan(length(ia), 8000);
lat = nan(length(ia), 1);
lon = nan(length(ia), 1);
dtstart = nan(length(ia), 1);
CTDkeys = {};

m = 1; % starting index of large data matrix
n = 0; % end index

for y = 1980:2018
   ff = strcat('C:\Owncloud_EOAS\Data\Margolin\CIOOS\ctd_data_v2.2_',num2str(y),'.mat');
   ctd = load(ff);
   ctd = ctd.data_v2;
   
   ii = find(ctdyears == y);
   [jj, ja, ~] = intersect( cellstr(ctd.filename), ctdkeys(ii) ) ;
   CTDkeys = cat(1, CTDkeys, jj);
   
   n = m + length(ja) -1;
   ndepth = ctd.depth(ja,:);
   ndepth(isnan(ndepth)) = 0;
   dcols = find( sum(ndepth,1) > 0 );
   
   depth(m:n, dcols) = ctd.depth(ja,dcols);
   pres(m:n, dcols)  = ctd.pres(ja,dcols);
   PT(m:n, dcols)  = ctd.t90(ja,dcols);
   SP(m:n, dcols)  = ctd.SP(ja,dcols);
   if ~isempty(ctd.oxyg) && nansum(nansum(ctd.oxyg(ja,:))) > 0
       oxyg(m:n, dcols)  = ctd.oxyg(ja,dcols);
   else
       disp(strcat("No oxygen: ",num2str(y)));
   end
   lat(m:n)  = ctd.lat(ja);
   lon(m:n)  = ctd.lon(ja);
   dtstart(m:n)  = ctd.datenum_start(ja);
   
   m = n + 1;
end


%% Interpolate to regular depths
% Set max depth 4800
rdepth = nan(length(ia),4800);
rpres = nan(length(ia),4800);
rPT = nan(length(ia),4800);
rSP = nan(length(ia),4800);
roxyg = nan(length(ia),4800);
rTScoverage = nan(length(ia),1);
rOxycoverage = nan(length(ia),1);

% 1-d linear interpolation - Need to loop individual tows
for i = 1:length(ia)
   ii = ~isnan(depth(i,:));
   [dep,ii,~] = unique(depth(i,ii));
   md = floor( max(dep) );
   
   rdepth(i,1:md) = 1:md;
   rpres(i,1:md) = interp1(dep,pres(i,ii),1:md);
   rPT(i,1:md) = interp1(dep,PT(i,ii),1:md);
   rSP(i,1:md) = interp1(dep,SP(i,ii),1:md);
   roxyg(i,1:md) = interp1(dep,oxyg(i,ii),1:md);
   rTScoverage(i) = sum(~isnan( rPT(i,1:md) ));
   rOxycoverage(i) =sum(~isnan( roxyg(i,1:md) ));
end

%% Calculate MLD from the original profiles (strict to d.B. Montegut et al)
% Code is based on code provided in:
% Holte, James and Lynne Talley, 2009: A new algorithm for finding mixed 
%  layer depths with applications to Argo data and Subantarctic Mode Water 
%  formation, Journal of Atmospheric and Oceanic Technology.
% This is based on the de Boyer Montegut et al. 2004 threshold criteria
% with 0.03 kg/m^3 or temperature difference of 0.2 degrees C. The value
% closest to 10 dbar is the reference value. Interpolation is done to
% exactly match the threshold criteria.
% This uses the teos 10 gsw toolkit rather than the SW toolkit

% Analyze by cast
mldepthdens = nan(length(CTDkeys),1);
mldepthdens05 = nan(length(CTDkeys),1);
mldepthdens125 = nan(length(CTDkeys),1);
mldepthptmp = nan(length(CTDkeys),1);

for i = 1:length(CTDkeys)
    presnow = pres(i,:);
    
    % SA, CT, rho
    SA = gsw_SA_from_SP(SP(i,:),presnow,lon(i),lat(i));
    CT = gsw_CT_from_t(SA,PT(i,:),presnow);
    rho = gsw_rho(SA,CT,presnow); %in-situ density
    potrho = gsw_rho(SA,CT,0); % potential density ref z = 0

    % if density can not be calculated because only temperature data is
    % present, mark MLD as NAN
    if isnan( max(potrho) )
       disp(strcat(num2str(i)," : ",datestr(dtstart(i))," no density profile"));
       continue
    end
    
    % Calculate the index of the reference value
    starti = min(find((presnow-10).^2==min((presnow-10).^2)));
    m = find(isnan(presnow(starti+1:end)),1,'first') + starti - 1;
    presnow = presnow(starti:m);
    sal = SA(starti:m);
    temp = CT(starti:m);
    pden = gsw_rho(sal,temp,0)-1000;
    
    % MLD can not be calculated when there are not enough depth readings
    if length(presnow) == 1
       continue 
    end
    
    % Remove rare instances of NaN in the middle of the profile
    oo = find(~isnan(pden));
    presnow = presnow(oo);
    sal = sal(oo);
    temp = temp(oo);
    pden = pden(oo);
    starti = 1;
    m = length(sal);
    
    % Threshold: (0.03 kg/m^3)
    % Search for the first level that exceeds the potential density threshold
    mldepthdens(i) = m;
    for j = starti:m
        if abs(pden(starti)-pden(j))>.03
            mldepthdens(i) = j;
            break;
        end
    end
    % Interpolate to exactly match the potential density threshold 
    clear pdenseg presseg presinterp pdenthreshold
    presseg = [presnow(mldepthdens(i)-1) presnow(mldepthdens(i))];
    pdenseg = [pden(starti)-pden(mldepthdens(i)-1) pden(starti) - pden(mldepthdens(i))];
    P = polyfit(presseg,pdenseg,1);
    presinterp = presseg(1):.25:presseg(2);
    pdenthreshold = polyval(P,presinterp);
    % The potential density threshold MLD value:
    mldepthdens(i) = presinterp(max(find(abs(pdenthreshold)<.03)));

    % Threshold: (0.05 kg/m^3)
    % Search for the first level that exceeds the potential density threshold
    mldepthdens05(i) = m;
    for j = starti:m
        if abs(pden(starti)-pden(j))>.05
            mldepthdens05(i) = j;
            break;
        end
    end
    % Interpolate to exactly match the potential density threshold 
    clear pdenseg presseg presinterp pdenthreshold
    presseg = [presnow(mldepthdens05(i)-1) presnow(mldepthdens05(i))];
    pdenseg = [pden(starti)-pden(mldepthdens05(i)-1) pden(starti) - pden(mldepthdens05(i))];
    P = polyfit(presseg,pdenseg,1);
    presinterp = presseg(1):.25:presseg(2);
    pdenthreshold = polyval(P,presinterp);
    % The potential density threshold MLD value:
    mldepthdens05(i) = presinterp(max(find(abs(pdenthreshold)<.05)));

    % Threshold: (0.125 kg/m^3)
    % Search for the first level that exceeds the potential density threshold
    mldepthdens125(i) = m;
    for j = starti:m
        if abs(pden(starti)-pden(j))>.125
            mldepthdens125(i) = j;
            break;
        end
    end
    % Interpolate to exactly match the potential density threshold 
    clear pdenseg presseg presinterp pdenthreshold
    presseg = [presnow(mldepthdens125(i)-1) presnow(mldepthdens125(i))];
    pdenseg = [pden(starti)-pden(mldepthdens125(i)-1) pden(starti) - pden(mldepthdens125(i))];
    P = polyfit(presseg,pdenseg,1);
    presinterp = presseg(1):.25:presseg(2);
    pdenthreshold = polyval(P,presinterp);
    % The potential density threshold MLD value:
    mldepthdens125(i) = presinterp(max(find(abs(pdenthreshold)<.125)));

    
    % Search for the first level that exceeds the temperature threshold
    % (0.2 C)
    mldepthptmp(i) = m;
    for j = starti:m
        if abs(temp(starti)-temp(j))>.2
            mldepthptmp(i) = j;
            break;
        end
    end
    
    % Interpolate to exactly match the temperature threshold
    clear tempseg presseg presinterp tempthreshold
    presseg = [presnow(mldepthptmp(i)-1) presnow(mldepthptmp(i))];
    tempseg = [temp(starti)-temp(mldepthptmp(i)-1) temp(starti) - temp(mldepthptmp(i))];
    P = polyfit(presseg,tempseg,1);
    presinterp = presseg(1):.25:presseg(2);
    tempthreshold = polyval(P,presinterp);
    
    % The temperature threshold MLD value:
    mldepthptmp(i) = presinterp(max(find(abs(tempthreshold)<.2)));
end

% Mark as NAN the impossibly large value
oo = find(mldepthdens > 1000);

mldepthdens(oo) = nan;
mldepthdens05(oo) = nan;
mldepthdens125(oo) = nan;
mldepthptmp(oo) = nan;

figure; 
subplot(2,2,1)
histogram(mldepthdens); xlabel("0.03 kg/m^3");
subplot(2,2,2)
histogram(mldepthptmp); xlabel("0.2 C");
subplot(2,2,3)
histogram(mldepthdens05); xlabel("0.05 kg/m^3");
subplot(2,2,4)
histogram(mldepthdens125); xlabel("0.125 kg/m^3");

figure;
subplot(1,3,1);
scatter(mldepthdens,mldepthdens05); xlabel("0.03 kg/m^3"); ylabel("0.05 kg/m^3");
subplot(1,3,2);
scatter(mldepthdens,mldepthdens125); xlabel("0.03 kg/m^3"); ylabel("0.125 kg/m^3");
subplot(1,3,3);
scatter(mldepthdens,mldepthptmp); xlabel("0.03 kg/m^3"); ylabel("0.2 C");

%% Derived CTD variables

% SA, CT, rho
SA = gsw_SA_from_SP(rSP,rpres,lon,lat);
CT = gsw_CT_from_t(SA,rPT,rpres);
rho = gsw_rho(SA,CT,rpres); %in-situ density
potrho = gsw_rho(SA,CT,0); % potential density ref z = 0

% Spiceness from 0 dbar, N2
% Warmer, saltier water is more spicy while cooler, less salty water is more minty
% Low N2, low stratification/mixing, high N2, more stratification
spice = gsw_spiciness0(SA, CT);
[N2, p_mid] = gsw_Nsquared(SA, CT, rpres, lat); % Negative N2 is unstable stratification

% Calculate N2 as (-g/rho)(drho/dz)
N2v2 = nan(length(lon),4800);
dRHO = nan(length(lon),4800);
for i = 1:length(lon)
    md = max(rdepth(i,:));
    
    drho = rho(i,1:md-1) - rho(i,2:md);
    dz = rdepth(i,2:md) - rdepth(i,1:md-1);
    N2v2(i,2:md) = (-9.7963 ./ rho(i,2:md) ) .* (drho./dz);
    dRHO(i,2:md) = -drho;
end

% use gsw_SA_CT_plot to view profiles
gsw_SA_CT_plot(SA,CT,0);

%% Vertical integration

% Create table of CTD Data
mldtable = cat(2,mldepthdens,mldepthdens05,mldepthdens125,mldepthptmp);

% Integrate CTD values until 1-10m, 1-50m, 1-net tow depth, 1-10/1-MLD, 
ctdinteg = [];
for i = 1:size(CT,1)
    imld = round( mldepthdens(i) );
    if isnan(imld) % For 4 instances when there is no MLD because no salinity data is given, use mean MLD of 15m
       imld =  15;
    end
    netd = nettowdepth(i);
    
    % if density can not be calculated because only temperature data is
    % present, mark MLD as NAN
    if isnan( max(potrho) )
       disp(strcat(num2str(i)," : ",datestr(dtstart(i))," no density profile"));
       continue
    end
    
    ctdinteg = cat(1,ctdinteg, ...
     [nanmean(CT(i,1:10)),nanmean(CT(i,1:50)),nanmean(CT(i,1:netd)),nanmean(CT(i,1:imld)),...
     nanmean(SA(i,1:10)),nanmean(SA(i,1:50)),nanmean(SA(i,1:netd)),nanmean(SA(i,1:imld)),... 
     nanmean(potrho(i,1:10)),nanmean(potrho(i,1:50)),nanmean(potrho(i,1:netd)),nanmean(potrho(i,1:imld)),...
     nanmean(spice(i,1:10)),nanmean(spice(i,1:50)),nanmean(spice(i,1:netd)),nanmean(spice(i,1:imld)),...
     nanmean(roxyg(i,1:10)),nanmean(roxyg(i,1:50)),nanmean(roxyg(i,1:netd)),nanmean(roxyg(i,1:imld)),...
     nanmean(N2v2(i,1:10-1)),nanmean(N2v2(i,1:30-1)), nanmean(N2v2(i,1:50-1)),...
     nanmean(N2v2(i,1:netd-1)),nanmean(N2v2(i,1:imld-1))] );
end


%% Create datatable of CTD values
% [ctdkeys, ...]

% Merge the MLD values
ctdtable = cat(2,mldtable,ctdinteg);
ctdtable = cat(2,ctdtable,[rTScoverage,rOxycoverage]);
ctdtable = cat(2, ctdkeys, num2cell(ctdtable));

ctdintvars = {'CTDKey',...
    'MLD003','MLD005','MLD125','MLDTemp',...
    'TempI10','TempI50','TempINet','TempIMLD',...
    'SaliI10','SaliI50','SaliINet','SaliIMLD',...
    'DensI10','DensI50','DensINet','DensIMLD',...
    'SpiceI10','SpiceI50','SpiceINet','SpiceIMLD',...
    'OxygI10','OxygI50','OxygINet','OxygIMLD',...
    'N2V2I10','N2V2I30','N2V2I50','N2V2INet','N2V2IMLD',...
    'CoverageTS','CoverageOxy'};

T2 = array2table(ctdtable,'VariableNames',ctdintvars);
writetable(T2,strcat(fol,'ctdnet_integrated_mld_10192021.xlsx'));



%%  Extract bottle data similar to CTD approach 
% note that there are 14 instances when the netkey was listed but not
% located in the bottle files. Thus, the size of the table needs to be
% adjusted and BOTkeys needs to be followed. This is done during the
% binning section in the next chunk.

% clear
% load('netctd_matching_03302021.mat','netmetable','netmetablehead') 

[botkeys, ib, ~] = unique(cellfun(@num2str,netmetable(:,41),'uni',0));
botkeys = botkeys(2:end);
ib = ib(2:end);
botyears =  cell2mat(netmetable(ib,43));

ntwdep = cell2mat(netmetable(ib,13));


% Aggregate all CTD tows to use first
nitr = nan(length(ib),25);
sili = nan(length(ib),25);
phos = nan(length(ib),25);
botPT = nan(length(ib), 25);
botSP = nan(length(ib), 25);
botoxyg = nan(length(ib), 25);
botdepth = nan(length(ib), 25);
botpres = nan(length(ib), 25);
botlat = nan(length(ib), 1);
botlon = nan(length(ib), 1);
botdtstart = nan(length(ib), 1);
BOTkeys = {};

m = 1; % starting index of large data matrix
n = 0; % end index

for y = 1980:2018
   ff = strcat('C:\Owncloud_EOAS\Data\Margolin\CIOOS\bot_data_v2.2_',num2str(y),'.mat');
   bot = load(ff);
   bot = bot.data_v2;
   
   ii = find(botyears == y);
   [jj, ja, ~] = intersect( cellstr(bot.filename), botkeys(ii) ) ;
   
   BOTkeys = cat(1, BOTkeys, jj);
   
   n = m + length(ja) -1;
   ndepth = bot.depth(ja,:);
   ndepth(isnan(ndepth)) = 0;
   dcols = find( sum(ndepth,1) > 0 );
   
   botdepth(m:n, dcols) = bot.depth(ja,dcols);
   botpres(m:n, dcols)  = bot.pres(ja,dcols);
   botPT(m:n, dcols)  = bot.t90(ja,dcols);
   botSP(m:n, dcols)  = bot.SP(ja,dcols);
   if ~isempty(bot.oxyg) && nansum(nansum(bot.oxyg(ja,:))) > 0
       botoxyg(m:n, dcols)  = bot.oxyg(ja,dcols);
   else
       disp(strcat("No oxygen: ",num2str(y)));
   end
   if ~isempty(bot.nitr) && nansum(nansum(bot.nitr(ja,:))) > 0
       nitr(m:n, dcols)  = bot.nitr(ja,dcols);
       sili(m:n, dcols)  = bot.sili(ja,dcols);
       phos(m:n, dcols)  = bot.phos(ja,dcols);
   else
       disp(strcat("No nutrients: ",num2str(y)));
   end
   botlat(m:n)  = bot.lat(ja);
   botlon(m:n)  = bot.lon(ja);
   botdtstart(m:n)  = bot.datenum_start(ja);
   
   m = n + 1;
end

%% Bin interpolation to 1m intervals 
botdep = 1:4500;

% only interpolate to round(max bottle tow depth)
rbotdepth = nan(length(BOTkeys), length(botdep) );
rbotpres = nan(length(BOTkeys), length(botdep) );
rbotPT = nan(length(BOTkeys), length(botdep) );
rbotSP = nan(length(BOTkeys), length(botdep) );
rbotoxyg = nan(length(BOTkeys), length(botdep) );
rnitr = nan(length(BOTkeys), length(botdep) );
rphos = nan(length(BOTkeys), length(botdep) );
rsili = nan(length(BOTkeys), length(botdep) );
botnettowdep = nan(length(BOTkeys),1);

c = 0;
emptybotfile = [];
% 1-d linear interpolation - Need to loop individual tows
for i = 1:length(BOTkeys)
    botnettowdep(i) = ntwdep(i);
        
    ii = find(~isnan(botdepth(i,:)));
    [~,ui,~] = unique(botdepth(i,ii));
    ii = ii(ui);
    if isempty(ii) % some bot files do not contain info at all, to exclude later
        emptybotfile = cat(1,emptybotfile,i);
        continue
    end
    md = find(botdep <= max(botdepth(i,ii)),1,'last');
    
    if length(ii)>1
        rbotdepth(i,1:md) = botdep(1:md);
        rbotpres(i,1:md) = interp1(botdepth(i,ii),botpres(i,ii),botdep(1:md));
        rbotPT(i,1:md) = interp1(botdepth(i,ii),botPT(i,ii),botdep(1:md));
        rbotSP(i,1:md) = interp1(botdepth(i,ii),botSP(i,ii),botdep(1:md));
        rbotoxyg(i,1:md) = interp1(botdepth(i,ii),botoxyg(i,ii),botdep(1:md));
        rnitr(i,1:md) = interp1(botdepth(i,ii),nitr(i,ii),botdep(1:md));
        rphos(i,1:md) = interp1(botdepth(i,ii),phos(i,ii),botdep(1:md));
        rsili(i,1:md) = interp1(botdepth(i,ii),sili(i,ii),botdep(1:md));
    else
        % if only a single value -- nearest neighbor the bottle data to nearest
        % depth interval (between 0-12m mostly)
        nnd = knnsearch(botdep',botdepth(i,ii));
        rbotdepth(i,nnd) = botdep(nnd);
        rbotpres(i,nnd) = botpres(i,ii);
        rbotPT(i,nnd) = botPT(i,ii);
        rbotSP(i,nnd) = botSP(i,ii);
        rbotoxyg(i,nnd) = botoxyg(i,ii);
        rnitr(i,nnd) = nitr(i,ii);
        rphos(i,nnd) = phos(i,ii);
        rsili(i,nnd) = sili(i,ii);
        c = c +1;
    end
end

BOTkeys(emptybotfile) = [];
botnettowdep(emptybotfile) = [];
rbotdepth(emptybotfile,:) = [];
rbotpres(emptybotfile,:) = [];
rbotPT(emptybotfile,:) = [];
rbotSP(emptybotfile,:) = [];
rbotoxyg(emptybotfile,:) = [];
rnitr(emptybotfile,:) = [];
rphos(emptybotfile,:) = [];
rsili(emptybotfile,:) = [];

%% Extract the most common bottle depths
% Intervals are based on inspection of the most commonly sampled depths
cbd = [5:5:25,30,40,50:25:200,250,300,400,500,600,800,1000,1250,1500:500:4000];

% A = rnitr(:,cbd);

%% Integrate across surface (0-50), MLD (0.03 threshold), and net tow depth

botnetmeta = [];
chemtable = [];
for i = 1:length(BOTkeys)
    netd = botnettowdep(i);
    
    % net tow depth and coverage for N, P, Si
    botnetmeta(i,:) = [netd, sum( ~isnan( rnitr(i,1:netd) )),... 
        sum( ~isnan( rphos(i,1:netd) )), sum( ~isnan( rsili(i,1:netd) ))];
    
    chemtable(i,:) = [nanmean(rnitr(i,1:50)),nanmean(rnitr(i,1:netd)),...
        nanmean(rphos(i,1:50)),nanmean(rphos(i,1:netd)),...
        nanmean(rsili(i,1:50)),nanmean(rsili(i,1:netd)),...
        nanmean(rbotoxyg(i,1:50)),nanmean(rbotoxyg(i,1:netd))];
end



%% Create datatable of CTD values
% [netkeys, coverage ...]

% Merge coverage and chem table
bottable = cat(2,botnetmeta,chemtable);

% there are rows with no integrated bottle data, maybe the bottle depths
% were deeper than the net tow depth
rs = find( nansum(bottable,2) > 0 );
bottable = cat(2, BOTkeys(rs), num2cell(bottable(rs,:)) );


botintvars = {'NetKey','NetTowDepth',...
    'CoverageN','CoverageP','CoverageSi',...
    'NitrI50','NitrINet',...
    'PhosI50','PhosINet',...
    'SiliI50','SiliINet',...
    'BotOxygI50','BotOxygINet'};

T3 = array2table(bottable,'VariableNames',botintvars);
writetable(T3,strcat(fol,'bottlenet_integrated_mld_10192021.xlsx'));


%% Save the workspace for plotting

save('ctd_bottle_profiles_integrated_10192021.mat')

%% Plot average depth profiles
% load ctd_bottle_profiles_integrated_10192021

% load bioregionalization and list of keys
% bioreg = readtable('Bioregion_phychm_used.csv');

% Keys to extract
A = bioreg.Key;

% Find which CTD keys match the net keys (kb)
B = keylist(:,[1:2]);
[kk,~, kb] = intersect( A, B(:,1)) ;

regions = {'Offshore','DeepShelf','Nearshore','DeepFjord'};
k4SA = nan(4,size(SA,2),4);
k4CT = nan(4,size(SA,2),4);
k4oxyg = nan(4,size(SA,2),4);
k4nitr = nan(4,size(rnitr,2),4);
k4phos = nan(4,size(rnitr,2),4);
k4sili = nan(4,size(rnitr,2),4);

k4deplimsts = nan(4,6);
k4deplimoxy= nan(4,6);
k4deplimnut= nan(4,6);
k4ctdnut = nan(4,3); %[T&S,O,N&P&Si]

% *** For T and S ***
figure;
for k = 1:4
    m = find(contains(bioreg.Bioregion, regions{k}));
    [nn,kc,~] = intersect(CTDkeys, keylist(kb(m),2)); % removes the repeated CTD files
    
    k4ctdnut(k,1) = length(kc);
    
    % calculate mean CTD depth
    zdep = rdepth(kc,:);
    deplims = max(zdep');
    deplims = [mean(deplims), std(deplims), min(deplims), max(deplims),...
        quantile(deplims,[0.01,0.99])];
    
    % extract the 1m binned T and S data
    % [mean, std, min, max]
    SA_mean = [nanmean(SA(kc,:)); nanstd(SA(kc,:)); nanmin(SA(kc,:)); nanmax(SA(kc,:))]';
    CT_mean = [nanmean(CT(kc,:)); nanstd(CT(kc,:)); nanmin(CT(kc,:)); nanmax(CT(kc,:))]';
    
    % aggregate information
    k4SA(k,:,:) = SA_mean;
    k4CT(k,:,:) = CT_mean;
    k4deplimsts(k,:) = deplims;
    
    % Average depth profiles
    deps = 1:400;
    subplot(4,3,1+(k-1)*3)
    plot(CT_mean(deps,1),deps,'r-','linewidth',2);
    hold on;
    plot([CT_mean(deps,1)-CT_mean(deps,2),CT_mean(deps,1)+CT_mean(deps,2)],deps,...
        'r--','linewidth',2);
    plot([CT_mean(deps,3),CT_mean(deps,4)],deps,'r:');
    set(gca, 'YDir','reverse')
    xlabel('Conservative temperature');
    ylabel('Depth (m)');
    xlim([2,20]);
    
    % Mark mean depth +- std (alternatively, can be vertical line)
    plot([3,3],[deplims(1),deplims(1)-deplims(2)],'k--+');
    plot([3,3],[deplims(1),deplims(1)+deplims(2)],'k--+');
    scatter(3, deplims(1),'ko','filled');
    hline(deplims(1),'k:');
    ylim([1,400]);
    
    subplot(4,3,2+(k-1)*3)
    plot(SA_mean(deps,1),deps,'b-','linewidth',2);
    hold on;
    plot([SA_mean(deps,1)-SA_mean(deps,2),SA_mean(deps,1)+SA_mean(deps,2)],deps,...
        'b--','linewidth',2);
    plot([SA_mean(deps,3),SA_mean(deps,4)],deps,'b:');
    set(gca, 'YDir','reverse')
    xlabel('Absolute Salinity');
    xlim([15,37]);
    hline(deplims(1),'k:');
    ylim([1,400]);
    
    
   
end

% *** For oxygen ***
for k = 1:4
   odata = cell2num(bioreg.Oxygen_surf);
    m = find(contains(bioreg.Bioregion, regions{k}) & ~isnan(odata));
    [nn,kc,~] = intersect(CTDkeys, keylist(kb(m),2)); % removes the repeated CTD files
    
    k4ctdnut(k,2) = length(kc);
    
    % calculate mean CTD depth
    zdep = rdepth(kc,:);
    deplims = max(zdep');
    deplims = [mean(deplims), std(deplims), min(deplims), max(deplims),...
        quantile(deplims,[0.01,0.99])]; 
    
    % extract the 1m binned data
    % [mean, std, min, max]
    oxyg_mean = [nanmean(roxyg(kc,:)); nanstd(roxyg(kc,:)); nanmin(roxyg(kc,:)); nanmax(roxyg(kc,:))]';
    % aggregate information
    k4oxyg(k,:,:) = oxyg_mean;
    k4deplimoxy(k,:) = deplims;
    
    
   subplot(4,3,3+(k-1)*3) 
   plot(oxyg_mean(deps,1),deps,'m-','linewidth',2);
   hold on;
   plot([oxyg_mean(deps,1)-oxyg_mean(deps,2),oxyg_mean(deps,1)+oxyg_mean(deps,2)],deps,...
       'm--','linewidth',2);
   plot([oxyg_mean(deps,3),oxyg_mean(deps,4)],deps,'m:');
   set(gca, 'YDir','reverse')
   xlabel('Dissolved Oxygen (umol/kg) 2001-2014');
   xlim([0 540]);
   hline(deplims(1),'k:');
   ylim([1,400]);
end

% ** For Nutrients **
for k = 1:4
    odata = cell2num(bioreg.Nitrate_surf);
    m = find(contains(bioreg.Bioregion, regions{k}) & ~isnan(odata));
    [nn,kc,~] = intersect(BOTkeys, keylist(kb(m),3)); % removes the repeated CTD files
    
    k4ctdnut(k,3) = length(kc);
    
    % calculate mean bottle tow depth
    zdep = rbotdepth(kc,:);
    deplims = max(zdep');
    deplims = [mean(deplims), std(deplims), min(deplims), max(deplims),...
        quantile(deplims,[0.01,0.99])]; 

    % extract the 1m binned data
    % [mean, std, min, max]
    k4nitr(k,:,:) = [nanmean(rnitr(kc,:)); nanstd(rnitr(kc,:)); nanmin(rnitr(kc,:)); nanmax(rnitr(kc,:))]';
    k4phos(k,:,:) = [nanmean(rphos(kc,:)); nanstd(rphos(kc,:)); nanmin(rphos(kc,:)); nanmax(rphos(kc,:))]';
    k4sili(k,:,:) = [nanmean(rsili(kc,:)); nanstd(rsili(kc,:)); nanmin(rsili(kc,:)); nanmax(rsili(kc,:))]';
    k4deplimnut(k,:) = deplims;
end

   
figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:4
   subplot(4,3,1+(k-1)*3)
   plot(k4nitr(k,:,1),botdep,'b-','linewidth',2);
   hold on;
   plot([k4nitr(k,:,1)-k4nitr(k,:,2);k4nitr(k,:,1)+k4nitr(k,:,2)],botdep,...
       'b--','linewidth',2);
   plot([k4nitr(k,:,4);k4nitr(k,:,4)],botdep,'b:');
   set(gca, 'YDir','reverse')
   xlabel('Nitrate');
   ylabel('Depth (m)');
   xlim([0, 50]);
   ylim([0, 400]);
   
   subplot(4,3,2+(k-1)*3)
   plot(k4phos(k,:,1),botdep,'r-','linewidth',2);
   hold on;
   plot([k4phos(k,:,1)-k4phos(k,:,2);k4phos(k,:,1)+k4phos(k,:,2)],botdep,...
       'r--','linewidth',2);
   plot([k4phos(k,:,4);k4phos(k,:,4)],botdep,'r:');
   set(gca, 'YDir','reverse')
   xlabel('Phosphate');
   ylabel('Depth (m)');
   xlim([0, 6]); % there are some outliers beyond 6
   ylim([0, 400]);
   
   subplot(4,3,3+(k-1)*3)
   plot(k4sili(k,:,1),botdep,'m-','linewidth',2);
   hold on;
   plot([k4sili(k,:,1)-k4sili(k,:,2);k4sili(k,:,1)+k4sili(k,:,2)],botdep,...
       'm--','linewidth',2);
   plot([k4sili(k,:,4);k4sili(k,:,4)],botdep,'m:');
   set(gca, 'YDir','reverse')
   xlabel('Silicate');
   ylabel('Depth (m)');
   xlim([0, 125])
   ylim([0, 400]);
end
%% Figure with all depth profiles
clrs4 = hex2rgb({'#004DFB','#FB16CA','#FFA500','#16FF00'});

% Combined average profiles across clusters - version 2, dashed beyond tow depth
% For nutrients, start at 5m depth
deps = 1:400;

figure('units','normalized','outerposition',[0 0 0.5 1]);
subplot(3,2,1);
hold on;
for k = 1:4
   mcd = round(k4deplimsts(k,1)); % mean ctd depth
   if mcd > 400
       mcd = 400;
   end
   plot(squeeze(k4CT(k,1:mcd,1)),deps(1:mcd),'-','linewidth',2,'color',clrs4(k,:));
end
set(gca, 'YDir','reverse')
xlabel('Temperature (°C)');
ylabel('Depth (m)');
xlim([4,16]);
ylim([0,400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

subplot(3,2,2);
hold on;
for k = 1:4
   mcd = round(k4deplimsts(k,1)); % mean ctd depth
   if mcd > 400
       mcd = 400;
   end
   plot(squeeze(k4SA(k,1:mcd,1)),deps(1:mcd),'-','linewidth',2,'color',clrs4(k,:));

end
set(gca, 'YDir','reverse')
xlabel('Salinity');
ylabel('Depth (m)');
xlim([20 35]);
ylim([0,400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

subplot(3,2,3);
hold on;
for k = 1:4
   mcd = round(k4deplimoxy(k,1)); % mean ctd depth
   if mcd > 400
       mcd = 400;
   end
   plot(squeeze(k4oxyg(k,1:mcd,1)),deps(1:mcd),'-','linewidth',2,'color',clrs4(k,:));
end
set(gca, 'YDir','reverse')
xlabel('Dissolved Oxygen (umol/kg)');
ylabel('Depth (m)');
xlim([0 350]);
ylim([0,400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

subplot(3,2,4);
hold on;
for k = 1:4
   mcd = knnsearch(botdep', round(k4deplimnut(k,1)) ); % get nearest bottle depth
   plot(squeeze(k4nitr(k,5:mcd,1)),botdep(5:mcd),'-','linewidth',2,'color',clrs4(k,:));
end
set(gca, 'YDir','reverse')
xlabel('Nitrate (umol/kg)');
ylabel('Depth (m)');
xlim([0, 50]);
ylim([0, 400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

subplot(3,2,5);
hold on;
for k = 1:4
   mcd = knnsearch(botdep', round(k4deplimnut(k,1)) ); % get nearest bottle depth
   plot(squeeze(k4phos(k,5:mcd,1)),botdep(5:mcd),'-','linewidth',2,'color',clrs4(k,:));
end
set(gca, 'YDir','reverse')
xlabel('Phosphate (umol/kg)');
ylabel('Depth (m)');
xlim([0, 3.2]); % there are some outliers beyond 6
ylim([0, 400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

subplot(3,2,6);
hold on;
for k = 1:4
   mcd = knnsearch(botdep', round(k4deplimnut(k,1)) ); % get nearest bottle depth
   plot(squeeze(k4sili(k,5:mcd,1)),botdep(5:mcd),'-','linewidth',2,'color',clrs4(k,:));
end
set(gca, 'YDir','reverse')
xlabel('Silicate (umol/kg)');
ylabel('Depth (m)');
xlim([0, 85])
ylim([0, 400]);
box on;
set(gca,'FontSize',14,'FontName','Arial');

% export_fig('supp_fig_waterprop_profiles.png','-transparent','-dpng','-r600');
%% Chla fitted splines
clrs4 = hex2rgb({'#004DFB','#FB16CA','#FFA500','#16FF00'});
regions = {'Offshore','DeepShelf','Nearshore','DeepFjord'};

netMonth = bioreg.Month;
chla = bioreg.Chlorophyll;

adjR2 = [];
rmse = [];

figure('units','normalized','outerposition',[0 0 0.5 0.4]);
hold on;
for i = 1:4
    ii = find(contains(bioreg.Bioregion, regions{i}));
    
    chlnow = chla(ii);
    jj = find(~isnan(chlnow));
    dmn = netMonth(ii);

    % fit a curve
    [fitobject,gof] = fit(dmn(jj),chlnow(jj),'smoothingspline');
    adjR2(i) = gof.adjrsquare;
    rmse(i) = gof.rmse;
    
    % confidence interval - need other fitting function for confidence
    % intervals ... or R
%     ci = confint(fitobject);
%     ci = csapsGCV(); ??
    
    if i == 1
        hL = plot(fitobject,netMonth,mat.Chla);
        set(hL,'color',[0.5 0.5 0.5],'linewidth',1)
        set(hL,'size',0);
        hL = plot(fitobject);
        set(hL,'color',clrs4(i,:),'linewidth',3)
    else
        hL= plot(fitobject);
        set(hL,'color',clrs4(i,:),'linewidth',3)
    end
end
% title('Satellite chlorophyll-a fitted splines','FontSize',14)
ylim([0 6]);
box on;
set(gca,'xtick',4:10,'xticklabel',{'Apr','May','Jun','Jul','Aug','Sep','Oct'},'FontSize',14);
ylabel('Chl-a (mgm^-^3)','FontSize',14)
xlabel('Month','FontSize',14)
xlim([3.75 10.25])

% export_fig('supp_fig_chlaprof.png','-transparent','-dpng','-r600');