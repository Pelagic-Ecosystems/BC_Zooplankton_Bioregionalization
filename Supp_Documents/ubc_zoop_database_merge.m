% UBC Zoop Database Merge
% August 1, 2019
% Revised April 21, 2020 to receive data files without comments
% Also added columns for day,month,year

% fol = 'C:\Users\patri\Desktop\DFO_Vertical_Nets\';
fol = 'D:\UBC Files\Sync\Zooplankton-data-files\orig_files_no_comments\';
files = ls(fol);


i_wcbc = [3,6,7,8,10,12,15];
i_sog = [4,9,11,14,16];
i_main = [5];

NUM_WCBC = [];
TXT_WCBC = [];
RAW_WCBC = [];
for i = 1:length(i_wcbc)
    f = strcat(fol,files(i_wcbc(i),:));
    disp(f);
    [NUM,TXT,RAW] = xlsread(f);
    
    NUM_WCBC = cat(1,NUM_WCBC,NUM);
    TXT_WCBC = cat(1,TXT_WCBC,TXT(2:end,:));
    RAW_WCBC = cat(1,RAW_WCBC,RAW(2:end,:));
    
end


NUM_SOG = []; TXT_SOG = []; RAW_SOG = [];
for i = 1:length(i_sog)
    f = strcat(fol,files(i_sog(i),:));
    disp(f);
    [NUM,TXT,RAW] = xlsread(f);
    
    NUM_SOG = cat(1,NUM_SOG,NUM);
    TXT_SOG = cat(1,TXT_SOG,TXT(2:end,:));
    RAW_SOG = cat(1,RAW_SOG,RAW(2:end,:));
end

f = strcat(fol,files(i_main(1),:));
disp(f);
[NUM_MAIN,TXT_MAIN,RAW_MAIN] = xlsread(f);

f = strcat(fol,files(17,:));
disp(f);
[NUM_GOA,TXT_GOA,RAW_GOA] = xlsread(f,'Sheet2');
TXT_GOA = TXT_GOA(2:end,:);
RAW_GOA = RAW_GOA(2:end,:);

% The line P data is already in the above files
% f = strcat(fol,files(18,:));
% disp(f);
% [NUM_P,TXT_P,RAW_P] = xlsread(f,'Sheet2');
% TXT_P = TXT_P(2:end,:);
% RAW_P = RAW_P(2:end,:);

%
% Concatenate all... MAIN has the headers
RAW_BC = cat(1,RAW_MAIN,RAW_SOG);
RAW_BC = cat(1,RAW_BC,RAW_WCBC);
RAW_BC = cat(1,RAW_BC,RAW_GOA);

TXT_BC = cat(1,TXT_MAIN,TXT_SOG);
TXT_BC = cat(1,TXT_BC,TXT_WCBC);
TXT_BC = cat(1,TXT_BC,TXT_GOA);

NUM_BC = cat(1,NUM_MAIN,NUM_SOG);
NUM_BC = cat(1,NUM_BC,NUM_WCBC);
NUM_BC = cat(1,NUM_BC,NUM_GOA);

% assign proper variable type to numerical values
NUM_BC2 = cell2num(TXT_BC(2:end,5:6)); % lon and lat
NUM_BC2(:,3) = datenum(TXT_BC(2:end,7),'ddmmyyyy'); % datenum
NUM_BC2(:,4) = NUM_BC(:,1); % time
NUM_BC2(:,5:6) = NUM_BC(:,4:5); % mesh size, net mouth diameter
NUM_BC2(:,7:10) = cell2num(TXT_BC(2:end,13:16)); % depth start, depth end,bottom depth, volume filtered
NUM_BC2(:,11) = NUM_BC(:,10); %ctd key
NUM_BC2(:,12:16) = NUM_BC(:,17:21); %abundance, biomass, num sp, stn div, stn equitability
NUM_HEADER = TXT(1,[5:8,11:17,24:28]);

NUM_BC2(:,17:19) = cell2num( split(TXT_BC(2:end,7)) ); % day month year
NUM_HEADER = [NUM_HEADER,{'Day','Month','Year'}];

% crop out text variables
TXT_BC2 = TXT_BC(:,[1:4,9:10,18:23]);
% separate genus, species, and stage

genus = {};
species = {};
stage = {};
for i = 2:length(TXT_BC2)
    m = split(TXT_BC2(i,12));
    genus(i-1) = m(1);
    species(i-1) = join(m(1:2));
    if length(m) > 2
        stage(i-1) = join(m(3:end));
    else
        stage(i-1) = '';
    end
end
TXT_BC2(1,13:15) = {'Genus','Species','Stage'};
TXT_BC2(2:end,13) = genus;
TXT_BC2(2:end,14) = species;
TXT_BC2(2:end,15) = stage;

%% Update Acartia info 
% August 8, 2020
% look for acartia samples
% extract all a tonsa and californiensis
% switch out
TXT_BC3 = TXT_BC2;
NUM_BC3 = NUM_BC2;

file = 'D:\UBC Files\Sync\Zooplankton-data-files\acartia_adjusted_tonsa_californiensis_07082020.xlsx';
[acartiaNUM,acartiaTXT,acartiaRAW] = xlsread(file);

keys = unique(acartiaRAW(2:end,1));
replace_count = 0;
unevenswitch = 0;
for i = 1:length(keys)
    jj = find(strcmp(acartiaRAW(:,1),keys{i}));
    spnew = strcat(acartiaRAW(jj,31),{' '},acartiaRAW(jj,32),{'  '},...
        acartiaRAW(jj,35));
    abubio = cell2num(acartiaRAW(jj,41:42));
    
    key = keys{i};
    key = key([1:7,9:end]);
    ii = find(strcmp(TXT_BC3(:,1), key) & ...
        contains(TXT_BC3(:,14),{'Acartia tonsa','Acartia californiensis'}));
    
    if ~isempty(ii)
       % remove old acartia 
        oldtaxa = TXT_BC3(ii,12);
        
        if ~isequal(spnew,oldtaxa) && (length(spnew) == length(oldtaxa) )
             TXT_BC3(ii,12) = spnew;
             TXT_BC3(ii,14) = strcat(acartiaRAW(jj,31),{' '},acartiaRAW(jj,32));
             TXT_BC3(ii,15) = acartiaRAW(jj,35);
             NUM_BC3(ii-1,12:13) = abubio;
        end
        if (length(spnew) ~= length(oldtaxa) )
            disp('more organisms than old version');
            break;
        end
       
    end
end


% Confirm by checking if total abundance is the same... total biomass can
% change a little because of small difference in A.tonsa and A.californiensis 
abu_change = sum(NUM_BC2(:,12)) - sum(NUM_BC3(:,12));
bio_change = sum(NUM_BC2(:,13)) - sum(NUM_BC3(:,13));
NUM_BC2 = NUM_BC3;
TXT_BC2 = TXT_BC3;
%% Save files
clearvars -except TXT_BC2 NUM_BC2 NUM_HEADER RAW_BC

outfol = 'D:\UBC Files\Sync\Zooplankton-data-files\';
save(strcat(outfol,'Zoop_BC_08122020.mat'),'TXT_BC2','NUM_BC2','NUM_HEADER');
% save(strcat(outfol,'Zoop_BC_RawSheet_08052019.mat'),'RAW_BC','-v7.3');
% takes 198 seconds

csvwrite_with_headers(strcat(outfol,'Zoop_BC_08122020_num.csv'),NUM_BC2,NUM_HEADER);

fileID = fopen(strcat(outfol,'Zoop_BC_08122020_text.csv'),'w+');
formatSpec = '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n';
[nrows,ncols] = size(TXT_BC2);
for row = 1:nrows
    fprintf(fileID,formatSpec,TXT_BC2{row,:});
end
fclose(fileID);
%% Print all unique keys with comments
% load('D:\PhD_Files\Data\DFO\DFO_Vertical_Nets\Zoop_BC_RawSheet_08052019.mat');

% % station subset
[C,IA,IC] = unique(RAW_BC(2:end,1));
BC_Stns = RAW_BC(IA+1,:);
BC_Stns_header = RAW_BC(1,:);

outfol = 'D:\PhD_Files\Data\DFO\DFO_Vertical_Nets\';
save(strcat(outfol,'BC_stations_first_instance_08062019.mat'),'C','IA','IC','BC_Stns','BC_Stns_header');

% load ('D:\PhD_Files\Data\DFO\DFO_Vertical_Nets\BC_stations_first_instance_08062019.mat');

notes = BC_Stns(:,18);
ii = find(~strcmp(notes,' '));
% comment list
cl = cat(2,BC_Stns(ii,[1,2,3,7]),notes(ii));
%%
flaglist = {'flow','shell','phyto','poor','fraction','paint','algae','coscino','rbr',...
    'bottom','sand','benthic','volume','expected','split','ctd','cod end','angle',...
    'speed','remove','took','release'};
ff = find(contains(notes,flaglist(21),'IgnoreCase',true ));
z = notes(ff);

%% Output unique stations
[C2,IA2,IC2] = unique(TXT_BC2(2:end,1));

sub_stns = [TXT_BC2(1,[1:4,8]),NUM_HEADER([1:4,11])];
sub_stns = cat(1,sub_stns,[TXT_BC2(IA2+1,[1:4,8]),num2cell(NUM_BC2(IA2,[1:4,11]))]);

save('D:\PhD_Files\Data\DFO\DFO_Vertical_Nets\Zoop_event_list.mat',sub_stns);
