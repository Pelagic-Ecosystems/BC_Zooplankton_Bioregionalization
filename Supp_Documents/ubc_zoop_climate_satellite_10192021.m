% Extract climate mode and satellite info for table of predictors
% October 19, 2021

load('netctd_matching_03302021.mat')

%% Calculate DayOfYear variable from date values 

yr = cell2mat( T.NetYear );
dt = cell2mat( T.NetDT );
doy = floor( dt - datenum(yr,0,0,0,0,0));

%% Get table of climate indices

% Add climate modes to envivariables
index_name = {'NPGO','PDO','ENSO','EMI','NPI'};
climate_i = [];
[climate_i(:,1),climate_d] = get_climate_index('npgo',1979,2018);
[climate_i(:,2),climate_d]  = get_climate_index('pdo',1979,2018);
[climate_i(:,3),climate_d]  = get_climate_index('nina34',1979,2018);
[climate_i(:,4),climate_d]  = get_climate_index('emi',1979,2018);
[climate_i(:,5),climate_d]  = get_climate_index('np',1979,2018);
ind_year = str2num( datestr( climate_d, 'yyyy' ) );
ind_month = str2num( datestr( climate_d, 'mm' ) );

env_climate = [];
for i = 1:9949
    cc = find( ind_year == cell2num(T.NetYear(i)) & ...
        ind_month == cell2num(T.NetMonth(i)) );
    env_climate(i,:) = climate_i(cc,:);
end

%% Load ACRI data - updated 03312021
% already extracted values for each net tow
load net_acri_09242020

%% Create new modified envi table for R - Updated 10182021
envitable = cat(2,doy,env_climate);
envitable = cat(2,envitable,netacri);
envitable = cat(2,netmetable(:,1), num2cell(envitable));

Tclimsat = array2table(envitable);
Tclimsat.Properties.VariableNames = [{'Key','DayOfYear'},index_name,acrihead];


writetable(Tclimsat,'C:\Sync\Zooplankton-data-files\net_climate_satellite_10192021.xlsx');


