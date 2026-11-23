%%MODIFED Code by Munir Gunes Kutlu and Banu Kutlu for Calipari Lab
% Contact gunes.kutlu@vanderbilt.edu
%%MODIFED Code by Munir Gunes Kutlu and Banu Kutlu for Calipari Lab
% Contact gunes.kutlu@vanderbilt.edu
% 
% clear all; close all;
% 
% low_speed = [];
% mid_speed = [];
% high_speed = [];
% full_trace_zscore = [];

clearvars; close all; %Clear all variables and close all figures
BlockDir=uigetdir('/Users/mariedoyle/Documents/Vanderbilt/Winder Lab/Operant behavior/Male and female c57 in STAR/Marie_Med-TDT files from GCaMP operant 030722/Organized files/First EtOH binge/1_1_032022'); %Get folder name
cd(BlockDir); %Change directory to photometry folder
[Tank,Block,~]=fileparts(cd); %List full directory for Tank and Block
data=TDTbin2mat(BlockDir); %Use TDT2Mat to extract data.


%% SYNAPSE- Extract Relevant Data from Data file
%Create Variables for each Photometry Channel and timestamps
Ch470=data.streams.x470A.data; %GCaMP
Ch405=data.streams.x405A.data; %Isosbestic Control
Ts = ((1:numel(data.streams.x470A.data(1,:))) / data.streams.x470A.fs)'; % Get Ts for samples based on Fs
StartTime=1; %Set the starting sample(recommend eliminating a few seconds for photoreceiver/LED rise time).
EndTime=length(Ch470); %Set the ending sample (again, eliminate some).
Fs=data.streams.x470A.fs; %Variable for Fs
Ts=Ts(StartTime:EndTime); % eliminate timestamps before starting sample and after ending.
Ch470=Ch470(StartTime:EndTime);
Ch405=Ch405(StartTime:EndTime);

%% Function to get DeltaF/F using regression fit
%Function to get DF/F for whole session
% Delta470=data.streams.x470A.data;
% Delta470=DeltaF(Ch470,Ch405);
% norm470 = smooth((Ch470./Ch405),0.0004,'lowess');
%norm470 = (Ch470-Ch405).^2;
%norm470 = Ch405.*Ch405+Ch470;
norm470 = smooth(Ch470,0.0004,'lowess');
Delta470 = (norm470-mean(norm470))/mean(norm470);

%% Taking TTL times and Delta470 values around TTLs
%data.epocs.Po6_.onset - Tone
data_to_use = data.epocs.LLr_.onset(2:end-1);
TTL_size = numel(data_to_use);
TTL_times = zeros(TTL_size:1);
TTL_temp = zeros(TTL_size:1);
interval_pre = -8000; %milliseconds
interval_post = 8000; %milliseconds
interval_count = abs(interval_pre) + abs(interval_post);
TTL_signal = zeros(interval_count, TTL_size);
for TTL_index = 1:TTL_size %% timestamp
    [c, ind] = min(abs(Ts-data_to_use(TTL_index, 1))); %returns the position best fit between T and Po0 onset
    
    position = 0;
    for interval_ind = interval_pre:interval_post
        position = position+1;
        TTL_signal(position, TTL_index) = Delta470(ind + interval_ind);
    end
    
    TTL_temp(TTL_index) = Delta470(ind); %writes the Delta490 value for the best T+Po0 position
    TTL_times(TTL_index) = ind;
end



%% Calculate z-score
% baseline_total = (TTL_times(1)-660000):(TTL_times(1)-420000); %time you want to use as baseline for the total trace z score calculation 
% b_mean = mean(Delta470(baseline_total)); % baseline period mean
% b_stdev = std(Delta470(baseline_total)); % baseline period stdev
% z_full_trace = (Delta470 - b_mean)/b_stdev;
% z_full_trace = z_full_trace((TTL_times(1)-660000):(TTL_times(1)+660000));

baseline_pre = -8000; %milliseconds
baseline_post = 0; %milliseconds
baseline_size = abs(baseline_pre) + abs(baseline_post);
baseline = zeros(size(TTL_signal,1));
z_all = zeros(size(TTL_signal));


for interval_index = 1:TTL_size
    baseline = TTL_signal(1:baseline_size, interval_index);
    b_mean = mean(baseline); % baseline period mean
    b_stdev = std(baseline); % baseline period stdev
    for TTL_index = 1:size(TTL_signal(:, interval_index)) % Z score per bin
        z_all(TTL_index, interval_index) = (TTL_signal(TTL_index, interval_index) - b_mean)/b_stdev;
    end
end



%% 
[bout_start_timestamps, bout_end_timestamps]=bout_analysis(TTL_times'/1000);


TTL_size = numel(bout_start_timestamps);
TTL_times = zeros(TTL_size:1);
TTL_temp = zeros(TTL_size:1);
interval_pre = -8000; %milliseconds
interval_post = 16000; %milliseconds
interval_count = abs(interval_pre) + abs(interval_post);
TTL_signal = zeros(interval_count, TTL_size);
for TTL_index = 1:TTL_size %% timestamp
    [c, ind] = min(abs(Ts-bout_start_timestamps(TTL_index, 1))); %returns the position best fit between T and Po0 onset
    
    position = 0;
    for interval_ind = interval_pre:interval_post
        position = position+1;
        TTL_signal(position, TTL_index) = Delta470(ind + interval_ind);
    end
    
    TTL_temp(TTL_index) = Delta470(ind); %writes the Delta490 value for the best T+Po0 position
    TTL_times(TTL_index) = ind;
end

baseline_pre = -8000; %milliseconds
baseline_post = 0; %milliseconds
baseline_size = abs(baseline_pre) + abs(baseline_post);
baseline = zeros(size(TTL_signal,1));
z_all_bouts = zeros(size(TTL_signal));


for interval_index = 1:TTL_size
    baseline = TTL_signal(1:baseline_size, interval_index);
    b_mean = mean(baseline); % baseline period mean
    b_stdev = std(baseline); % baseline period stdev
    for TTL_index = 1:size(TTL_signal(:, interval_index)) % Z score per bin
        z_all_bouts(TTL_index, interval_index) = (TTL_signal(TTL_index, interval_index) - b_mean)/b_stdev;
    end
end
%% 

% plot(z_all)
% tracesize(jj) = size(z_full_trace,1);
% 
% low_speed = [low_speed,z_all(:,1)];
% mid_speed = [mid_speed,z_all(:,2)];
% high_speed = [high_speed,z_all(:,3)];
% full_trace_zscore = [full_trace_zscore,z_full_trace(1:tracesize(1))];

%Save TTL_signal to CSV file
% csvwrite('Zhaker_events.csv',TTL_signal);
% csvwrite('peaks(Zscore).csv',z_all);


%plot(TTL_signal,'DisplayName','TTL_signal')
%hold on


% figure
% plot(Delta470)
% hold
% x=(1:length(Delta470));
% z = zeros(length(x),1) + -.01;
% p= plot(x,z,'-^','MarkerIndices',TTL_times)
% p.MarkerSize= 10
% p.MarkerFaceColor= 'green'


% figure 
% plot(mean(low_speed,2))
% hold on
% plot(mean(mid_speed,2))
% plot(mean(high_speed,2))
% 
% remainder = rem(size(low_speed,1),1000); % truncate the data to have equal number of datapoints in each bin
% low_speed_binned = bin_data(1000,low_speed(1:end-remainder,:),@mean);
% mid_speed_binned = bin_data(1000,mid_speed(1:end-remainder,:),@mean);
% high_speed_binned = bin_data(1000,high_speed(1:end-remainder,:),@mean);
% 
% remainder = rem(size(full_trace_zscore,1),1000); % truncate the data to have equal number of datapoints in each bin
% full_trace_zscore_binned = bin_data(1000,full_trace_zscore(1:end-remainder,:),@mean);
% 
% 
% csvwrite('low_speed.csv',low_speed);
% csvwrite('mid_speed.csv',mid_speed);
% csvwrite('high_speed.csv',high_speed);
% csvwrite('full_trace_zscore_binned.csv',full_trace_zscore_binned);


