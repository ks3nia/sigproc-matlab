%% Create DATA structure with PSTH and Zscore
% updated to work with concatenated outputs.
% function DATA = createDATAstruct(path0) 
% path to kilosort and phy outputs (should be raw data location too)
%Edited by KV 8/21/23 for PF optrode recordings
%Edited by KV in October 2023 to add waveform properties, check sorting
%quality, added scrolling psth/raster plots, add auto search for responders, some quantification of firing rate etc

clear; clc;
addpath(genpath('Z:\Shared\MATLAB scripts\Ksenia\in vivo analysis\Kilosort post-processing'));
% %for all. must have folder with all folders containing KS output, both
%.events files, acqTimes.mat, and 'trialOrder.mat' file that is a cell array in proper format
path0 = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231130_DCNtest4_LPnolaser';
datafolders = dir(path0); 
datafolders = datafolders(~contains({datafolders(:).name},{'.','..','Measures','DATA'})); % if it includes any of the {}, it will exclude them. Because all filenames include .mat, it keeps only folders w/o DATA or Measures in their name.
[DATA,DATAlab] = createDATAstruct_v2_KV(path0,datafolders);
save(fullfile(path0,['DATA_' datestr(now,'YYYYmmddhhMM')]),'DATA','DATAlab');

%% Plot All Rasters
%can either run this right after creating new datastruct above or load a
%previous one in the folder of interest

m.krnl_size=16;
krnl = gausswin(m.krnl_size)/sum(gausswin(m.krnl_size));

close all;
% for f = 1:size(DATA.th.mW10,1) %this is the total num of "good" units
%     figure;
% %     subplot(3,1,1); plot(DATA.th.mW10{f,5},DATA.th.mW10{f,6}); ylim([0 40]);... %this is the raster
% subplot(3,1,1); plot(DATA.th.mW10{f,5},DATA.th.mW10{f,6});... %this is the raster
% %         subplot(3,1,2);plot(DATA.th.mW10{f,3}, conv2(DATA.th.mW10{f,2}/m.bin_size,krnl','valid')); %this looks odd.. 
%             subplot(3,1,2);plot(DATA.th.mW10{f,3}, DATA.th.mW10{f,2}); %this looks odd.. 
%         subplot(3,1,3); plot(DATA.th.mW10{f,3},DATA.th.mW10{f,4});
% %         ylim([-2 2]);    % this is the z scored FR (smoothed bins,
% %         smoothed and z scored psth) original
% %     subplot(3,1,3); plot(DATA.th.mW10{f,1},DATA.th.mW10{f,4}); ylim([-2 2]);    % this is the z scored FR (smoothed bins, smoothed and z scored psth)
% end

%ok something going on with the psth calculations in createDATAstruct. and
%in psthAndBa used later.. 
%note: sorted this out with AD's help 10/30-10/31. 

for f=1:size(DATA.th.mW10,1)
    figure;
    subplot(3,1,1); plot(DATA.th.mW10{f,5},DATA.th.mW10{f,6}); ylim([0 40]);... %this is the raster
    subplot(3,1,2);
%     plot(conv(m.timeAxis,ones(10,1),'valid')/10,conv(DATA.th.mW10{f,2},ones(10,1),'valid')/10);title(f); %smoothed psth
%     subplot(3,1,3);
%     plot(conv(m.timeAxis,ones(10,1),'valid')/10,conv(DATA.th.mW10{f,9},ones(10,1),'valid')/10);title(f); %smoothed z scored psth
    plot(DATA.th.mW10{f,3},conv(DATA.th.mW10{f,2},ones(10,1),'valid')/10);title(f); %smoothed psth
    subplot(3,1,3);
    plot(DATA.th.mW10{f,3},conv(DATA.th.mW10{f,9},ones(10,1),'valid')/10);title(f); %smoothed z scored psth
end

% 
% group1=cat(1,DATA.th.mW10{:,2}); 
% FR=group1./m.bin_size;


% Notes on outputs:
% - psth can be plotted with plot(bins, psth);
% - rasters can be plotted with plot(rasterX, rasterY);
% - spikeCounts is nEvents x 1, where each entry is the number of spikes
% that occurred within the window around that event. 

%DATAstruct outputs layout:
% DATA.(reg).(stim){rC.(stim),1} = bins; %the time bins based on the window
% and bin size set in parameters (this is the raw output from the
% psthRasterandCounts function)
% DATA.(reg).(stim){rC.(stim),2} = psth; %the psth for all good units
% within that time window (this is the raw output from the
% psthRasterandCounts function)
% DATA.(reg).(stim){rC.(stim),3} = binsSM; % smoothed bins (and spikes/sec)
% DATA.(reg).(stim){rC.(stim),4} = ZpsthSM; %z scored and smoothed FR
% (spikes/sec)
% DATA.(reg).(stim){rC.(stim),5} = rasterX;
% DATA.(reg).(stim){rC.(stim),6} = rasterY;
%7 = event times that have been adjusted for the correct timestamp, should
%hopefully work with the rest of the github scripts? it does!
%DATA.(reg).(stim){rC.(stim),8} = psthSM; smoothed psth in spikes/sec
%DATA.(reg).(stim){rC.(stim),9} = Zpsth; z scored psth in spikes/sec *not
%working
%  DATA.(reg).(stim){rC.(stim),10} = FR;  this is the raw FR (psth in
%  spikes/sec) &working? 

% DATAlab.(reg).(stim){rC.(stim),1} = datafolders(f).name;
% DATAlab.(reg).(stim){rC.(stim),2} = spike_id;



%%KV note: 102423; some confusion still here about which variables to use
%%in the find responders function below. leave for now. 


%% Find responders: adapted from Collective_1v14_KVopto4
%1. go into folder  (ex: Z:\Shared\DATA\Electrophysiology\In vivo\KV\in
%vivo\PF optrode\z20230929_PF optrode\09292023_40trialsLP)
%2. load the most recent data and measures files for the recording of
%interest (can run this section independently) - note: for 1 second baseline load
%data and measures from 10/16.***********

%most useful to run this and then look at the psths/rasters of each cluster
%to see if it's picking out reasonable "responders". To do that, run
%section 3 in SUMMARY DATA ANALYSIS and look at NameRespI_clusters printed
%here to match

% KV edit 10/13/23: trying to figure out different windows (24 seconds, 20
% seconds, 2 seconds)

%KV edit 10/17/23: tried a different way to figure out the window
%(uncommented). since there are periodic fluctuations in the signal with ketamine (~1sec intervals), changed the z scoring in createDATAstruct_v2_KV (edited
%baseline epoch to be 2.3 seconds instead of 1) but that doesn't seem to make a big difference in terms of ease of finding responders (makes heatmap for 0928 animal better though). play around with response
%window width too. go back to original (1 second before baseline)?

clear all; close all;
%%%For ArchT's:
load('DATA_202310241701.mat'); %9/15
% load('DATA_202310241706.mat'); %9/20
% load('DATA_202310241705.mat'); %9/18
% load('DATA_202310241708.mat'); %9/29
% load('DATA_202310241709.mat'); %10/02

load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230912_PF optrode\09122023_40trialsLP\thMeasures202310241604.mat'); %REMEMBER TO LOAD THE RIGHT MEASURES FILE!

%%%for Controls:
% load('DATA_202310241631.mat'); %0913 gfp female
% load('DATA_202310241707.mat'); %0928 gfp male
% load('DATA_202310241604.mat'); %912 female no virus


Zcrit = 2; %was 3 for BLA
% ZcritInh = -0.5; %was -1 for the BLA recordings (played around with higher values too, gets rid of what looks like real responders if lower)
% respWinWidth = 0.1; % response window (in seconds);  %80 ms; 100 ms
ZcritInh = -0.75; 
% ZcritInh = -0.9; 
respWinWidth =5; % response window (in seconds);  % tried 80 ms; 100 ms
% respWinWidth = 10; % response window (in seconds);  
% respWinWidth = 20; % response window (in seconds);
respWin = find(m.timeAxis>0 & m.timeAxis<=respWinWidth); % response window e.g. first 10 seconds
respWinI = respWin;

group1=cat(1,DATA.th.mW10{:,4}); %z scored smoothed FR from data structure
% group1=cat(1,DATA.th.mW10{:,9}); %z scored FR from data structure
% %updated 1024 so this isnt the smoothed data... confusing 

% EA: Is DATA.th.mW10{:,4} a row vector or a column vector? To concat in 1st dimension it must be a row vector. If it's a column instead, you need to transpose: group1=cat(1,DATA.th.mW10{:,4}');

%original:
% %Rspnder.(stimArea).PLowS = LowS.(stimArea)(any(colZ_shell(LowS.(stimArea),respWin)>Zcrit,2));
Rspnder_totalEx=[];
Rspnder_totalInh=[];

respwintotal=[];
Rspnder_sum=[];

%trying to do consecutive timebins for inhibitory responders 
for a=1:size(group1,1)
    unit=group1(a,:);
    RspnderEx = any(unit(respWin)>Zcrit,2); %first instance of threshold crossing %yes/no does this unit cross threshold during this time?
%     RspnderInh = any(unit(respWinI)<ZcritInh & unit(respWinI+1)<ZcritInh & unit(respWinI+2)<ZcritInh,2); %classic three consecutive time bins
    %EA: By shifting response window, we examine whether consequtive samples fall below ZcritInh
    %EA: E.g., if sample 100 is less than ZcritInh and you then shift respwin to the right by 1 point, then the new sample 100 was originally sample 101. Then you shift again by one place and so on.
    %EA: In order for the original sample 100 to meet the inhibitory response criterion, it means that samples 101 and 102 are also less than ZcritInh

    %stuff I tried to see if there's inhibition for the full 24 seconds:
   % RspnderInh = any(unit(respWinI)<ZcritInh &unit(respWinI+respWinWidth)<ZcritInh,2) %if reswinwidth=12, then this gives you 24 seconds? no
%     RspnderInh = any(unit(respWinI)<ZcritInh &unit(respWinI+respWinI)<ZcritInh,2) % this gives you 24 seconds if respWinWidth=10
%EA: To identify response that lasts that long the above approach is not practical, because you'd have the create a very long any(...) with 24 ANDs
%EA: An easier alternative is to do diff(find(unit(respWinI)>=ZcritInh)). The diff will show the distance between samples that are greater than Zcrit. If this distance exceeds 23, it means you have 24 consecutive samples that are less than Zcrit.
%EA: You can try max(diff(find(unit(respWinI)>=ZcritInh))) to gather the max duration of inhibition for each cell, but its upper bound will be the respWinI, which would underestimate inhibitions that start late and outlast the response window.

%EA: Based on the above, to test inhibition over X consecutive bins, you can use the following:
% RspnderInh = any(diff(find(unit(respWinI)>=ZcritInh))>=X); %EA: Where X can be any integer number of samples up to the number of samples in respWinI. Making the response window very long would lead to inclusion of responses that start very late.  
%EA: You can try out different values of X, changing in big chunks, like X=3 (our old method), X=10, X=20. Once you find an X that shows less than 5% responders in GFP cells, then evaluate the % responders in ArchT cells. You can then adjust the value of X in finer steps, until you get optimal % responders in ArchT cells and less than 5% responders in the GFP cells.
%EA: In Prism you can compare the proportion of responders in GFP vs. ArchT using a Chi-Square test for proportions. 

max(diff(find(unit(respWinI)>=ZcritInh)))

RspnderInh = any(diff(find(unit(respWinI)>=ZcritInh))>=20); %works with10

% RspnderInh = any(diff(find(unit(respWinI)>=ZcritInh))>=15);



    %another dumb way to do this (tried for loop earlier but am not figuring it out):
    %works if respwin=1? this is 24 consecutive timebins but not =24
    %seconds in the timeaxis
%     RspnderInh = any(unit(respWinI)<ZcritInh & unit(respWinI+1)<ZcritInh & unit(respWinI+2)<ZcritInh & unit(respWinI+3)<ZcritInh & unit(respWinI+4)<ZcritInh & unit(respWinI+5)<ZcritInh...
%         & unit(respWinI+6)<ZcritInh & unit(respWinI+7)<ZcritInh & unit(respWinI+8)<ZcritInh & unit(respWinI+9)<ZcritInh & unit(respWinI+10)<ZcritInh & unit(respWinI+11)<ZcritInh...
%         & unit(respWinI+12)<ZcritInh & unit(respWinI+13)<ZcritInh & unit(respWinI+14)<ZcritInh & unit(respWinI+15)<ZcritInh & unit(respWinI+16)<ZcritInh & unit(respWinI+17)<ZcritInh...
%         & unit(respWinI+18)<ZcritInh & unit(respWinI+19)<ZcritInh & unit(respWinI+20)<ZcritInh & unit(respWinI+21)<ZcritInh & unit(respWinI+22)<ZcritInh & unit(respWinI+23)<ZcritInh,2);

   %RspnderInh = any(unit(respWinI)<ZcritInh & unit(respWinI+1)<ZcritInh &
   %unit(respWinI+2)<ZcritInh & unit(respWinI+3)<ZcritInh,2); %this gives
   %you 24 seconds if respwin=20? no

  % RspnderInh = any(unit(respWinI)<ZcritInh & unit(respWinI+1)<ZcritInh & unit(respWinI+2)<ZcritInh & unit(respWinI+3)<ZcritInh & unit(respWinI+4)<ZcritInh & unit(respWinI+5)<ZcritInh,2); %this gives you 24 seconds if respwin=4?

    Rspnder_totalEx=[Rspnder_totalEx,RspnderEx];
    Rspnder_totalInh=[Rspnder_totalInh,RspnderInh];
end


%get the ones that are 1 vs 0

NameRespI=find(Rspnder_totalInh==1) %this gives you the index in data array to find the cluster
for ii=1:length(NameRespI)
[NameRespI_clust]=DATAlab.th.mW10{(NameRespI(ii)),2} %name of the cluster to check if that looks reasonable in the psth/raster plots later (section 3)
end

NumRespI=numel(NameRespI)%num of responders
NameRespE=find(Rspnder_totalEx==1)
NumRespE=numel(NameRespE)%num of responders
NumTotal = length(Rspnder_totalInh) %num total units for this intensity grouping
figure; X=[NumRespI/NumTotal]
labels = {'Inhibitory Responders'};
pie(X)

A=[NumRespE NumRespI (NumTotal-(NumRespI+NumRespE))] %this gives the relevant numbers of units

p=pie(A)
labels = {'Excite Resp';'Inh Resp';'NonResp'};
% unitstring=sprintf('total num units %f',NumTotal)

% Another way to do this and save actual values (# of neurons):
figure;
value = compose(['Num units: %0.0f'], A);
percents = compose([newline '%.2f%%'], A/sum(A)*100);
pie(A, strcat(value, percents))


lgd = legend(labels);
% colormap([0.30,0.75,0.93; 0.49,0.18,0.56; 0.65,0.65,0.65; 1 0 0]);
colormap([0.30, 0.75, 0.93; 1 0 0; 0.65,0.65,0.65]);



%% Make Heatmap of all cells
%load('Z:\Shared\DATA\Electrophysiology\In vivo\AD\DualLaserExps\PF\GoodPF\coreMeasures202305251048.mat')
%load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\Sorted PF\thMeasures202308231702.mat')
load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\thMeasures202310161441.mat'); %don't necessarily need to change this, just need timeaxis to be correct
group1=cat(1,DATA.th.mW10{:,4}); %this loads the z scored psths for all units

HTMP = figure;

% subplot(2,1,1)
imagesc(m.timeAxis,1:size(group1,1),group1);axis xy;colorbar;
colormap(viridis);
xlabel('Time (s) from stimulus onset');
ylabel('Cells');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005]);
set(HTMP,'color','w');
% clim([-1, 4]);
% xlim([-0.5,0.75]);
% xlim([-100,300]);
title('10mW inhibition');
xline(0,'color','w'); % marking stim on/off
hold on;
xline(24,'color','w');

%% Some notes: 
%m.timeAxis = -m.preStimWin:m.bin_size:m.postStimWin-m.bin_size; % time to be referenced to stimulus onset
%m.bsln_s = find(m.timeAxis>-1.1 & m.timeAxis<-0.1); % timeAxis in samples from -1s to stimulus onset
%load measures file for the file of interest
%m.bsln_s instead of bsln_s
% bsln_s = find(timeA>-1.1 & timeA<-0.1); % samples of last ~1s before stimulus
respWinWidth = 5; % response window (in seconds); 500 ms 
respWin = find(m.timeAxis>0.005 & m.timeAxis<=respWinWidth); % response window    

%in the BLA code we found the max and min FR for responders using the
%collective PSTH structure (non z scored): from collective_1v14_KVopto4
% group2=cat(1,DATA.th.mW10{:,8});%I think this? each column is a unit over time (rows)

group2= cat(1,DATA.th.mW10{:,2}); %this is spike/sec within the timebin (per 10ms timebin)


%not sure if this is correct yet, don't think its calculating fr correctly
%based on how kilosort is doing it... 
for i=1:size(group2,1) %for the # of units
    unit=group2(i,:); %go through unit of interest
    [mxS(i), mxS_smp(i)] = max(unit(respWin),[],2); %mxS =max FR, mxS_smp=sample of the max FR - this computes the maximum firing rate of each unit within the response window 
    [bslS(i)] = mean(unit(m.bsln_s),2,'omitnan'); %baseline FR for each unit
    [minSI(i), minSI_smp(i)] = min(unit(respWin),[],2); %min FR
%     [mxS_tm(i)] = m.timeAxis(respWin(1)+mxS_smp(i))'; % peak latency in time (at which time the maximum FR happens) 
%     [minSI_tm(i)] = m.timeAxis(respWin(1)+minSI_smp(i))'; % peak latency in time (at which time the min FR happens)
end
    
[mean(mxS) ,std(mxS)/sqrt(length(mxS))] %this gives you the mean and SEM (spikes/sec)
[mean(minSI) ,std(minSI)/sqrt(length(minSI))] %this gives you the mean and SEM (spikes/sec)

%note: for latencies, need to run the data code at a finer binning (1 ms as
%before)
%maybe create csv files to plot this all in graphpad easier.. 
loc=horzcat(DATAlab.th.mW10{:,2}); %this is the number of the cluster
a=vertcat(mxS,bslS,minSI,mxS_tm,minSI_tm,loc)';
csvwrite('summary stats.csv',a);

%other cluster info from phy output:
filename =  'cluster_info.tsv' %open and copy paste as text into a new csv in excel (havent fixed a function to do this automatically yet)

% 
% % csvwrite('clusterinfo.csv',filename)
% [cids, amp, ch, depth, fr, group, n_spikes] = readClusterInfoCSV(filename)
% 


stats=readtable('summary stats2.csv');

%% Try another way to find responders without explicit thresholding 

% 
% colPSTH= cat(1,DATA.th.mW10{:,2}); %this is spike/sec within the timebin (per 10ms timebin)

%KV edit 102323: can't find a responsewindow/timewindow that works for all
%animals to figure out responders. Try just firing rate difference (take
%baseline window and compare FR to a response window)

clear all
close all;
% path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20231002_PF optrode\100223_40trialsLP\z2023-10-02_13-53-54\'; %worked more or less well for this (with 5 sec windows)
path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230929_PF optrode\09292023_40trialsLP\z2023-09-29_14-00-01\';
load(fullfile(path,'trialOrder20230820_v2.mat'),'trialInd','trialInt');%
spikeStruct = loadKSdir(path); %need this to load the next lines
% preStimWin=0
% postStimWin=1 %500 ms after onset

preStimWin=0
postStimWin=5 %5 sec after onset

window = [-preStimWin postStimWin];
bin_size   = .01; %10 ms

goodCs = find(spikeStruct.cgs==2); %the good clusters
FR = struct([]);

[~, events_ts, events_info] = load_open_ephys_data_faster(fullfile(path,'all_channels.events')); %load stim times
[~,t1] = load_open_ephys_data_faster(fullfile(path,'Recording\100_CH1.continuous')); % to find stim onsets from the open Ephys time vector
tsStart = get_OpenEphys_recStart(path); %this uses the messages.events file; may be used to get the correct offset for spikes and events to align
pulseOnsetsOld = events_ts(events_info.eventId==1);% eventId is 1 for onset and 0 for offset of each pulse
% pulseOnsets = events_ts(events.id==1);% eventId is 1 for onset and 0 for offset of each pulse
trainIndx = [1; 1+ find(diff(pulseOnsetsOld)>0.1)]; %all train onsets %this isn't working for non trains... KV edit 8/23/23
%         trainIndx = [1; 1 + find(diff(pulseOnsets)>6)];
trainOnsets = pulseOnsetsOld(trainIndx);
eventTimes = trainOnsets;

for c = 1:length(goodCs) %go through one "good"-labeled cluster at a time
    spike_id = spikeStruct.cids(goodCs(c)); %the cluster name for this recording
    stInd = find(spikeStruct.clu==spike_id); %spike time indexing for this particular cluster
    st = spikeStruct.st(stInd); %sike times for this cluster

    eventInd = eventTimes(trialInd(1,:));
    % sync = tsStart(find(eventInd(1)>tsStart,1,'last'));
    sync = t1(1); % EA aligns to first timestamp of recording
    eventInd = eventInd-sync;%KV edit here 10/11/23: since this works, maybe can export the event times


    [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, eventInd, window, bin_size);

    FRrespwin(c).psth = psth;
    FRrespwin(c).bins = bins;
    FRrespwin(c).rasterX = rasterX;
    FRrespwin(c).rasterY=rasterY;
    FRrespwin(c).spikeCounts=spikeCounts;

end

%to get the FR divide the psth by the timebins

%Got the FR during the response window of interest for each trial for each good unit (in FR(c).spikeCounts)
%now need average FR across trials for each unit. save that. then calculate
%the FR for the baseline, take average, save. compare and find signficant
%units

% 
% preStimWin=1.1
% postStimWin=-0.1 %1 min of baseline before stim onset
%try 5 seconds like in rasters


preStimWin=4.9
postStimWin=-0.1 %1 min of baseline before stim onset

window = [-preStimWin postStimWin];

for c = 1:length(goodCs) %go through one "good"-labeled cluster at a time
    spike_id = spikeStruct.cids(goodCs(c)); %the cluster name for this recording
    stInd = find(spikeStruct.clu==spike_id); %spike time indexing for this particular cluster
    st = spikeStruct.st(stInd); %sike times for this cluster

    eventInd = eventTimes(trialInd(1,:));
    % sync = tsStart(find(eventInd(1)>tsStart,1,'last'));
    sync = t1(1); % EA aligns to first timestamp of recording
    eventInd = eventInd-sync;%KV edit here 10/11/23: since this works, maybe can export the event times


    [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, eventInd, window, bin_size);

    FRblsn(c).psth = psth;
    FRblsn(c).bins = bins;
    FRblsn(c).rasterX = rasterX;
    FRblsn(c).rasterY=rasterY;
    FRblsn(c).spikeCounts=spikeCounts;

end


%take the average across trials:
meanFR_total_respwin=[];
SEM_total_respwin=[];

for i=1:length(goodCs) %same as size of FR.spikecounts

    meanFR=mean(FRrespwin(i).psth/bin_size) %this is average across 40 trials for this unit
    meanFR_total_respwin=[meanFR_total_respwin, meanFR];

    SEM_respwin = std(FRrespwin(i).psth/bin_size)./sqrt(length(FRrespwin(i).psth)); %%edit this!
    SEM_total_respwin=[SEM_total_respwin, SEM_respwin];

end

clear meanFR
meanFR_total_blsn=[];
SEM_total_blsn=[];
for i=1:length(goodCs) %same as size of FR.spikecounts

    meanFR=mean(FRblsn(i).psth/bin_size) %this is average across 40 trials for this unit
    meanFR_total_blsn=[meanFR_total_blsn, meanFR];

    SEM_blsn = std(FRblsn(i).psth/bin_size)./sqrt(length(FRblsn(i).psth)); %%edit this!
    SEM_total_blsn=[SEM_total_blsn, SEM_blsn];

end

a=vertcat(meanFR_total_blsn, SEM_total_blsn, meanFR_total_respwin, SEM_total_respwin)';
csvwrite('summary blsn vs respwin FR.csv',a);
clear a


%split up the psth and organize by trials to do multiple t tests or whatever other stats comparison - not sure how to do this.. 
% would need to recalculate the psth for each trial? 




% for i=1:length(a)
% [h,p,ci,stats] =ttest(a(i,1), a(i,2)) %this doesnt work 
% end


%rarrange the cellaraays: save and export to prism for analysis
% blsnFR_sum=cat(1,FRblsn(:).psth);
% respFR_sum=cat(1,FRrespwin(:).psth); 

%% try this a third way (slightly different, taken from the psthViewer function)


sp=loadKSdir(myKsDir)
calcWindow = [-5 0]; % look at spike times from 3 sec before each event to  sec after
p.startRange = calcWindow(1);
p.stopRange = calcWindow(2);

eventTimes=DATA.th.mW10{1, 7};  %this is from the DATA structure of the recording of interest (is the same for all clusters)
trialGroups = ones(size(eventTimes)); 
binSize=0.01;
st=sp.st;

goodCs = find(sp.cgs==2); %the good clusters

for c = 1:length(goodCs) 

 [psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(st, eventTimes, calcWindow, binSize);

    FRblsn(c).psth = psth;
    FRblsn(c).bins = bins;
    FRblsn(c).rasterX = rasterX;
    FRblsn(c).rasterY=rasterY;
    FRblsn(c).spikeCounts=spikeCounts;
    FRblsn(c).ba=ba;



end
%%
% trGroupLabels = myData.trGroupLabels;
% nGroups = myData.nGroups;
inclRange = bins>p.startRange & bins<=p.stopRange;
spikeCounts2 = sum(ba(:,inclRange),2)./(p.stopRange-p.startRange);
meanspikes=mean(spikeCounts2)

tuningCurve = zeros(nGroups,2);
for g = 1:nGroups
    theseCounts = spikeCounts(myData.trGroups==trGroupLabels(g));
    tuningCurve(g,1) = mean(theseCounts);
    tuningCurve(g,2) = std(theseCounts)./sqrt(length(theseCounts));
end


%% SUMMARY HEATMAPS: plot heatmap of all units across all animals? try arranging anterior to posterior (y axis): 10/04/23 edit more 
%update 10/16/23: will need to change these data files (updated sorting)
%Anterior most: 9/15/23
% then :9/20, 9/18, 9/29/ and 10/02
%use the lowpass filtered data with the most units, add all to path first

data0915 = load('DATA_202310241047.mat'); %9/15
data0920 = load('DATA_202310241056.mat'); %9/20
data0918 = load('DATA_202310241054.mat'); %9/18
data0929 = load('DATA_202310241103.mat'); %9/29
data1002 = load('DATA_202310241112.mat'); %10/02

data_total = [data0915, data0920, data0918, data0929, data1002];

load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20231002_PF optrode\100223_40trialsLP\thMeasures202310241112.mat'); %REMEMBER TO LOAD THE RIGHT MEASURES FILE!
group1_total=[];

for i=1:length(data_total) %5 recordings
    group1=cat(1,data_total(i).DATA.th.mW10{:,4});
    group1_total =[group1_total;group1];
end


data0913 = load('DATA_202310241045.mat'); %0913 gfp female
data0928 = load('DATA_202310241036.mat'); %0928 gfp male
data0912 = load('DATA_202310241043.mat'); %912 female no virus

clear group1
load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230912_PF optrode\09122023_40trialsLP\thMeasures202310241043.mat'); %REMEMBER TO LOAD THE RIGHT MEASURES FILE!
group1_totalgfp=[];
data_totalgfp = [data0913, data0928, data0912];

for i=1:length(data_totalgfp) %5 recordings
    group1=cat(1,data_totalgfp(i).DATA.th.mW10{:,4});
    group1_totalgfp =[group1_totalgfp;group1];
end


HTMP = figure;
% subplot(2,1,1)
% imagesc(m.timeAxis,1:size(group1_totalgfp,1),group1_totalgfp);axis xy;colorbar;
imagesc(m.timeAxis,1:size(group1_total,1),group1_total);axis xy;colorbar;
colormap(viridis);
xlabel('Time (s) from stimulus onset');
ylabel('Cells');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005]);
set(HTMP,'color','w');
clim([-1, 4]);
xline(0,'color','w'); % marking stim on/off
hold on;
xline(24,'color','w');
% xlim([-0.5,0.75]);
% xlim([-5,25]);
title('10mW inhibition: ArchT');


HTMP2 = figure;
% subplot(2,1,1)
imagesc(m.timeAxis,1:size(group1_totalgfp,1),group1_totalgfp);axis xy;colorbar;
% imagesc(m.timeAxis,1:size(group1_total,1),group1_total);axis xy;colorbar;
colormap(viridis);
xlabel('Time (s) from stimulus onset');
ylabel('Cells');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005]);
set(HTMP2,'color','w');
clim([-1, 4]);
xline(0,'color','w'); % marking stim on/off
hold on;
xline(24,'color','w');
% xlim([-0.5,0.75]);
% xlim([-5,25]);
title('10mW inhibition: Control');


%% Plot Waveforms 
%KV edit 92623 this works now
%Further KV edits 10/08-10/10, no better way to plot this as of yet 
%KV edit 10/18/23: finally figured out how to just plot the waveform (for
%now without the SEM)!

% path = 'Z:\Shared\DATA\Electrophysiology\In
% vivo\AD\DualLaserExps\PF\CB_PF_BLA20230330\z2023-03-30_10-55-51'; %path
% to recording data (original)
%path = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230803_PF_optrode\z2023-08-03_15-31-47'; 
% path = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230915_PF optrode\09152023_40trials\z2023-09-15_12-26-49\';
%path = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230915_PF optrode\09152023_40trials\z2023-09-15_12-26-49\';
path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\z2023-09-20_13-17-48\';

gwfparams.dataDir = path;    % KiloSort/Phy output folder
% apD = dir(fullfile(path, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = dir([path,'/*.dat']);   
gwfparams.fileName = gwfparams.fileName.name; % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out

sp = loadKSdir(path); %need this to load the next lines
% 
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% 
% [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(path);
% gwfparams.spikeTimes=spikeTimes;
% gwfparams.spikeClusters=sp.clu; %maybe?
% 
% % gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

gwfparams.spikeTimes = ceil(sp.st(sp.clu==38)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters %this allows you to pick the unit of interest... 
gwfparams.spikeClusters = sp.clu(sp.clu==38);

wf = getWaveForms(gwfparams); %all of these values are in samples

% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform


figure;
imagesc(squeeze(wf.waveFormsMean)) %this is a heatmap of where the unit is located on the array
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); 
% caxis([-1 1]*max(abs(caxis()))/2); 
box off;
%%
figure; 
meanWf=squeeze(wf.waveFormsMean(:,28,:)); %the second dimension, pick the channel with the biggest signal based on the heatmap (should be a better way to do this)
plot(meanWf) %this works!! 
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('voltage (uV??)'); %whats the scale here?

%now figure out standard error? not working

% wf.waveForms(:,:,21,:)

% spkStd = sqrt((spkSum2/indx) - (spkMean.^2)); % standard deviation
% shadedErrorBar(1:length(spkMean),spkMean,spkStd);% plot average wave w/ error
% 
% shadedErrorBar(1:length(meanWf),meanWf,{@mean,@std},'-r',1); 
% shadedErrorBar(82,meanWf,std(meanWf),'g');


%% PLOT ALL RELEVANT INFO FOR ALL UNITS: 
%NOTE: this function closes all previous figs; plots heatmap of location,
%waveform, spike rate etc
%useful but doesnt seem to plot all units...

% path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\z2023-09-20_13-17-48\';
path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231130_DCNtest4_LP\z2023-11-30_13-22-16\';
sp = loadKSdir(path); %need this to load the next lines
load('rez.mat'); %manually load the rez file from the relevant .phy folder 
ops=rez.ops;
clusts=sp.cids; %if you've loaded the structure above, this will work (gives the indices of the good clusters)
plot_waveforms2(rez,clusts) %this works but only seems to plot some of the good clusters.. 

%other issue is the scaling for the waveforms? doens't seem right, see
%method above for better way?
% 
% filename =  'cluster_group.tsv'
% [cids, cgs] = readClusterGroupsCSV(filename)


%% Sorting quality analysis:
%KV added 10/16/23

clear all; close all;

%Note 10/16/23: the gfp/non virus data shows "responders", tried manually
%resorting to better get rid of putative multiunit activity (looked at
%correlogram/ISIview) but still some that pop up and those tend to be "responders"
%since they cross the z threshold (likely synchronized fluctuations due to ketamine? this is most obvious in 09/28 animal) - to check
%this in a more unbiased way, use an actual function instead of potentially
%cherry-picking units in phy. 

%Notes on sorting: some units marked as "noise" in certain recordings are just outside the
%bounds of the PF on the electrode too so not truly noise. as a habit,
%havent really been marking things as multiunit for ease. 

%update 10/16 - not sure how best to use this info.. gives all clusters
%(except 0). right now, manually checking for ISI violations in phy
%(anything <~10 spikes in the 0 is still allowed in with the caveat that it's likely not just one unit, most units at least visually have no
%spikes at 0 ms). 

%following this code: https://github.com/cortex-lab/sortingQuality/tree/master


% params = readKSparams(filename) %changed this function to hardcode the
% params

resultsDirectory='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230929_PF optrode\09292023_40trialsLP\z2023-09-29_14-00-01\'
isiV=isiViolations(resultsDirectory);

[clusterIDs, unitQuality, contaminationRate] = maskedClusterQuality(resultsDirectory)

filename =  'cluster_group.tsv';
[cids, cgs] = readClusterGroupsCSV(filename);

% cgs=clusterIDs
uQ=unitQuality;
cR=contaminationRate;

plotAllMeasures(cgs,uQ,cR,isiV);



% spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
% spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');
% spikeTimesPath= fullfile(resultsDirectory,'spike_times.npy');
% paramsPath= fullfile(resultsDirectory,'params.py');

% 
% Nt = length(spikeTrain); % total spike count
% D = max(spikeTrain); % duration of recording
% isis = diff(spikeTrain);
% refDur = 0.002;
% minISI =0.001;
% 
% numViolations = sum(isis <= refDur & isis>minISI); % number of observed violations
% 
% fpRate =  1 - sqrt(1 - numViolations*D/(Nt^2*(refDur-minISI)));
% 
% [fpRate, numViolations] = ISIViolations(spikeTrain, minISI, refDur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUMMARY DATA ANALYSIS: try other functions to get all the info we need from phy directly: 10/08-10/12/23 KV edits

% Following this tutorial: https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods#loading-raw-waveforms-or-spike-triggered-lfp-eg
%various functions downloaded from here:
%https://github.com/cortex-lab/spikes/blob/master/analysis/templatePositionsAmplitudes.m
%(added all to Kilosort analysis folder in matlabscripts folder on the server)
%Note: some of this works, could potentially be useful. Users should follow
%sections in order to load proper variables. 

% myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF
% optrode\z20230929_PF optrode\09292023_40trialsLP\z2023-09-29_14-00-01\'; 

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231107_DCNtest2_LP\z2023-11-07_15-58-17\';

sp = loadKSdir(myKsDir)

%sp.st are spike times in seconds, and sp.clu are cluster identities.
%sp.cgs are the "cluster groups", i.e. the labels that you gave to the clusters during manual sorting in phy (1=MUA, 2=Good, 3=Unsorted). 
% %Spikes from clusters labeled "noise" have already been omitted. 
% sp.cids tells you which cluster corresponds to each entry of sp.cgs, e.g. sp.cgs(sp.cids==943) will give you the cluster group for cluster 943. 

%Analyze drift over time: 
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; 
spikeYpos=spikeDepths; %maybe this is right? 

plotDriftmap(spikeTimes, spikeAmps, spikeYpos); %not sure if this is right

%Note that the channel map file used here has channel 1 at y-position=0, and since channel 1 is the site nearest the tip of the probe, this plot goes from the tip of the probe at the bottom to the most superficial part at the top.

%% 2. Quantify spiking amplitudes: see where on the probe different amplitudes were recorded and plot a colormap of the spikes across depth and amp

depthBinSize = 1; % in units of the channel coordinates, in our case this is channel number (ycoords)
depthBins = min(spikeDepths):depthBinSize:max(spikeDepths); %taken from psthByDepth code
ampBins = 0:5:min(max(spikeAmps),max(spikeAmps)); %not sure about this but seems reasonable
recordingDur = sp.st(end); %this is correct (~21 min)

% [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
% plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

ksDir=myKsDir;
[pdfs, cdfs, ampBins, depthBins] = computeAndPlotWFampHist(ksDir,ampBins,depthBins) %this updated function seems to work


[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
  templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


%% 3. Individual psths and rasters: can scroll through these to see all the clusters easily
% 
% can also load phy output from scratch here first (change the data struct so that it has events)
%need to go into datafolder of interest each time probably

%%% LOAD ONE OF THESE SETS OF DATA FILES FOR ARCHT ANIMAL OF INTEREST:
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230915_PF optrode\09152023_40trialsLP\z2023-09-15_12-26-49\';
% %load data:
load('DATA_202310161430.mat'); %9/15 %checked sorting again for ISIs

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230918_PF optrode\091823_40trialsLP\z2023-09-18_12-53-21\';
% %load data:
load('DATA_202310161435.mat'); %9/18 %checked sorting again for ISIs

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\z2023-09-20_13-17-48\';
% %load data:
load('DATA_202310161441.mat'); %9/20  %checked sorting again for ISIs

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230929_PF optrode\09292023_40trialsLP\z2023-09-29_14-00-01\';
% %load data:
load('DATA_202310161447.mat'); %9/29  %checked sorting again for ISIs

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20231002_PF optrode\100223_40trialsLP\z2023-10-02_13-53-54\';
%load data:
load('DATA_202310161455.mat'); %10/02  %checked sorting again for ISIs


%%% LOAD DATA FOR GFP:
myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230912_PF optrode\09122023_40trialsLP\z2023-09-12_14-22-34\';
load('DATA_202310161407.mat'); %0912 no virus female %checked sorting again for ISIs

myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV \in vivo\PF optrode\z20230913_PF optrode\091323_40trialsLP\z2023-09-13_12-31-59\';
load('DATA_202310161417.mat'); %0913 gfp female %checked sorting again for ISIs

myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230928_PF optrode\09282023_40trialsLP\z2023-09-28_13-50-23\'
load('DATA_202310161359.mat'); %0928 gfp male %checked sorting again for ISIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for DCN recordings (Nov 2023):
myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231130_DCNtest4_LPnolaser\z2023-11-30_13-58-49\'
load('DATA_202312011511.mat')
%%%%%%%%% PLOT PSTHS AND RASTERS FOR EACH EXPERIMENT
%for this, need to have all clusters marked as either good or noise (skips
%the noise)
%can skip above loading steps if loaded already from previous sections
sp=loadKSdir(myKsDir)
window = [-5 15]; % look at spike times from 3 sec before each event to  sec after
eventTimes=DATA.th.mW10{1, 7};  %this is from the DATA structure of the recording of interest (is the same for all clusters)
trialGroups = ones(size(eventTimes)); 
psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups); %this works 


%only plot the good clusters instead of all, this doesnt work:
% cg = [2] %good cluster
% goodCs = find(sp.cgs==cg)'
% goodCs_int=int32(goodCs)

% psthViewer(sp.st, goodCs_int, eventTimes, window, trialGroups); %this
% doesnt work

%% 4. Plot the stimulus aligned PSTH heatmap for all spikes at all depths

%basically same thing as our heatmap in the second section (more
%flexibility with the scaling since you dont have to recalculate the
%datastruct each time)

depthBinSize = 1; % in units of the channel coordinates, in this case its 1-32 (with 32 being the top of the probe)
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, eventTimes, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
