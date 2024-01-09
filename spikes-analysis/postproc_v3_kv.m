%% Post processing for sorted single unit data from Kilosort/phy. Used for thalamic and cerebellar in vivo recordings (OpenEphys; 32 channel NNX probe).

%Create DATA structure with PSTH and Zscore
% updated to work with concatenated outputs.
% function DATA = createDATAstruct(path0) 
% path to kilosort and phy outputs (should be raw data location too)
%Edited by KV 8/21/23 for PF optrode recordings
%Edited by KV in October 2023 to add waveform properties, check sorting quality, added scrolling psth/raster plots, add auto search for responders, some quantification of firing rate etc

clear; clc;
addpath(genpath('Z:\Shared\MATLAB scripts\Ksenia\in vivo analysis\Kilosort post-processing'));
% %for all. must have folder with all folders containing KS output, both
%.events files, acqTimes.mat, and 'trialOrder.mat' file that is a cell array in proper format
path0 = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231130_DCNtest4_LPnolaser';
datafolders = dir(path0); 
datafolders = datafolders(~contains({datafolders(:).name},{'.','..','Measures','DATA'})); % if it includes any of the {}, it will exclude them. Because all filenames include .mat, it keeps only folders w/o DATA or Measures in their name.
[DATA,DATAlab] = createDATAstruct_v2_KV(path0,datafolders);
save(fullfile(path0,['DATA_' datestr(now,'YYYYmmddhhMM')]),'DATA','DATAlab');

%% 1. Plot All Rasters
%can either run this right after creating new datastruct above or load a
%previous one in the folder of interest

m.krnl_size=16;
krnl = gausswin(m.krnl_size)/sum(gausswin(m.krnl_size));

close all;
for f=1:size(DATA.th.mW10,1)
    figure;
    subplot(3,1,1); plot(DATA.th.mW10{f,5},DATA.th.mW10{f,6}); ylim([0 40]);... %this is the raster
    subplot(3,1,2); plot(DATA.th.mW10{f,3},conv(DATA.th.mW10{f,2},ones(10,1),'valid')/10);title(f); %smoothed psth
    subplot(3,1,3); plot(DATA.th.mW10{f,3},conv(DATA.th.mW10{f,9},ones(10,1),'valid')/10);title(f); %smoothed z scored psth
end

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


%% 3. Make Heatmap of all cells
load('Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\thMeasures202310161441.mat'); %don't necessarily need to change this, just need timeaxis to be correct
group1=cat(1,DATA.th.mW10{:,4}); %this loads the z scored psths for all units

HTMP = figure;
imagesc(m.timeAxis,1:size(group1,1),group1);axis xy;colorbar;
colormap(viridis);
xlabel('Time (s) from stimulus onset');
ylabel('Cells');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005]);
set(HTMP,'color','w');
title('10mW inhibition');
xline(0,'color','w'); % marking stim on/off
hold on;
xline(24,'color','w');

%% SUMMARY HEATMAPS: plot heatmap of all units across all animals. try arranging anterior to posterior (y axis)
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
title('10mW inhibition: ArchT');

HTMP2 = figure;
imagesc(m.timeAxis,1:size(group1_totalgfp,1),group1_totalgfp);axis xy;colorbar;
colormap(viridis);
xlabel('Time (s) from stimulus onset');
ylabel('Cells');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005]);
set(HTMP2,'color','w');
clim([-1, 4]);
xline(0,'color','w'); % marking stim on/off
hold on;
xline(24,'color','w');
title('10mW inhibition: Control');


%% Plot Waveforms 
path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\z2023-09-20_13-17-48\';
gwfparams.dataDir = path;    % KiloSort/Phy output folder
gwfparams.fileName = dir([path,'/*.dat']);   
gwfparams.fileName = gwfparams.fileName.name; % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
sp = loadKSdir(path); %need this to load the next lines

gwfparams.spikeTimes = ceil(sp.st(sp.clu==38)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters %this allows you to pick the unit of interest... 
gwfparams.spikeClusters = sp.clu(sp.clu==38);
wf = getWaveForms(gwfparams); %all of these values are in samples

figure;
imagesc(squeeze(wf.waveFormsMean)) %this is a heatmap of where the unit is located on the array
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); 
box off;

figure; 
meanWf=squeeze(wf.waveFormsMean(:,28,:)); %the second dimension, pick the channel with the biggest signal based on the heatmap 
plot(meanWf) 
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('voltage (uV??)'); %whats the scale here?

%% PLOT ALL RELEVANT INFO FOR ALL UNITS: 
%NOTE: this function closes all previous figs; plots heatmap of location, waveform, spike rate etc
%useful but doesnt seem to plot all units...

path='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231130_DCNtest4_LP\z2023-11-30_13-22-16\';
sp = loadKSdir(path); %need this to load the next lines
load('rez.mat'); %manually load the rez file from the relevant .phy folder 
ops=rez.ops;
clusts=sp.cids; %if you've loaded the structure above, this will work (gives the indices of the good clusters)
plot_waveforms2(rez,clusts) %this works but only seems to plot some of the good clusters.. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUMMARY DATA ANALYSIS: try other functions to get all the info we need from phy directly: 10/08-10/12/23 KV edits

% Following this tutorial: https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods#loading-raw-waveforms-or-spike-triggered-lfp-eg
%various functions downloaded from here: https://github.com/cortex-lab/spikes/blob/master/analysis/templatePositionsAmplitudes.m

myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\Test DCN\20231107_DCNtest2_LP\z2023-11-07_15-58-17\';
sp = loadKSdir(myKsDir)

%sp.st are spike times in seconds, and sp.clu are cluster identities.
%sp.cgs are the "cluster groups", i.e. the labels that you gave to the clusters during manual sorting in phy (1=MUA, 2=Good, 3=Unsorted). 
% %Spikes from clusters labeled "noise" have already been omitted. 
% sp.cids tells you which cluster corresponds to each entry of sp.cgs, e.g. sp.cgs(sp.cids==943) will give you the cluster group for cluster 943. 

%%%%%%%%%%%%%%%%%%%%
% 1. Analyze drift over time: 
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; 
spikeYpos=spikeDepths; %maybe this is right? 
plotDriftmap(spikeTimes, spikeAmps, spikeYpos); %not sure if this is right

%Note that the channel map file used here has channel 1 at y-position=0, and since channel 1 is the site nearest the tip of the probe, this plot goes from the tip of the probe at the bottom to the most superficial part at the top.

%%%%%%%%%%%%%%%%%%%%
%% 2. Quantify spiking amplitudes: see where on the probe different amplitudes were recorded and plot a colormap of the spikes across depth and amp

depthBinSize = 1; % in units of the channel coordinates, in our case this is channel number (ycoords)
depthBins = min(spikeDepths):depthBinSize:max(spikeDepths); %taken from psthByDepth code
ampBins = 0:5:min(max(spikeAmps),max(spikeAmps)); %not sure about this but seems reasonable
recordingDur = sp.st(end); %this is correct (~21 min)
ksDir=myKsDir;
[pdfs, cdfs, ampBins, depthBins] = computeAndPlotWFampHist(ksDir,ampBins,depthBins) %this updated function seems to work
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
  templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

%%%%%%%%%%%%%%%%%%%%
%% 3. Individual psths and rasters: can scroll through these to see all the clusters easily
% can also load phy output from scratch here first (change the data struct so that it has events)
%need to go into datafolder of interest each time probably

%%% LOAD ONE OF THESE SETS OF DATA FILES FOR ARCHT ANIMAL OF INTEREST:
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230915_PF optrode\09152023_40trialsLP\z2023-09-15_12-26-49\';
% %load data:
load('DATA_202310161430.mat'); %9/15 %checked sorting again for ISIs
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230918_PF optrode\091823_40trialsLP\z2023-09-18_12-53-21\';
% %load data:
load('DATA_202310161435.mat'); %9/18 
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230920_PF optrode\092023_40trialsLP\z2023-09-20_13-17-48\';
% %load data:
load('DATA_202310161441.mat'); %9/20  
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230929_PF optrode\09292023_40trialsLP\z2023-09-29_14-00-01\';
% %load data:
load('DATA_202310161447.mat'); %9/29 
myKsDir = 'Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20231002_PF optrode\100223_40trialsLP\z2023-10-02_13-53-54\';
%load data:
load('DATA_202310161455.mat'); %10/02  

%%% LOAD DATA FOR GFP:
myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230912_PF optrode\09122023_40trialsLP\z2023-09-12_14-22-34\';
load('DATA_202310161407.mat'); %0912 no virus female %checked sorting again for ISIs
myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV \in vivo\PF optrode\z20230913_PF optrode\091323_40trialsLP\z2023-09-13_12-31-59\';
load('DATA_202310161417.mat'); %0913 gfp female %checked sorting again for ISIs
myKsDir='Z:\Shared\DATA\Electrophysiology\In vivo\KV\in vivo\PF optrode\z20230928_PF optrode\09282023_40trialsLP\z2023-09-28_13-50-23\'
load('DATA_202310161359.mat'); %0928 gfp male %checked sorting again for ISIs

%Plot psths and rasters for each experiment:
%for this, need to have all clusters marked as either good or noise (skips the noise)
%can skip above loading steps if loaded already from previous sections
sp=loadKSdir(myKsDir)
window = [-5 15]; % look at spike times from 3 sec before each event to  sec after
eventTimes=DATA.th.mW10{1, 7};  %this is from the DATA structure of the recording of interest (is the same for all clusters)
trialGroups = ones(size(eventTimes)); 
psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups); %this works 

%%%%%%%%%%%%%%%%%%%%
%% 4. Plot the stimulus aligned PSTH heatmap for all spikes at all depths
%basically same thing as our heatmap in the second section (more flexibility with the scaling since you dont have to recalculate the datastruct each time)

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
