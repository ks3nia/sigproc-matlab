%%This is collaborative code created by EA, AD, and KV for multiunit analysis of in vivo recordings (OpenEphys; 32 channel neuronexus probe) 
%4/18/2021: edited by EA from OpEp_spikes_v9b, to analyze KV's opto recordings. 
%6/3/2021: v1_2 edited by AD to plot waveform after stim pulses
%6/28/21: edited by KV to figure out latency of response, fixed errors
%7/12/21: KV edit to add referencing (average across all channels)
%7/14/21: AD edits to fix channel naming and adjust which figures are plotted
%9/23/21 edits to make bins 5ms and adjut kernel accordingly. 
%Summer/Fall 2023: edited a bunch by AD and then KV for optrode recordings using neuronexus probes. 

%%%%%%%%%
clear all; close all;
path0 = uigetdir('Select the experimental session to analyze');% user points to a folder that contains all recorded channels for that stim & recording location     
%NOTE: all channels must be recorded as module 100 (e.g., %100_CH4.continuous). If not, then update the module variable below
module = '100';
files1 = dir([path0,'\',module,'_CH*.continuous']);    

% All temporal parameters below are in seconds:
bin_size   = 0.01; % width of time bin for instantaneous spike rate %10 ms 
Spkrefrac  = 0.002; % minimum accepted time difference between two consecutive spikes 
preStimWin= 5; % duration of post-stimulus response window 
postStimWin= 10; % duration of post-stimulus response window 
minIBI = 20
krnl_size  = 20; % num of bins for smoothing spike rate this should = 100ms
spkThres=-3.5; %z scored FR threshold (check each recording)

% load event timestamps (common for all channels in that folder):
try
    [~, events_ts, events_info] = load_open_ephys_data_faster([path0,'\','all_channels.events']);
catch
    disp('ERROR: Cannot load event timestamps');
    return;
end
pulseOnsets = events_ts(find(events_info.eventId));% eventId is 1 for onset and 0 for offset of each pulse

pulseOffsets = events_ts(find(events_info.eventId==0));
if ~length(pulseOnsets)==length(pulseOffsets) % a simple catch of trouble 
    disp('ERROR: Problem with pulse timestamps');
    return;
end

% Demarcate trials, so long as there's at least 10-second interval btwn trains:
% Because diff removes one value, we need to shift the indices below by one, and include the first trial as well:
trainIndx = [1; 1 + find(diff(pulseOnsets)>minIBI)]; 
trainOnsets = pulseOnsets(trainIndx)
numPulses = mean(diff(trainIndx)); % num pulses per train
trainDuration = 1; % only the first train, assuming they're all like the first %was 24 before
trainFrq = numPulses/trainDuration; % train frequency in pulses/sec

clear measures;
measures.path = path0;
measures.bin_size = bin_size;
measures.spkRefrac = Spkrefrac;
measures.krnl_size = krnl_size;
measures.stimTrainInterval = mean(diff(trainOnsets)); % inter-train (or inter-trial) interval
measures.TrainOnsets = trainOnsets; % timestamps of trains (same as trials) 
measures.preStimWin = preStimWin;
measures.postStimWin = postStimWin;
measures.spkThresZ = spkThres;
measures.sigFilterHi = 'H'; % H for High-pass @300Hz using 4th order for reduced artifact distortion
measures.sigFilterHiOrder = 4; % order of filter
measures.sigFilterHiCutoff = 300; % Hi-pass cutoff frequency
measures.sigFilterLo = 'L'; % Low-pass for LFP signal using 4th order for reduced artifact distortion
measures.sigFilterLoOrder = 4; % order of lo-pass filter
measures.sigFilterLoCutoff = 200; % Lo-pass cutoff frequency

% Start channel-specific processing:
load('Z:\Shared\DATA\Electrophysiology\In vivo\AD\CB_NODE_eStim\chanMapNNX32b.mat');
spkSig_total =[];
spkSig_total_LFP =[];

for channel=1:32 %for all channels
    FHCi=find(ismember(chanMap,channel)); % corresponding index in FHC numbering
    if ~isempty(FHCi) % for the 12/16 channels that are connected to the array
        chnI=chanMap(FHCi)+1; % to start from 1 rather than 0
    else
        disp(['ERROR: cannot find channel ',num2str(channel)]);
        continue;
    end
    
  % load openephys data:
    clear data LFP0 spkSig sig* t1
    fileName = dir([path0,'/*_CH',num2str(channel),'.continuous']);
    [data, t1, info] = load_open_ephys_data_faster([path0,'\',fileName.name]);
    measures.timeOfAcquisition = info.header.date_created; % should be the same for all channels - just save the last channel analyzed
    Fs = info.header.sampleRate; % in Hz by default
  % to isolate LFP:
    lfpSig = neuralFilter_v2(data,Fs,measures.sigFilterLo,measures.sigFilterLoOrder,measures.sigFilterLoCutoff);
    lfp2 = downsample(lfpSig,30);% downsampling to 1KHz, assuming it was 30KHz originally
    t1_lfp = downsample(t1,30);
  % to isolate spikes:   
    spkSig = neuralFilter_v2(data,Fs,measures.sigFilterHi,measures.sigFilterHiOrder,measures.sigFilterHiCutoff);
    spkSig=spkSig/1000; %convert to mV
    lfp2=lfp2/1000; %also convert to mV
    %7/12/21: KV add referencing here before thresholding and converting to z score:
    %spkSig is for all trials, per channel as it goes through the for loop..need to save it by channel though. 
    spkSig_total =[spkSig_total,spkSig]; 
    spkSig_total_LFP = [spkSig_total_LFP,lfp2];
end

%%%%%%%%%
% CONTINUE HERE and go channel by channel:
RecLocation = inputdlg('enter recording location');   
StimIntensity = inputdlg('enter stim intensity');
ArchT = inputdlg('ArchT? 1-yes 0-no');
if isequal(ArchT{1},'1')
    inhDur = inputdlg('Duration of Arch');
    measures.ArchDur = inhDur;
    inhAmp = inputdlg('Amplitude of Arch');
    measures.ArchAmp = inhAmp;
end
trialNumbrs = 1:40;
stims=trainOnsets;
measures.whichTrls = stims;
measures.ArchT = ArchT;
measures.recLocation = RecLocation;
measures.mWlaserPower = StimIntensity;

% Now go channel by channel and plot LFP/spectrogram, raster, psth, and z scored FR:
for ch11= 1:32
    channel = chanMap(ch11);
    spkSigref = spkSig_total;
    spkSig1 = spkSigref(:,ch11);
    bslnEpoch = 10*Fs; % number of samples in 10 seconds
    bslnMu = mean(spkSig1(1:bslnEpoch)); 
    bslnSD = std(spkSig1(1:bslnEpoch));
    sig2 = (spkSig1-bslnMu)./bslnSD;  % converting to z score based on initial baseline voltage for each file. Keeping spkSig in voltage and sig2 in z

% Thresholding spikes:
    if spkThres > 0
        spks1 = t1(sig2 >= spkThres);
    else
        spks1 = t1(sig2 <= spkThres);
    end
    spksisi = diff(spks1); % diff will shift data by one position and we must correct for that
    spks1 = [spks1(1);spks1(find(spksisi>Spkrefrac)+1)]';% keeping only threshold crossings that are sufficiently spaced

% KV edit 091923: add a spectrogram here (per channel): 
    spkSigref0=spkSig_total_LFP;
    spkSig0=spkSigref0(:,ch11);
    movingwin=[2,0.05]; % set the moving window dimensions
    params.Fs=1000; % sampling frequency
    params.fpass=[0 100]; % frequencies of interest
    params.tapers=[5 9]; % tapers
    params.trialave=0; % average over trials
    params.err=0; % no error computation
    [S1,t,f]=mtspecgramc(spkSig0,movingwin,params); % compute spectrogram
    
    figure 
    plot_matrix(S1,t,f);
    xlabel([]); % plot spectrogram
    colormap(jet);
    figure;
    p=0.05;
    params.err=[1 p] %with p=0.05, and then execute
    [S,f,Serr]=mtspectrumc(spkSig0,params);
    plot_vector(S,f,[],Serr);
    hold on;
    
    params.tapers=[5 9]
    params.fpass=[0 100]
    params.pad=2
    params.trialave=1 %(average spectrum shown in black)
    params.err=[1 .05] %in red
    params.err=[2 .05] %in red

    params.err=[2 p];
    [S,f,Serr]=mtspectrumc(spkSig0,params); % S: 1st dimension is time, f: 2nd is freqency, Serr: power for the electrode
    plot(f,10*log10(Serr(1,:)),'r');
    plot(f,10*log10(Serr(2,:)),'r');
    xlim([0, 100]);
    title('LFP spectrum average with error bars')
    ylabel('Power (dB)')
    xlabel('Frequency (Hz)')

% KV 6/15/21: to plot the signal with spike thereshold for confirmation: and look at raw voltage trace for artifacts: 
    Cnfrm = figure('WindowState','maximized','Visible','on'); % update visibility to confirm it or not -- it is not saved currently
    ThresX = axes('Parent',Cnfrm,'Position', [0.1 0.1 0.8 0.8]); 
    plot(ThresX,t1,sig2);xlabel('time (s)');ylabel('z score');
    hold on;
    for i=1:length(pulseOffsets)
        xline(pulseOffsets(i), 'k');
        hold on
    end
    hold on;
    for i=1:length(pulseOnsets)
        xline(pulseOnsets(i), 'r','Linewidth',2);
        hold on
    end
    hold on;
    line(ThresX,[t1(1) t1(end)],[spkThres spkThres ],'Color','g');  % spike threshold  
    title([path0,'\',fileName.name]);

 % EA RESULT FIGURE: to plot results in 3-panel figure
    Rslt =  figure('WindowState','maximized','Visible','on');
    RasterX = axes('Parent',Rslt,'Position',[0 0.64 1 0.31]); 
    PSTHX   = axes('Parent',Rslt,'Position',[0 0.32 1 0.31]); 
    ZtraceX = axes('Parent',Rslt,'Position',[0 0 1 0.31]); 

 % A. Rasterplot:        
    subplot(3,1,1,RasterX);
    clear dots
    instFR=[];
    rasterYoffset=0;
      for trl=1:length(stims)
          dots{trl} = spks1(spks1>=stims(trl)-preStimWin & spks1<=stims(trl)+postStimWin); % find spike timestamps within a time window of stim onset
          line([dots{trl}; dots{trl}]-stims(trl),rasterYoffset+      [0.1*ones(length(dots{trl}),1)-0.04,0.1*ones(length(dots{trl}),1)+0.04]','Color','k');%,'LineWidth',2);
          hold on;
          bins = stims(trl)-preStimWin:bin_size:stims(trl)+postStimWin;
          [N,edges] = histcounts(spks1,bins);
          instFR = cat(1,instFR,N/bin_size); % scaled to spks/sec -- (trials x time)
          rasterYoffset = rasterYoffset+0.15; % a little extra separation across the different intensities
      end
    ylabel('Trials');
    set(gca,'Xlim',[-preStimWin,postStimWin],'fontName','arial','fontSize',12,'box','off','Ytick',[0,rasterYoffset],'Yticklabel',{'1','20'},'Ycolor','k','Xcolor','k','TickDir','out'); % 'TickLength',[0.02,0.02]
    title(['Filtered and referenced data', num2str(channel)]);
    measures.channels(channel).instFRperSecond = instFR;
    
% B. PSTH: 
    subplot(3,1,2,PSTHX);
    timeAxis = -preStimWin:bin_size:postStimWin-bin_size;    
    meanFR = mean(instFR,1);
    bar(timeAxis,meanFR,'FaceColor','k');
    ylabel('Mean spikes/sec');    

% C. Mean, baseline-corrected, smoothed firing rate: 
    % 1C. Calculate average across trials and smooth:        
    krnl = hanning(krnl_size)/sum(hanning(krnl_size)); % hanning window - the kernel must be normalized to its size so that it sums to 1
    tAxis_smoothed = conv(mean(timeAxis,1),krnl,'valid');
    instFR_smoothed = conv2(instFR,krnl','valid'); % the direction of kernel ensures smoothing only in time dimension
    bslnWin = find(timeAxis<=-0.01); % baseline activity up to 10 ms before first pulse       
    bslnFR = instFR(:,bslnWin);
    bslnFR_smoothed = conv2(bslnFR,krnl','valid'); % the direction of kernel ensures smoothing only in time dimension

    % 2C. Baseline correction using resampling technique:    
    FRzScore = [];
    permutedFR = [];
    for trl = 1:length(stims) % each train
        permutedFR(trl,:,1) = instFR_smoothed(trl,:); % the observed values
        for iiy = 2:1001 % 1000 permutations of baseline activity particular to the trial
            permutedFR(trl,:,iiy) = randsample(bslnFR_smoothed(trl,:),length(instFR_smoothed(trl,:)),true); % w/ replacement
        end
        dummy = zscore(permutedFR(trl,:,:),[],3); % z transform along dimension 3 (permutations)
        FRzScore(trl,:) = squeeze(dummy(:,:,1)); % keeping every observed trial's z score            
    end
    
    % 3C. For aggregate z:
    %permutedFR: trials x time x permutations (1st is ovserved value)
    MeanTrlZscore = zscore(squeeze(mean(permutedFR,1)),[],2); % z score of trial average relative to permutations
    measures.channels(channel).aggregateZscore = MeanTrlZscore(:,1); % keeping only obsverved permutation's z score
    measures.channels(channel).trialZscore = FRzScore; % all trial z scores are saved
    
    % 4C. Quantify response:
    ResponseWin = find(tAxis_smoothed>0 & tAxis_smoothed<=postStimWin); % the acceptable response window in samples after time zero
    % Positive peak of each trial:
    [MaxZ,Latency] = max(FRzScore(:,ResponseWin),[],2); 
    Latency = tAxis_smoothed(ResponseWin(Latency)); % going from the sample of peak response to the corresponding timepoint
    measures.channels(channel).peakResp = MaxZ;
    measures.channels(channel).peakLatency = Latency;
    % Negative peak of each trial:
    [MinZ,NegLatency] = min(FRzScore(:,ResponseWin),[],2);
    NegLatency = tAxis_smoothed(ResponseWin(NegLatency)); % going from the sample of peak response to the corresponding timepoint
    measures.channels(channel).NegPeakResp = MinZ;
    measures.channels(channel).NegPeakLatency = NegLatency;
    measures.PSTHtimeAxis   = timeAxis;
    measures.zScoreTimeAxis = tAxis_smoothed;
    
    % 5C. Plot z scored FR:
    subplot(3,1,3,ZtraceX);
    plot(tAxis_smoothed,measures.channels(channel).aggregateZscore,'k');
    xlabel('Time (s) from stimulus onset');
    ylabel({'Firing Rate';'(Z score relative to baseline)'}); hold on;
    set(gca,'Xlim',[-preStimWin,postStimWin],'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out'); % 'TickLength',[0.02,0.02]
    set(gca,'color','w');
    expDate = ['2023',extractBefore(extractAfter(path0,'2023'),'_')];
    clear newName
    newName = [path0,'\',expDate,'_CH',num2str(channel),'_',RecLocation,'_',StimIntensity,'mW_',datestr(now,'yyyymmddhhMM')];
    disp(['newname: ',newName]);
    clear permutedFR FRzScore MeanTrlZscore
    
% Save overall figure with raster, psth, and z scored FR for each channel: 
        savefig(Rslt,[num2str(channel) '_' '_V1_wr.fig']);
        saveas(Rslt,[num2str(channel) '_' '__V1_wr.jpeg']);

% KV edit 6/29: add a heatmap of the rasterplot; need zscored-firing rate by trial 
    HTMP_raster = figure;
    preStimWin1= 5;
    postStimWin1=10;
    imagesc(tAxis_smoothed,1:size(measures.channels(channel).trialZscore,1),measures.channels(channel).trialZscore);axis xy;
    colorbar;
    caxis([-5, 30]); %zscore
    ylabel('Trials');
    xlabel('Time from stimulus onset (sec');
    set(gca,'Xlim',[-preStimWin1,postStimWin1],'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out'); % 'TickLength',[0.02,0.02]
    title(['Heatmap raster',num2str(channel)]);
    set(HTMP_raster,'color','w');
    saveas(HTMP_raster,[num2str(channel),'_heatmapraster_wr.jpeg']);

    clear newName
end % end of channel-specific processing

%%%%%%%%%
% Save the results:
save([path0,'\','_',measures.mWlaserPower{1},'mW_','l',RecLocation{1},'_',datestr(now,'yyyymmddhhMM'),'_V1.mat'],'measures');

%%%%%%%%%
%% KV edit to plot heatmap of all leads on the array: 
dist2tip = chanMap;
[X,I] = sort(dist2tip); % ordering the dist2tip and getting indices I
group = cat(2,measures.channels(:).aggregateZscore)';
groupOrdered = group(I,:); % concatenated aggr. Z scores ordered by dist2tip
HTMP = figure;
imagesc(measures.zScoreTimeAxis,1:size(groupOrdered,1),groupOrdered);axis xy;colorbar;
caxis([-3 5]);
xlabel('Time (s) from stimulus onset');
ylabel('Distance (mm) from tip of array');
set(gca,'fontName','arial','fontSize',12,'box','off','Ycolor','k','Xcolor','k','TickDir','out','TickLength',[0.005,0.005],'Ytick',1:size(groupOrdered,1),'YtickLabel',X);
set(HTMP,'color','w');
line([0,24;0,24],[1,1;24,24],'color','w'); % marking CS on/off
% ylabel('Trials');xlabel('Time (s) from CS onset');
% set(gca,'tickdir','out','fontname','arial','fontsize',12);
% colormap(plasma);
title(extractAfter(measures.path,'\'));
xlim([-5,10]);
saveas(HTMP,[measures.path,'_heatmaprefd_' num2str(spkThres) '_' num2str(bin_size) '.jpeg']);
savefig(HTMP,[measures.path,'_heatmaprefd' num2str(spkThres) '_' num2str(bin_size) '.fig']);
