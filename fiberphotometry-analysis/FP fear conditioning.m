%Written by KV 6/23/22 to analyze fp recordings during fear conditioning and extinction sessions. 
%added in code from Robinson 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6819083/
%methods also in Lerner 2015. 
%1/19/23: made the timesteps automated so you don't have to manually copy paste anything 
%1/25/23: standardize baseline correction for z transform 

%% 1. Load your data: 
clear all; close all;
path0 = uigetdir('Select the experimental session to analyze');% user points to a folder that contains all recorded channels for that stim & recording location     
files1 = dir([path0,'\','Fluorescence.csv']);  %finds the raw csv data (make sure naming is consistent)
A=readtable(files1.name,  'Format', 'auto'); %check this to make sure it loads correctly before continuing 

%% 2. Load the indx of relevant timestamps: 
%this is hardwired into the program through med associates, if the SG-231 TTL parameters ever change, this might also need to change. TTL set to "low" now so 0=onset, 1=offset. 
%Input1*2*0; SG onset (tone start)
%Input1*2*1; offset (shock end if training, tone end if extinction) - to get the shock onset need to subtract 2 seconds
%from this timestep later. 

indx_onset=find(strcmp('Input1*2*0;',A{:,2})) 
%onset is always going to have one more element because it signals when the session ends too, so need to delete:
indx_onset(end)=[];
indx_offset=find(strcmp('Input1*2*1;',A{:,2})) 

%% 3. Find timestamps based on indx values:
%Adjust number of trials depending on whether this is a training or extinction session: 
% N=5; %for training 
N=40; %for extinction - remember to comment this in/out depending on whether analyzing training or extinction! 

toneon_ts=[]; %tone onset timestamp
for i=1:N
    tone(i)=A.TimeStamp(indx_onset(i))
    toneon_ts=[toneon_ts, tone(i)];
end
clear tone

% shock_ts=[]; %shock onset timestep (subtract 2 seconds)
% for i=1:N
%     shock(i)=A.TimeStamp(indx_offset(i))- 2000;
%     shock_ts=[shock_ts, shock(i)];
% end

%% 4. Find the calcium and isosbestic signals: 
timestamps = table2array(A(:,1));
iso_signal1=table2array(A(:,3)); %for Channel 1
ca_signal1=table2array(A(:,4)); %for channel 1 
% iso_signal2=table2array(A(:,5)); %for Channel 2 (if you have multiple) 
% ca_signal2=table2array(A(:,6)); %for channel 2 

%Plot raw data to verify everything looks correct: 
CAISO = figure(1)
subplot(211)
plot(timestamps,iso_signal1);
title(['CH1: CB Isosbestic 410']);
xlabel('time(sec)'); ylabel('fluorescence');

subplot(212)
plot(timestamps,ca_signal1,'g');
title(['CH1: CB Calcium 470']);
xlabel('time(sec)'); ylabel('fluorescence');

saveas(CAISO,'Calcium and isosbestic signal.jpeg'); 
savefig(CAISO,'Calcium and isosbestic signal.fig');

%%%%%%%%%%%%%%%%
%% 5. Calculate dF/F and visualize:
% To calculate deltaF/F, a least-squares linear fit was applied to the 405 nm signal to align it to the 470 nm signal,
%producing a fitted 405 nm signal that was used to normalize the 470 nm as follows:
%deltaF/F = (470 nm signal - fitted 405 nm signal)/fitted 405 nm signal.
% From: https://github.com/GradinaruLab/dLight1/blob/master/FP_Session_Processing2.m

% Process FP Data:
Data= ca_signal1; %ch1 calcium signal
% use polyfit to fit least squares linear fit to the 405: 
%405 iso (data2)
%470 calcium (data1)
P = polyfit(iso_signal1,ca_signal1,1);
a = P(1); b = P(2);
fit400 = a.*iso_signal1 + b;
dF = 100*detrend((ca_signal1-fit400)./fit400); % Get dF/F in %

% Visualize "demeaned" 470 and 405 signals: 
MEANFIT = figure; subplot(2,1,1);
plot(timestamps,ca_signal1-mean(ca_signal1)); hold on; plot(timestamps,iso_signal1-mean(iso_signal1),'k');
set(gca,'FontSize',10); xlim([0 timestamps(end)]);
legend('GCaMP6 Signal, 470 nm','Iso Signal, 405 nm');
xlabel('Time (sec)','FontSize',12);
ylabel('Fluorescence (au)','FontSize',12);
title('De"mean"ed 470 nm and 405 nm signals','FontSize',12);

% 470 and fit signal:
subplot(2,1,2);
plot(timestamps,ca_signal1); hold on; plot(timestamps,fit400,'k');
set(gca,'FontSize',10); xlim([0 timestamps(end)]);
legend('470 nm','Fitted 410 nm');
xlabel('Time (sec)','FontSize',12);
ylabel('Fluorescence (au)','FontSize',12);
title('470 nm and "fitted" 405 nm signals','FontSize',12);

% dF/F calcualted using fitted 405 signal:
DFF = figure;
plot(timestamps,dF); hold on;
set(gca,'FontSize',10); xlim([0 timestamps(end)]);
xlabel('Time (sec)','FontSize',12);
ylabel('dF/F (%)','FontSize',12);
title('Processed GCaMP6 dF/F Signal','FontSize',12);

for i=1:length(toneon_ts) %plot lines for tone and shock onset timestamps (comment in shock_ts if you have them)
    xline(toneon_ts(i), 'r');
%     hold on
%     xline(shock_ts(i),'g');
end

%%%%%%%%%%%%%%%%
%% 6. Align signal to stim (tone or shock):
time=timestamps/1000;% we need timestamps in seconds first. divide everything by 1000.
dt=mean(diff(time)); %You can estimate dT based on your data as, for instance, the average/median or mode value of the differences between the time stamps.
F=1/dt; %Hz

%adjust these values as needed for your session:
preStimWin = 10*F; % duration of pre-stimulus baseline in sec. 
postStimWin= 40*F; % duration of post-stimulus response window 
sig470 = [time,dF]; %time in seconds, calcium sig in 2nd column.

tst= toneon_ts/1000; %tone onset timestamps in seconds, this is exact
tst=tst'
%for tone: 
ts_T=[];
clear a; clear y;
for a=1:length(tst)
    [ts_T(a),y(a)]=find(sig470==tst(a)); 
end
ts_T=ts_T'; 

for trl=1:length(ts_T) %time of tone onset (place in array)
tmT(trl,:) = time(ts_T(trl)-preStimWin:ts_T(trl)+postStimWin)-time(ts_T(trl));
deltaF_T(trl,:) = dF(ts_T(trl)-preStimWin:ts_T(trl)+postStimWin);
end

%%%%%%%%%%%%%%%%
%% 7. PLOT RESPONSES ALIGNED TO TONE ONSET:%%%%
% Plot the mean df/F: 
deltaF_Tmu = mean(deltaF_T,1); %mean across time
timeAxis = tmT(1,:);
AVDF=figure;plot(timeAxis,deltaF_Tmu,'r');
shadedErrorBar(timeAxis,deltaF_Tmu, std(deltaF_T,[],1)/sqrt(size(deltaF_T,1)),'lineprops','-b','patchSaturation',0.33);
ylabel('DF/F');xlabel('time(s) from tone onset');
title('Average calcium response to tone onset');
%Note: comment this in to save automatically: 
saveas(AVDF,'Average calcium response to tone onset.jpeg'); 
savefig(AVDF,'Average calcium response to tone onset.fig');

% Plot indiv trials (df/F): 
n = N; %5 trials in training, 40 for extinction
colors = parula(n);
DF10=figure;
h = plot(timeAxis, deltaF_T);
set(h, {'color'}, num2cell(colors, 2));
ylabel('DF/F');xlabel('time(s) from tone onset');
title('Trial by trial calcium response to tone onset');
lgd = legend;
saveas(DF10,'Trial by trial calcium response to tone onset.jpeg');
savefig(DF10,'Trial by trial calcium response to tone onset.fig');

% New Z-transform: z tranform each trial first and then take average: 
%baseline period: 
bsln = find(timeAxis<-0.2); %samples of signal a little before the tone onset
dfZtotal = [];
% for trl = 1:N
%     dfZ= deltaF_T(trl,:)-mean(deltaF_T(trl,:))./std(deltaF_T(trl,:)); %old
%     %version
%     dfZtotal=[dfZtotal;dfZ]; %each trial's z score
% end

% for all trials:
for trl = 1:N
    dfZ= deltaF_T(trl,:)-mean(deltaF_T(trl,bsln))./std(deltaF_T(trl,bsln)); %take baseline mean 
    dfZtotal=[dfZtotal;dfZ]; %each trial's z score
end

%%% Note: if you want to plot sets of trials (ex: first 10, 11-20, etc), recalculate the z score here and then replot below: 
% %for first 10 only:
% for trl = 1:10
%     dfZ= deltaF_T(trl,:)-mean(deltaF_T(trl,bsln))./std(deltaF_T(trl,bsln)); %take baseline mean 
%     dfZtotal=[dfZtotal;dfZ]; %each trial's z score
% end

% n = 10; %use this for groups of 10 trials, otherwise will be whatever the
% default N is for the full session (aka 5 or 40) 

n=N; %for the full session of however many trials
colors = parula(n); 
DFZ10=figure;
% h = plot(timeAxis, dfZtotal); %plots each indiv trial's z scored df/F
% 
% set(h, {'color'}, num2cell(colors, 2));
% for i=1:length(h)
%     h(i).Color = [h(i).Color 0.4];  % make each line transparent
% end
% ylabel('z score');xlabel('time(s) from tone onset');
% title('Trial by trial calcium response to tone onset');
% hold on

meanZ= mean(dfZtotal);
shadedErrorBar(timeAxis,meanZ, std(dfZtotal,[],1)/sqrt(size(dfZtotal,1)),'lineprops','-r','patchSaturation',0.33); 
% lgd = legend('trial1','trial2','trial3','trial4','trial5','average');
xline(0, '-.',{'tone onset'})
xline(20, '-.',{'omitted shock onset'})
yline(0,'.-b');
xlim([-5,30]);ylabel('z score');xlabel('time(s) from tone onset');
title('Trial by trial calcium response to tone onset');
saveas(DFZ10,'Trials and average(z score).jpeg');
savefig(DFZ10,'Trials and average (z score).fig');

%%%%%%%%%%%%%%%%
%% 8. Plot some heatmaps: 2/3/23
%trials on y axis, x axis = z scored df/F or reg df/F
HTMP = figure;
subplot(1,2,1);
imagesc(timeAxis,1:size(dfZtotal), dfZtotal);axis xy;colorbar;title('Z scored df/F');
line([0,20;0,20],[1,1;40,40],'color','w'); % marking CS on/off
ylabel('Trials');xlabel('Time (s) from CS onset');
set(gca,'tickdir','out','fontname','arial','fontsize',12);
colormap(plasma);  
