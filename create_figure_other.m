%% create_figure_other
% This script reads smartphone - AEP EEG data.
% Code was used to create various individual data figures.
% Dependencies: eeglab

%% Preparation

clearvars
close all
clc

cd(fileparts(which(mfilename)))

% Figure color props
LineColors = [0 0 255; 255 0 0; 0 255 0]./255;
AreaColors = [164 227 255; 255 137 137; 137 255 156]./255;
cmap = colormap(jet);

% Create time vector
fs = 500;                   %sample freq
PreStim = 0.2;              % time before stim [s]
PostStim = 0.5;             % time after stim [s]
t = [-PreStim:1/fs:PostStim-1/fs]*1000;     % time vec [s]

% Load chanlocs
load('.\Data2022\LIMO_onesample_Constant\Outdoor\LIMO.mat');
chanlocs = LIMO.data.chanlocs;

%% Fig 1a: Outline of experiment with insets of data examples

% Audio tone 
NrPulses = 30;
Fs = 40000;
IPI = 2*Fs*(rand(NrPulses,1))+Fs*1;
TotTime = cumsum(IPI/Fs);

% Inter tap intervals
a = 100; % [ms]
b = 5000; % [ms]
ITI = normrnd(300,1000,1000,1);     % Normal random distribution of random numbers (based on Ghosh data)
Indices = ITI < 100 | ITI > 5000;   % Normal range of inter tap intervals 100-5000 ms
ITI(Indices) = [];
ITIsum = cumsum(ITI)/1000;

figure;
plot(TotTime,0.25*ones(size(TotTime)),'or','MarkerSize',6,'MarkerFaceColor','r')
hold on
stem([ITIsum],0.2*ones(size(ITIsum)),'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 2,'Color','k')
xlim([0 30])

% Inset sound pulse
N = 2000;
tsound = (1:N)*(1/Fs);
freq = 100;

sound_pulse = sin(2*pi*freq*tsound).';
sound_array = [zeros(2500,1);sound_pulse;zeros(2500,1)];
figure;plot(sound_array)

% Inset FOOF power spectrum
load(['.\Data2022_FOOF\RW13\RegOutdoor\Yr.mat']);
load('.\Data2022_FOOF\LIMO_onesample_Constant\Outdoor\LIMO.mat');

f = LIMO.Freq;              % freq vec [Hz]
figure;plot(f,mean(Yr(1,:,:),3)')
xlim([1 40]) 

% Inset AEP
load(['.\Data2022\RW43\RegOutdoor\Yr.mat']);

figure;plot(t,mean(Yr(1,:,:),3)')
ylim([-2 2])

%% Fig 1b: Example of smartphone use across time

load(['.\Data2022\RW43\RegOutdoor\LIMO.mat']);

figure; 
plot([1:length(LIMO.data.Cont(1,:))],LIMO.data.Cont(1,:))
xticks([])
xlabel('Number of smartphone touches per min')
ylabel('\surd Count','Interpreter','tex')

%% Fig 2b: Inset smartphone use vs spectral power

load(['.\Data2022_FOOF\RW40\RegOutdoor\LIMO.mat']);
load(['.\Data2022_FOOF\RW40\RegOutdoor\Yr.mat']);

mdl = fitlm([LIMO.data.Cont(1,:)',LIMO.data.Cont(2,:)'],squeeze(Yr(1,101,:)),'RobustOpts','on');

figure; 
plotAdjustedResponse(mdl,'x1','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

%% Fig 3a / Supp Fig 3a: Indv example of smartphone usage over time

load(['.\Data2022_FOOF\RW25\RegOutdoor\LIMO.mat']);

figure; 
plot([1:length(LIMO.data.Cont(1,:))],LIMO.data.Cont(1,:))
xticks([])
xlabel('Number of smartphone touches per min')
ylabel('\surd Count','Interpreter','tex')

%% Fig 3b / Supp Fig 3b: examples exp and offset [Regressor: Usage]

load(['.\Data2022_FOOF\RW25\RegOutdoor\ExponentShiftStats.mat']);

% Exponent
Fig_ExpTime = figure;
plot(Yr_e(1,:))

% Offset
Fig_ShiftTime = figure;
plot(Yr_s(1,:))

%% Fig 3 / Supp Fig 2-4: Insets with correlations

% Exponent
mdl = fitlm([LIMO.data.Cont(1,:)',LIMO.data.Cont(2,:)'],squeeze(Yr_e(1,:)),'RobustOpts','on');

% Fig 3b [Regressor: Usage]
Fig_CorrExpUse = figure; 
plotAdjustedResponse(mdl,'x1','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

% Supp Fig 2a [Regressor: time]
Fig_CorrExpTime = figure; 
plotAdjustedResponse(mdl,'x2','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

% Offset
mdl = fitlm([LIMO.data.Cont(1,:)',LIMO.data.Cont(2,:)'],squeeze(Yr_s(1,:)),'RobustOpts','on');

% Supp Fig 3b [Regressor: Usage]
Fig_CorrShiftUse = figure; 
plotAdjustedResponse(mdl,'x1','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

% Supp Fig 4a [Regressor: rime]
Fig_CorrShiftTime = figure; 
plotAdjustedResponse(mdl,'x2','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

%% Fig 4b: Inset smartphone use vs AEP amplitude

load(['.\Data2022\RW40\RegOutdoor\LIMO.mat']);
load(['.\Data2022\RW40\RegOutdoor\Yr.mat']);

mdl = fitlm([LIMO.data.Cont(1,:)',LIMO.data.Cont(2,:)'],squeeze(Yr(1,149,:)),'RobustOpts','on');

figure; 
plotAdjustedResponse(mdl,'x1','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

%% Fig 3 / Supp Fig 2-4: Topo grand mean exponent and offset

% Load individual exponent and offset data
Folders = dir('Data2022_FOOF');
for ii = 7:length(Folders)

load(['.\Data2022_FOOF\',Folders(ii).name,'\RegOutdoor\ExponentShiftStats.mat']);

Exp(:,ii-6) = mean(Yr_e(:,:),2);
Shift(:,ii-6) = mean(Yr_s(:,:),2);
end

% Fig 3c, SuppFig 2b: Mean exponent
ExpMean = figure;
subplot(2,3,1)
topoplot(mean(Exp,2),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin','colormap',cmap(end/2+1:end,:));
title(['Exponent'])

FigHandles = gcf; 
for jj = 1:4
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Supp Fig 3c, Supp Fig 4b: Mean offset
ShiftMean = figure;
subplot(2,3,2)
topoplot(mean(Shift,2),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin','colormap',cmap(1:end/2,:));
title(['Shift'])

FigHandles = gcf; 
for jj = 1:4
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Store colorbar in separate figure
figure;
ax(1) = subplot(2,3,1);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(mean(Exp,2)),max(mean(Exp,2))])
colormap(ax(1),cmap(end/2+1:end,:))

ax(2) = subplot(2,3,2);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(mean(Shift,2)),max(mean(Shift,2))])
colormap(ax(2),cmap(1:end/2,:))

%% Fig 3b: Individual beta for exponent
% (data for RW25)

load('.\Data2022_FOOF\LIMO_onesample_foofExponentShift\Outdoor\foofexpstats.mat');

% Fig 3b [regressor: usage]
ExpTime = figure;
subplot(231)
topoplot(squeeze(Par_yr_raw_exp(:,8,1)),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin','colormap',cmap(1:end/2,:));
hold on
topoplot([],chanlocs,'headrad',0.43,'plotchans',1,'style','blank','electrodes','on','emarker',{'x','k',[],2},'colormap',cmap(1:end/2,:));
title(['Usage'])

FigHandles = gcf;
for jj = 1:8
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Supp Fig 2a [regressor: time]
ExpUse = figure;
subplot(232)
topoplot(squeeze(Par_yr_raw_exp(:,8,2)),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin');
hold on
topoplot([],chanlocs,'headrad',0.43,'plotchans',1,'style','blank','electrodes','on','emarker',{'x','k',[],2});
title(['Time'])

FigHandles = gcf;
for jj = 1:8
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Store colorbar in separate figure
figure;
ax(1) = subplot(2,3,1);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(squeeze(Par_yr_raw_exp(:,8,1))),max(squeeze(Par_yr_raw_exp(:,8,1)))])
colormap(ax(1),cmap(1:end/2,:))

ax(2) = subplot(2,3,2);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(squeeze(Par_yr_raw_exp(:,8,2))),max(squeeze(Par_yr_raw_exp(:,8,2)))])
colormap(ax(2),jet)

%% Supp Fig 3b/4a: Individual beta for offset
% (data for RW25)

load('.\Data2022_FOOF\LIMO_onesample_foofExponentShift\Outdoor\foofexpstats.mat');

% Supp Fig 3b [regressor: usage]
ShiftUse = figure;
subplot(231)
topoplot(squeeze(Par_yr_raw_shift(:,8,1)),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin','colormap',cmap(1:end/2,:));
hold on
topoplot([],chanlocs,'headrad',0.43,'plotchans',1,'style','blank','electrodes','on','emarker',{'x','k',[],2},'colormap',cmap(1:end/2,:));
title(['Usage'])

FigHandles = gcf;
for jj = 1:8
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Supp Fig 4a [regressor: time]
ShiftTime = figure;
subplot(232)
topoplot(squeeze(Par_yr_raw_shift(:,8,2)),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin','colormap',cmap(end/2+1:end,:));
hold on
topoplot([],chanlocs,'headrad',0.43,'plotchans',1,'style','blank','electrodes','on','emarker',{'x','k',[],2},'colormap',cmap(end/2+1:end,:));
title(['Time'])

FigHandles = gcf;
for jj = 1:8
    FigHandles.Children(1).Children(jj).LineWidth = 1;
end

% Store colorbar in separate figure
figure;
ax(1) = subplot(2,3,1);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(squeeze(Par_yr_raw_shift(:,8,1))),max(squeeze(Par_yr_raw_shift(:,8,1)))])
colormap(ax(1),cmap(1:end/2,:))
ax(2) = subplot(2,3,2);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([min(squeeze(Par_yr_raw_shift(:,8,2))),max(squeeze(Par_yr_raw_shift(:,8,2)))])
colormap(ax(2),cmap(end/2+1:end,:))

