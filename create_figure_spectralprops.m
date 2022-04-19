%% create_figure_spectralprops
% This script reads the spectral properties of the smartphone - AEP EEG
% data.
% Code was used to create Fig 2, 3 and Supp Fig 1-4
% Dependencies: eeglab

%% Preparation

clearvars
close all
clc

cd(fileparts(which(mfilename)))

% Includes paths
addpath(genpath('..\FSW\'))
cd '.\Data2022_FOOF\'

% Figure color props
LineColors = [0 0 255; 255 0 0; 0 255 0]./255;
AreaColors = [164 227 255; 255 137 137; 137 255 156]./255;
cmap = colormap(jet);
cmap(end/2+1,:) = [1,1,1]; % 0 is white

% Plot details
Powchans = [1];                 % EEG channels to use for Pow Spectrum/Beta plots
topoFreqs = [5,10,15,20,25,32]; % Freqs used for topoplots [Hz]

% Load chanlocs
load('.\LIMO_onesample_Constant\Outdoor\LIMO.mat');

% Find indices in freq vector for topoplot frews
topoFreqIndices=zeros(size(topoFreqs));
for ii = 1:length(topoFreqs)
    [~,Freq] = min(abs(LIMO.Freq-topoFreqs(ii)));
    topoFreqIndices(ii) = Freq;
end

chanlocs = LIMO.data.chanlocs;      % EEG channel locations
f = LIMO.Freq;                      % freq vec [Hz]

%% Load required data
% Load power spectra
Yr_Param_Outdoor = load('.\LIMO_onesample_Constant\Outdoor\Yr.mat');

% Load betas
Yr_Beta_Outdoor = load('.\LIMO_onesample\Outdoor\Yr.mat');
Yr_Beta_Outdoor_time = load('.\LIMO_onesample_time\Outdoor\Yr.mat');

% Load one sample t-test results
Ttest_Outdoor = load('.\LIMO_onesample\Outdoor\one_sample_ttest_parameter_1.mat');
Ttest_Outdoor_time = load('.\LIMO_onesample_time\Outdoor\one_sample_ttest_parameter_1.mat');

% Load masks
Mask_Outdoor = load('.\LIMO_onesample\Outdoor\mask.mat');
Mask_Outdoor_time = load('.\LIMO_onesample_time\Outdoor\mask.mat');

% FOOOF details & statistics
Stats = load('.\LIMO_onesample_foofExponentShift\Outdoor\foofexpstats.mat');

%% Figa 2a: power spectra

% Mean power spectrum (same for all params)
MeanPow = transpose(mean(Yr_Param_Outdoor.Par_yr_raw,3));
MinPoworg = min(min(MeanPow));
MaxPoworg = max(max(MeanPow));

PlotMinMaxOrg = max(abs([MinPoworg,MaxPoworg]));

% Figure with power spectrum for selected channels
Powfig = figure;
Powfig.Name = 'Pow Spectrum';

Ax = gca;
for chan = 1:length(Powchans)
    hold on
    PlotPow = transpose(squeeze(Yr_Param_Outdoor.Par_yr_raw(Powchans(chan),:,:)));

    options.handle     = Powfig;
    options.color_area = AreaColors(chan,:);
    options.color_line = LineColors(chan,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = f;

    plot_areaerrorbar(PlotPow, options)
    clear PlotPow
end

xlabel('Frequency [Hz]')
ylabel('[uV^2/Hz]', 'Interpreter', 'tex')
grid on
h = findobj('Type','Line');
legend([flip(h)],[num2str(Powchans')])

ylim([[MinPoworg,MaxPoworg]])
xlim([1 40])

% Figure indicating channel locations (inset Fig 2a)
PlotChanLoc = figure;
topoplot([],chanlocs,'plotchans',Powchans,'style','blank','electrodes','on','emarker',{'o','r',[],10},'headrad',0.55,'plotrad',0.7);
hold on
topoplot([],chanlocs,'style','blank','electrodes','on','plotchans',1:60,'headrad',0.55,'plotrad',0.7);

%% Fig 2c: Topos of power spectra

topoPow = figure;
topoPow.Name = 'Pow topos';

% For all defined frequencies
for ii = 1:length(topoFreqIndices)

    subplot(2,ceil(length(topoFreqIndices)/2),ii)
    topoplot(MeanPow(topoFreqIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'maplimits',[-PlotMinMaxOrg,PlotMinMaxOrg],'electrodes','off','numcontour',0);
    title([num2str(round(f(topoFreqIndices(ii)))),' Hz'])

    FigHandles = gcf;
    for jj = 1:4 % for all children
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

end % frequencies

% Store colorbar in separate figure
ColorBarPowtopo = figure;
ColorBarPowtopo.Name = 'Colorbar Pow topo';
subplot(2,length(topoFreqIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([-PlotMinMaxOrg,PlotMinMaxOrg])
colormap jet

%% Beta and F-values [Regressor: smartphone use]

% Determine color scale ranges for topos
% Betas
Outdoor_Beta = transpose(mean(Yr_Beta_Outdoor.Par_yr_raw,3));

MinBeta1 = min(min(Outdoor_Beta));
MaxBeta1 = max(max(Outdoor_Beta));
PlotMinMax = max(abs([MinBeta1,MaxBeta1])); % Plot scale betas

% Calculate F-statistic
Fstat = transpose(squeeze(Ttest_Outdoor.one_sample(:,:,4)).^2);
Fstat(~transpose(Mask_Outdoor.mask)) = 0;   % Set F-statistic to zero when not significant

MinF = min(min(Fstat(topoFreqIndices,:)));
MaxF = max(max(Fstat(topoFreqIndices,:)));
PlotFMinMax = max(abs([MinF,MaxF]));    % Plot scale F-statistics


%% Fig 2b: timeseries of beta [Regressor: smartphone use]

Betafig = figure;
Betafig.Name = 'Beta Smartphone timeseries';

Ax = gca;
for chan = 1:length(Powchans)
    hold on

    % Plot beta values
    PlotBeta = transpose(squeeze(Yr_Beta_Outdoor.Par_yr_raw(Powchans(chan),:,:)));

    options.handle     = Betafig;
    options.color_area = AreaColors(chan,:);    % Blue theme
    options.color_line = LineColors(chan,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = f;

    subplot(1,length(Powchans),chan)
    plot_areaerrorbar(PlotBeta, options)
    ylim([-PlotMinMax,PlotMinMax])
    xlabel('Time [ms]')
    ylabel('\beta', 'Interpreter', 'tex')
    grid on
    xlim([1 40])
    clear PlotBeta

    % Mark significant time periodes with transparant background
    CurrFig = gcf;
    SigIndex = find(Mask_Outdoor.mask(Powchans(chan),:)==1);    % find indices were a sig. difference exist

    if ~isempty(SigIndex)
        ClustersInd = find(diff(SigIndex)>1);                       % check where transition from non-sig to sig (or vice versa) occures
        ClustersInd = [ClustersInd,length(SigIndex)];               % Add last cluster index = length number of signifcant indices

        % Create significant cluster blocks [start end] index
        for kk = 1:length(ClustersInd)
            if kk == 1
                Clusters(kk,:) = [1,ClustersInd(1)];
            elseif kk ~= 1  &&  kk < length(ClustersInd)
                Clusters(kk,:) =  [ClustersInd(kk-1)+1,ClustersInd(kk)];
            else
                Clusters(kk,:) = [ClustersInd(kk-1)+1,ClustersInd(end)];
            end
        end

        % Draw patch on figure for each significant time period
        for jj = 1:size(Clusters,1)

            if diff(Clusters(jj,:))
                v = [f(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); f(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); f(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2); f(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2)];
                face = [1 2 3 4];
                P = patch('Faces',face,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','none');
            else
                P=xline(f((SigIndex(Clusters(jj,1)))),'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.5);
            end
            uistack(P,'bottom')

        end
    end

    clear Clusters*

end

%% Fig 2d: Plot topos for beta and F-statistic [Regressor: smartphone use]

Param_BetaFig = figure;
Param_BetaFig.Name = 'Beta Param Smartphone topos';
Param_FstatFig = figure;
Param_FstatFig.Name = 'Fstat Param Smartphone topos';

for ii = 1:length(topoFreqIndices)

    % Beta topos
    figure(Param_BetaFig)
    subplot(2,length(topoFreqIndices)/2,ii)
    topoplot(Outdoor_Beta(topoFreqIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0);
    title([num2str(round(f(topoFreqIndices(ii)))),' Hz'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

    % F-statistic topo
    figure(Param_FstatFig)
    subplot(2,length(topoFreqIndices)/2,ii)
    topoplot(Fstat(topoFreqIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','maplimits',[0,PlotFMinMax],'pmask',Mask_Outdoor.mask(:,topoFreqIndices(ii)),'colormap',cmap(end/2+1:end,:));
    title([num2str(round(f(topoFreqIndices(ii)))),' Hz'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

end

% Store colorbar in separate figure
ColorBarBetatopo = figure;
ColorBarBetatopo.Name = 'Colorbar Beta Use topo';
subplot(2,length(topoFreqIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([-PlotMinMax,PlotMinMax])
colormap jet

ColorBarFstatTopo = figure;
ColorBarFstatTopo.Name = 'Colorbar Fstat Use topo';
subplot(2,length(topoFreqIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,PlotFMinMax])
colormap(cmap(end/2+1:end,:))

%% Beta and F-values [Regressor: time]

% Determine color scale ranges for topos
% Betas
Outdoor_Beta_time = transpose(mean(Yr_Beta_Outdoor_time.Par_yr_raw,3));

MinBeta1 = min(min(Outdoor_Beta_time));
MaxBeta1 = max(max(Outdoor_Beta_time));
PlotMinMax = max(abs([MinBeta1,MaxBeta1])); % Plot scale betas

% Calculate F-statistic
Fstat_time = transpose(squeeze(Ttest_Outdoor_time.one_sample(:,:,4)).^2);
Fstat_time(~transpose(Mask_Outdoor_time.mask)) = 0; % Set F-statistic to zero when not significant

MinF = min(min(Fstat_time(topoFreqIndices,:)));
MaxF = max(max(Fstat_time(topoFreqIndices,:)));
PlotFMinMax = max(abs([MinF,MaxF]));    % Plot scale F-statistics


%% Fig 2b: Plot topos for beta and F-statistic [Regressor: time]

Betafig_time = figure;
Betafig_time.Name = 'Beta Time timeseries';

Ax = gca;
for chan = 1:length(Powchans)
    figure(Betafig)
    hold on

    % Plot beta values
    PlotBeta = transpose(squeeze(Yr_Beta_Outdoor_time.Par_yr_raw(Powchans(chan),:,:)));

    options.handle     = Betafig;
    options.color_area = AreaColors(2,:);    % Blue theme
    options.color_line = LineColors(2,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = f;

    subplot(1,length(Powchans),chan)
    plot_areaerrorbar(PlotBeta, options)
    ylim([-PlotMinMax,PlotMinMax])
    xlabel('Time [ms]')
    ylabel('\beta', 'Interpreter', 'tex')
    grid on
    xlim([1 40])
    clear PlotBeta

    % Mark significant time periodes with transparant background
    CurrFig = gcf;
    SigIndex = find(Mask_Outdoor_time.mask(Powchans(chan),:)==1);    % find indices were a sig. difference exist

    if ~isempty(SigIndex)
        ClustersInd = find(diff(SigIndex)>1);                       % check where transition from non-sig to sig (or vice versa) occures
        ClustersInd = [ClustersInd,length(SigIndex)];               % Add last cluster index = length number of signiciant indices

        % Create significant cluster blocks [start end] index
        for kk = 1:length(ClustersInd)
            if kk == 1
                Clusters(kk,:) = [1,ClustersInd(1)];
            elseif kk ~= 1  &&  kk < length(ClustersInd)
                Clusters(kk,:) =  [ClustersInd(kk-1)+1,ClustersInd(kk)];
            else
                Clusters(kk,:) = [ClustersInd(kk-1)+1,ClustersInd(end)];
            end
        end

        % Draw patch on figure for each significant time period
        for jj = 1:size(Clusters,1)

            if diff(Clusters(jj,:))
                v = [f(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); f(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); f(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2); f(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2)];
                face = [1 2 3 4];
                P = patch('Faces',face,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','none');
            else
                P=xline(f((SigIndex(Clusters(jj,1)))),'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.5);
            end
            uistack(P,'bottom')

        end
    end
    clear Clusters*

end

%% Supp Fig 1: Plot topos for beta and F-statistic [Regressor: time]
ParamTime_BetaFig = figure;
ParamTime_BetaFig.Name = 'Beta Param Time topos';
ParamTime_FstatFig = figure;
ParamTime_FstatFig.Name = 'Fstat Param Time topos';

for ii = 1:length(topoFreqIndices)

    % Beta topos
    figure(ParamTime_BetaFig)
    subplot(2,ceil(length(topoFreqIndices)/2),ii)
    topoplot(Outdoor_Beta_time(topoFreqIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0);
    title([num2str(round(f(topoFreqIndices(ii)))),' Hz'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

    % F-statistic topo
    figure(ParamTime_FstatFig)
    subplot(2,ceil(length(topoFreqIndices)/2),ii)
    topoplot(Fstat_time(topoFreqIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','maplimits',[0,PlotFMinMax],'pmask',Mask_Outdoor.mask(:,topoFreqIndices(ii)),'colormap',cmap(end/2+1:end,:));
    title([num2str(round(f(topoFreqIndices(ii)))),' Hz'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

end

% Store colorbar in separate figure
ColorBarBetaTimetopo = figure;
ColorBarBetaTimetopo.Name = 'Colorbar Beta Time topo';
subplot(2,length(topoFreqIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([-PlotMinMax,PlotMinMax])
colormap jet

ColorBarFstatTimeTopo = figure;
ColorBarFstatTimeTopo.Name = 'Colorbar Fstat Time topo';
subplot(2,length(topoFreqIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,PlotFMinMax])
colormap(cmap(end/2+1:end,:))

%% Fig 3d,e / Supp Fig 2c,d / Supp Fig 3d,e & Supp Fig 4c,d: FOOOF details & stats

% Exponent for regressor smartphone use
BetaExp_Use = Stats.B_exp_usage_onesample(:,1);
FstatExp_Use = Stats.B_exp_usage_onesample(:,2).^2;

% Offset for regressor smartphone use
BetaOffset_Use = Stats.B_shift_usage_onesample(:,1);
FstatOffset_Use = Stats.B_shift_usage_onesample(:,2).^2;

% Exponent for regressor time
BetaExp_time = Stats.B_exp_time_onesample(:,1);
FstatExp_time = Stats.B_exp_time_onesample(:,2).^2;

% Offset for regressor time
BetaOffset_time = Stats.B_shift_time_onesample(:,1);
FstatOffset_time = Stats.B_shift_time_onesample(:,2).^2;

% Plot topos for beta and F-statistic
for aa = 1:4
    ParamFOOF_BetaFig(aa) = figure;
    ParamFOOF_BetaFig(aa).Name = 'Beta FOOOF Param Use topos';
end
ParamFOOF_FstatFig = figure;
ParamFOOF_FstatFig.Name = 'Fstat FOOOF Param Use topos';

BetaExpRangeUsage = [min([BetaExp_Use]),max([BetaExp_Use])];
BetaOffsetRangeUsage = [min([BetaOffset_Use]),max([BetaOffset_Use])];
BetaExpRangeTime = [min([BetaExp_time]),max([BetaExp_time])];
BetaOffsetRangeTime = [min([BetaOffset_time]),max([BetaOffset_time])];

% Beta topos (in separate figs due to colorbar differences)
figure(ParamFOOF_BetaFig(1))
ax(1) = subplot(2,3,1);
topoplot(BetaExp_Use',chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits','maxmin');
title('Beta Exp Use')
colormap(ax(1),cmap(1:end/2,:))
figure(ParamFOOF_BetaFig(2))
ax(2) = subplot(2,3,2);
topoplot(BetaOffset_Use',chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits',BetaOffsetRangeUsage);
title('Beta Offset Use')
colormap(ax(2),cmap(1:end/2,:))
figure(ParamFOOF_BetaFig(3))
ax(3) = subplot(2,3,3);
topoplot(BetaExp_time',chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits',BetaExpRangeTime);
title('Beta Exp Time')
colormap(ax(3),cmap(end/2+1:end,:))
figure(ParamFOOF_BetaFig(4))
ax(4) = subplot(2,3,4);
topoplot(BetaOffset_time',chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0,'maplimits',BetaOffsetRangeTime);
title('Beta Offset Time')
colormap(ax(4),cmap(end/2+1:end,:))

for ii = 1:4
    FigHandles = ParamFOOF_BetaFig(ii);
    for jj = 1:4
        FigHandles.Children(1).Children(jj).LineWidth = 1;
    end
end

% F-statistic topo
figure(ParamFOOF_FstatFig)
subplot(2,3,1)
topoplot(FstatExp_Use,chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','colormap',cmap(end/2+1:end,:),'maplimits',[0,max([FstatExp_Use])],'pmask',Stats.B_exp_usage_onesample_mask);
title('Fstat Exp Use')
subplot(2,3,2)
topoplot(FstatOffset_Use,chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','colormap',cmap(end/2+1:end,:),'maplimits',[0,max([FstatOffset_Use])],'pmask',Stats.B_shift_usage_onesample_mask);
title('Fstat Offset Use')
subplot(2,3,3)
topoplot(FstatExp_time,chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','colormap',cmap(end/2+1:end,:),'maplimits',[0,max([FstatExp_time])],'pmask',Stats.B_exp_time_onesample_mask);
title('Fstat Exp Time')
subplot(2,3,4)
topoplot(FstatOffset_time,chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','colormap',cmap(end/2+1:end,:),'maplimits',[0,max([FstatOffset_time])],'pmask',Stats.B_shift_time_onesample_mask);
title('Fstat Offset Time')

FigHandles = gcf;
for ii = 1:4
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end
end

% Store colorbar in separate figure
ColorBarBetaFOOF = figure;
ColorBarBetaFOOF.Name = 'Colorbar Beta FOOF';
ax(1) = subplot(2,3,1);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis(BetaExpRangeUsage)
colormap(ax(1),cmap(1:end/2,:))
ax(2) = subplot(2,3,2);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis(BetaOffsetRangeUsage)
colormap(ax(2),cmap(1:end/2,:))
ax(3) = subplot(2,3,3);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis(BetaExpRangeTime)
colormap(ax(3),cmap(end/2+1:end,:))
ax(4) = subplot(2,3,4);
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis(BetaOffsetRangeTime)
colormap(ax(4),cmap(end/2+1:end,:))

ColorBarFstatFOOF = figure;
ColorBarFstatFOOF.Name = 'Colorbar Fstat FOOF';
subplot(2,3,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,max([FstatExp_Use])])
colormap(cmap(end/2+1:end,:))
subplot(2,3,2)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,max([FstatOffset_Use])])
colormap(cmap(end/2+1:end,:))
subplot(2,3,3)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,max([FstatExp_time])])
colormap(cmap(end/2+1:end,:))
subplot(2,3,4)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,max([FstatOffset_time])])
colormap(cmap(end/2+1:end,:))
