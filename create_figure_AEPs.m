%% create_figure_AEPs
% This script reads the AEP properties of the smartphone - AEP EEG
% data.
% Code was used to create Fig 4 and Supp Fig 5 and 6
% Dependencies: eeglab

%% Preparation

clearvars
close all
clc

cd(fileparts(which(mfilename)))

% Includes paths
addpath(genpath('..\FSW\'))
cd '.\Data2022\'

% Figure color props
LineColors = [0 0 255; 255 0 0; 0 255 0]./255;
AreaColors = [164 227 255; 255 137 137; 137 255 156]./255;
cmap = colormap(jet);
cmap(end/2+1,:) = [1,1,1]; % 0 is white

% Create time vector
fs = 500;                   %sample freq
PreStim = 0.2;              % time before stim [s]
PostStim = 0.5;             % time after stim [s]
t = [-PreStim:1/fs:PostStim-1/fs]*1000;     % time vec [s]

% Plot details
AEPchans = [1];       % EEG channels to use for ERP/Beta plots
topoIndices = [66, 130, 149, 163, 181, 243]; % time points used for topos (index in time vector)

% Load chanlocs
load('.\LIMO_onesample_Constant\Outdoor\LIMO.mat');
chanlocs = LIMO.data.chanlocs;

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

Ttest_AEP = load('.\LIMO_onesample_Constant\Outdoor\one_sample_ttest_parameter_1.mat');
Mask_AEP = load('.\LIMO_onesample_Constant\Outdoor\mask.mat');

%% Fig 4a: AEP

% Mean AEP (same for all params)
MeanAEP = transpose(mean(Yr_Param_Outdoor.Par_yr_raw,3));
MinAEPorg = min(min(MeanAEP))*1.2;
MaxAEPorg = max(max(MeanAEP));
PlotMinMaxOrg = max(abs([MinAEPorg,MaxAEPorg]));

% Figure with AEPs for selected channels
AEPfig = figure;
AEPfig.Name = 'AEP EEG';

Ax = gca;
for chan = 1:length(AEPchans)
    hold on
    PlotAEP = transpose(squeeze(Yr_Param_Outdoor.Par_yr_raw(AEPchans(chan),:,:)));

    options.handle     = AEPfig;
    options.color_area = AreaColors(chan,:);
    options.color_line = LineColors(chan,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = t;

    plot_areaerrorbar(PlotAEP, options)
    clear PlotAEP
end

xlabel('Time [ms]')
ylabel('[uV^2/ms]', 'Interpreter', 'tex')
grid on
h = findobj('Type','Line');
legend([flip(h)],[num2str(AEPchans')])
ylim([[MinAEPorg,MaxAEPorg]])

% Figure indicating channel locations (inset Fig 4a)
PlotChanLoc = figure;
topoplot([],chanlocs,'plotchans',AEPchans,'style','blank','electrodes','on','emarker',{'o','r',[],10},'headrad',0.55,'plotrad',0.7);
hold on
topoplot([],chanlocs,'style','blank','electrodes','on','plotchans',1:60,'headrad',0.55,'plotrad',0.7);

%% Supp Fig 5: AEP including significance patches

AEPstatfig = figure;
AEPstatfig.Name = 'AEP stat timeseries';

Ax = gca;

for chan = 1:length(AEPchans)
    hold on

    % Plot beta values
    PlotAEP = transpose(squeeze(Yr_Param_Outdoor.Par_yr_raw(AEPchans(chan),:,:)));

    options.handle     = AEPstatfig;
    options.color_area = AreaColors(chan,:);    % Blue theme
    options.color_line = LineColors(chan,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = t;

    subplot(1,length(AEPchans),chan)
    plot_areaerrorbar(PlotAEP, options)
    xlabel('Time [ms]')
    ylabel('[uV]', 'Interpreter', 'tex')
    grid on
    h = findobj('Type','Line');
    ylim([[MinAEPorg,MaxAEPorg]])

    % Mark significant time periodes with transparant background
    CurrFig = gcf;
    SigIndex = find(Mask_AEP.mask(AEPchans(chan),:)==1);    % find indices were a sig. difference exist

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
                v = [t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2); t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2)];
                face = [1 2 3 4];
                P = patch('Faces',face,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','none');
            else
                P=xline(t((SigIndex(Clusters(jj,1)))),'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.5);
            end
            uistack(P,'bottom')

        end
    end

    clear Clusters*

end

%% Fig 4c: Topos of AEP amplitudes

topoAEP = figure;
topoAEP.Name = 'AEP topos';

% For all defined time points
for ii = 1:length(topoIndices)
    subplot(2,ceil(length(topoIndices)/2),ii)
    topoplot(MeanAEP(topoIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'maplimits',[-PlotMinMaxOrg,PlotMinMaxOrg],'electrodes','off','numcontour',0);
    title([num2str((t(topoIndices(ii)))),' ms'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end
end

% Store colorbar in separate figure
ColorBarAEPtopo = figure;
ColorBarAEPtopo.Name = 'Colorbar AEP topo';
subplot(2,length(topoIndices)/2,1)
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
Fstat(~transpose(Mask_Outdoor.mask)) = 0; % Set F-statistic to zero when not significant

MinF = min(min(Fstat(topoIndices,:)));
MaxF = max(max(Fstat(topoIndices,:)));
PlotFMinMax = max(abs([MinF,MaxF]));    % Plot scale F-statistics

%% Fig 4b: timeseries of beta [Regressor: smartphone use]

Betafig = figure;
Betafig.Name = 'Beta Smartphone timeseries';

Ax = gca;
for chan = 1:length(AEPchans)
    hold on

    % Plot beta values
    PlotBeta = transpose(squeeze(Yr_Beta_Outdoor.Par_yr_raw(AEPchans(chan),:,:)));

    options.handle     = Betafig;
    options.color_area = AreaColors(chan,:);    % Blue theme
    options.color_line = LineColors(chan,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = t;

    subplot(1,length(AEPchans),chan)
    plot_areaerrorbar(PlotBeta, options)
    ylim([-PlotMinMax,PlotMinMax])
    xlabel('Time [ms]')
    ylabel('\beta', 'Interpreter', 'tex')
    grid on
    clear PlotBeta

    % Mark significant time periodes with transparant background
    CurrFig = gcf;
    SigIndex = find(Mask_Outdoor.mask(AEPchans(chan),:)==1);    % find indices were a sig. difference exist

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
                v = [t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2); t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2)];
                face = [1 2 3 4];
                P = patch('Faces',face,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','none');
            else
                P=xline(t((SigIndex(Clusters(jj,1)))),'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.5);
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

for ii = 1:length(topoIndices)

    % Beta topos
    figure(Param_BetaFig)
    subplot(2,length(topoIndices)/2,ii)
    topoplot(Outdoor_Beta(topoIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0);
    title([num2str((t(topoIndices(ii)))),' ms'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

    % F-statistic topo
    figure(Param_FstatFig)
    subplot(2,length(topoIndices)/2,ii)
    topoplot(Fstat(topoIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','maplimits',[0,PlotFMinMax],'pmask',Mask_Outdoor.mask(:,topoIndices(ii)),'colormap',cmap(end/2+1:end,:));
    title([num2str((t(topoIndices(ii)))),' ms'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

end

% Store colorbar in separate figure
ColorBarBetatopo = figure;
ColorBarBetatopo.Name = 'Colorbar Beta Use topo';
subplot(2,length(topoIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([-PlotMinMax,PlotMinMax])
colormap jet

ColorBarFstatTopo = figure;
ColorBarFstatTopo.Name = 'Colorbar Fstat Use topo';
subplot(2,length(topoIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,PlotFMinMax])
colormap(cmap(end/2+1:end,:))

%% Beta and F-values [Regressor: time]

% Determine color scale ranges for topos
Outdoor_Beta_time = transpose(mean(Yr_Beta_Outdoor_time.Par_yr_raw,3));

MinBeta1 = min(min(Outdoor_Beta_time));
MaxBeta1 = max(max(Outdoor_Beta_time));
PlotMinMax = max(abs([MinBeta1,MaxBeta1])); % Plot scale betas

% Calculate F-statistic
Fstat_time = transpose(squeeze(Ttest_Outdoor_time.one_sample(:,:,4)).^2);
Fstat_time(~transpose(Mask_Outdoor_time.mask)) = 0; % Set F-statistic to zero when not significant

MinF = min(min(Fstat_time(topoIndices,:)));
MaxF = max(max(Fstat_time(topoIndices,:)));
PlotFMinMax = max(abs([MinF,MaxF]));    % Plot scale F-statistics


%% Supp Fig 6a: Plot topos for beta and F-statistic [Regressor: time]

Betafig_time = figure;
Betafig_time.Name = 'Beta Time timeseries';

Ax = gca;
for chan = 1:length(AEPchans)
    figure(Betafig)
    hold on

    % Plot beta values
    PlotBeta = transpose(squeeze(Yr_Beta_Outdoor_time.Par_yr_raw(AEPchans(chan),:,:)));

    options.handle     = Betafig;
    options.color_area = AreaColors(2,:);    % Blue theme
    options.color_line = LineColors(2,:);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'c95';
    options.x_axis     = t;

    subplot(1,length(AEPchans),chan)
    plot_areaerrorbar(PlotBeta, options)
    ylim([-PlotMinMax,PlotMinMax])
    xlabel('Time [ms]')
    ylabel('\beta', 'Interpreter', 'tex')
    grid on
    clear PlotBeta

    % Mark significant time periodes with transparant background
    CurrFig = gcf;
    SigIndex = find(Mask_Outdoor_time.mask(AEPchans(chan),:)==1);    % find indices were a sig. difference exist

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
                v = [t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(1); t(max(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2); t(min(SigIndex(Clusters(jj,:)))) CurrFig.Children(1).YLim(2)];
                face = [1 2 3 4];
                P = patch('Faces',face,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','none');
            else
                P=xline(t((SigIndex(Clusters(jj,1)))),'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.5);
            end
            uistack(P,'bottom')

        end
    end

    clear Clusters*

end

%% Supp Fig 6b,c: Plot topos for beta and F-statistic [Regressor: time]

ParamTime_BetaFig = figure;
ParamTime_BetaFig.Name = 'Beta Param Time topos';
ParamTime_FstatFig = figure;
ParamTime_FstatFig.Name = 'Fstat Param Time topos';

for ii = 1:length(topoIndices)

    % Beta topos
    figure(ParamTime_BetaFig)
    subplot(2,ceil(length(topoIndices)/2),ii)
    topoplot(Outdoor_Beta_time(topoIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','numcontour',0);
    title([num2str((t(topoIndices(ii)))),' ms'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

    % F-statistic topo
    figure(ParamTime_FstatFig)
    subplot(2,ceil(length(topoIndices)/2),ii)
    topoplot(Fstat_time(topoIndices(ii),:),chanlocs,'headrad',0.55,'plotchans',1:60,'electrodes','off','maplimits',[0,PlotFMinMax],'pmask',Mask_Outdoor.mask(:,topoIndices(ii)),'colormap',cmap(end/2+1:end,:));
    title([num2str((t(topoIndices(ii)))),' ms'])

    FigHandles = gcf;
    for jj = 1:4
        FigHandles.Children(end-((ii-1))).Children(jj).LineWidth = 1;
    end

end

% Store colorbar in separate figure
ColorBarBetaTimetopo = figure;
ColorBarBetaTimetopo.Name = 'Colorbar Beta Time topo';
subplot(2,length(topoIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([-PlotMinMax,PlotMinMax])
colormap jet

ColorBarFstatTimeTopo = figure;
ColorBarFstatTimeTopo.Name = 'Colorbar Fstat Time topo';
subplot(2,length(topoIndices)/2,1)
colorbar('FontSize',8)
set(gca,'XColor', 'none','YColor','none','Color','none')
axis square
caxis([0,PlotFMinMax])
colormap(cmap(end/2+1:end,:))




