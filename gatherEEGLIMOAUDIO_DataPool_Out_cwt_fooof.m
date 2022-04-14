%% Gather data for fooof analysis from non-foofed power spectra 

%% Output dir 
% for aperiodic signal 
dir_foofa = 'F:\Smartphone_AEP_2021\AllSubjectscwt\fooof\A'; 


%% General limo 
% General LIMO
load('\\Human_electrodes\expected_chanlocs')
load('\\channeighbstructmat')


%% input dir 

dir_cwt = '\Smartphone_AEP_2021\AllSubjectscwt'; 

dirRW = dir([dir_cwt, '\RW*']);

for d = 13:length(dirRW)

    if exist([dirRW(d).folder, '\',dirRW(d).name, '\RegOutdoor']) == 7

        % tmp_name
        tmp_name = dirRW(d).name;
        tmp_folder = [dirRW(d).folder, '\',dirRW(d).name, '\RegOutdoor'];

        % load the LIMO file
        load([tmp_folder,'\LIMO.mat']);
        var_all = LIMO.data.Cont;

        % load the Yr file
        load([tmp_folder,'\Yr.mat'])

        clear Yr_* fooffreq

        % create new Yr foof file
        ff = 0; 
        for f = 1:size(Yr,3)
            try
            tic
            parfor ch = 1:64

                settings = struct();
                fooof_results = fooof(LIMO.Freq, squeeze(Yr(ch,:,f)).^2, [0.1 40],settings, true);
                Yr_s(ch,f) = fooof_results.aperiodic_params(1);
                Yr_e(ch,f) = fooof_results.aperiodic_params(2);
                Yr_ap(ch,:,f) = fooof_results.power_spectrum-fooof_results.ap_fit;
                Yr_power(ch,:,f) = fooof_results.power_spectrum;
                Yr_ap_fit(ch,:,f) = fooof_results.ap_fit;

                fooffreq(ch,:) = fooof_results.freqs;
            end
            toc
            catch
            ff = ff+1; 
                RejIdx(ff) = f; 
                display('no fit found eliminating window')
            end
        end

        %% trim data 
        if ff> 0
            Yr_s(:,RejIdx) = [];
            Yr_e(:,RejIdx) = [];
            Yr_ap(:,:,RejIdx) = [];
            Yr_power(:,:,RejIdx) = [];
            Yr_ap_fit(:,:,RejIdx) = [];
            var_all(:,RejIdx) = [];
        end

        %% Outdoor LIMO model Aperiodic
        clear tmp_var_all tmp_EEG_all d_id LIMO
        tmp_EEG_all = Yr_ap;
        tmp_var_all = [zscore(var_all,[],2)'];

        tmp_dir_LIMO = [dir_foofa,'\',tmp_name,'\RegOutdoor\'];
        mkdir(tmp_dir_LIMO);
        cd(tmp_dir_LIMO);
        LIMO.dir = tmp_dir_LIMO;
        LIMO.Analysis = 'Regression';
        LIMO.Type = 'Channels';
        LIMO.Level = 1;
        LIMO.Freq = fooffreq(1,:);
        LIMO.data.chanlocs = expected_chanlocs;
        LIMO.data.neighbouring_matrix = channeighbstructmat;
        LIMO.bootstrap = 0;
        LIMO.design.bootstrap =0;
        LIMO.design.tfce = 1;
        LIMO.design.fullfactorial = 0;
        LIMO.data.Cont = tmp_var_all';
        save LIMO LIMO
        save Yr_power Yr_power
        save Yr_ap_fit Yr_ap_fit

        limo_random_robust(4,tmp_EEG_all,tmp_var_all',[1],0,0);
        close all; clear tmp_var_all tmp_EEG_all



        %% Outdoor regression model shift & exp model
        clear tmp_var_all tmp_EEG_all d_id LIMO
        tmp_var_all = [zscore(var_all,[],2)'];

        for s = 1:size(Yr_s,1)

             Stats{s} = fitlm(tmp_var_all, Yr_s(s,:), 'RobustOpts', 'on');

            State{s} = fitlm(tmp_var_all, Yr_e(s,:), 'RobustOpts', 'on');

        end

        save('ExponentShiftStats.mat', 'Stats', 'State', "Yr_e", "Yr_s", '-mat');



    end


end
