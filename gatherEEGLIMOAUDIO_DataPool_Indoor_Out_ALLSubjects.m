%% For each EEG data, gather the 1st level analysis (for inside only). 
% App category: Non-social (0) or Social (1)
% Interaction Density Continous 


% Directory settings 
data_in = '\Feb_2018_2020_RUSHMA_ProcessedEEG'; 
dat_out = '\Smartphone_AEP_2021\AllSubjects\';

% General LIMO 
load('\Human_electrodes\expected_chanlocs')
load('\Human_electrodes\channeighbstructmat')

% Load ICA config set by Mark 
load('ICinfo_audio_bilat.mat')

% Load RW excel configured by Roderick and Mark 
[~,~,rawsheet] = xlsread('\DataDetailList_MvdR.xlsx');
rawsheet(14,:) = []; % remove the RW28 recording which did not work out (no data saved on 29/11/2019 on deposit folder). 

ICinfo([13 19 21]) = [];% rejections (5) or data error 

%% Now loop through to perform LIMO regressions 
for s = 1:length(ICinfo) % 13 19 21
    %try
    % Gather times of EEG and trim eeg_name to EEG elevant eeg_name
    [~,tmp_name] = fileparts(fileparts(ICinfo(s).path));
    
    % Get EEG data
    tmp_EEG_all = [];
    tmp_var_all = [];
    
    tmp_dir =  strcat(data_in,'\', tmp_name);
    EEG = pop_loadset(strcat(ICinfo(s).processed_name,'.set'), tmp_dir) ;
    EEG = gettechnincallycleanEEG(EEG);

    
    % Onset of EEG in UTC
    dvec_t = EEG.urevent(1).bvtime;
    dvec =[dvec_t{1,1}(1:4),'\' , dvec_t{1,1}(5:6), '\', dvec_t{1,1}(7:8), '\', dvec_t{1,1}(9:10), '\', dvec_t{1,1}(11:12), '\', dvec_t{1,1}(13:14)];
    dvec_ms = str2num(dvec_t{1,1}(15:end))./1000;
    onset = posixtime(datetime(dvec,'InputFormat','yyyy\MM\dd\HH\mm\ss','TimeZone', 'Europe/Amsterdam'))+(dvec_ms/1000);
    
    %% loop through in 5 minute windows to create an erp matrix ch x t x window, and with the properties noted in the window (indoor, outdoor
    minu = 1; % number of min
    tmp_sample = EEG.srate;% samplerate
    tmp_idx = strcmp({EEG.event.type},'M  1');
    tmp_stimstmap = [EEG.event(tmp_idx).latency];
    dur = [tmp_stimstmap(end)-tmp_stimstmap(1)]./(EEG.srate*60*minu);   % number of 1 minutes
    
    g = 0; tmdatamatrix = []; Description = [];
    for m = 1:(floor(dur)-1)
        % List of events
        tstamp = tmp_stimstmap(1)+[(minu*1000*60*m-1)*(EEG.srate/1000)];
        sidx = [];
        sidx = and(tmp_stimstmap>tstamp, tmp_stimstmap < (tstamp+[(minu*1000*60)*(EEG.srate/1000)])) ;
        
        if sum(sidx) > 10 % a minimum of 10 trials needed to proceed
            g = g+1;
            % EEG data epoched
            [epochdat, indexes] = epoch(EEG.data, tmp_stimstmap(sidx), [-.2 .5].*EEG.srate);
            % remove baseline
            epochdatm = epochdat - mean(epochdat(:,[1:(abs(indexes(1)))],:),2);
            
            % obtain trimmed mean of the EEG epoched data (20%)
            [~,tmdatamatrix(:,:,g),~,~,~,~,~]=limo_trimci(epochdatm);
            
            % Number of trials
            Description(g,1) = sum(sidx);
            % Number of touches
            scount = and(EEG.Aligned.Phone.Blind{1,1}(:,2)>tstamp, EEG.Aligned.Phone.Blind{1,1}(:,2) < (tstamp+[(minu*1000*60)*(EEG.srate/1000)])) ;
            Description(g,2) = sum(scount);
            
            % Typical gap between touches in ms
            if sum(scount)>0
                Description(g,3) = median(diff(EEG.Aligned.Phone.Blind{1,1}(scount,1) ));
            else
                Description(g,3) = minu*60*1000;
            end
            % In(1) or out(0) determined based on the onset sample
            currentTime = onset+(tstamp*[1000/EEG.srate]./1000); % Current time
            
                tmp_nidx = find(strcmp(rawsheet(:,1),tmp_name) ==1);
                Indoor_onset = posixtime(datetime(rawsheet{tmp_nidx,27},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
                Indoor_offset = posixtime(datetime(rawsheet{tmp_nidx,28},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
                Dummy_onset = posixtime(datetime(rawsheet{tmp_nidx,24},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
                Dummy_offset = posixtime(datetime(rawsheet{tmp_nidx,25},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
                
                Description(g,4) =  and(currentTime>Indoor_onset,currentTime<Indoor_offset);
                Description(g,5) =  and(currentTime>Dummy_onset,currentTime<Dummy_offset);

                
                
            
        end
    end
    
    
    
    t = 1:g;
    
    
        
        %% Outdoor LIMO model
        clear d_id; 
        d_id = and([Description(:,4) ==0],[Description(:,5) ==0]); 
        tmp_EEG_all = tmdatamatrix(:,:,d_id);
        tmp_var_all = [sqrt(Description(d_id,2)), t(d_id)'];
        
        if length(tmp_var_all) > 10
            tmp_dir_LIMO = [dat_out,'\',tmp_name,'\RegOutdoor\'];
            mkdir(tmp_dir_LIMO);
            cd(tmp_dir_LIMO);
            LIMO.dir = tmp_dir_LIMO;
            LIMO.Analysis = 'Regression';
            LIMO.Type = 'Channels';
            LIMO.Level = 1;
            LIMO.data.chanlocs = expected_chanlocs;
            LIMO.data.neighbouring_matrix = channeighbstructmat;
            LIMO.bootstrap = 0;
            LIMO.design.bootstrap =0;
            LIMO.design.tfce = 1;
            LIMO.design.fullfactorial = 0;
            LIMO.data.Cont = tmp_var_all';
            save LIMO LIMO
            limo_random_robust(4,tmp_EEG_all,tmp_var_all',[1],0,0);
            close all; clear tmp_var_all tmp_EEG_all
        end
        
        
            clear EEG

    end
    
            
            
        
    
    
