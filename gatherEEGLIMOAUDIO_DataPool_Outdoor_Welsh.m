%% For each EEG data, gather the 1st level analysis (for inside only). 
% App category: Non-social (0) or Social (1)
% Interaction Density Continous 


% Directory settings 
data_in = '\Feb_2018_2020_RUSHMA_ProcessedEEG'; 
dat_out = '\Smartphone_AEP_2021\AllSubjectswelch\';

% General LIMO 
load('\Human_electrodes\expected_chanlocs')
load('\EEG\Human_electrodes\channeighbstructmat')

% Load ICA config set by Mark 
load('\Leiden_CODELAB\Paper_specific_codes\AEP_RWSeries_2021\ICinfo_audio_bilat.mat')

% Load RW excel configured by Roderick and Mark 
[~,~,rawsheet] = xlsread('\Leiden_CODELAB\Paper_specific_codes\AEP_RWSeries_2021\DataDetailList_MvdR.xlsx');

ICinfo([5 19 21]) = [];% rejections (5) or data error 
%% Now loop through to perform LIMO regressions
for s = 1:length(ICinfo) % 
    %try
    % Gather times of EEG and trim eeg_name to EEG elevant eeg_name
    [~,tmp_name] = fileparts(fileparts(ICinfo(s).path));
    
    % Get EEG data
    tmp_EEG_all = [];
    tmp_var_all = [];
    
    tmp_dir =  strcat(data_in,'\', tmp_name);
    EEG = pop_loadset(strcat(ICinfo(s).processed_name,'.set'), tmp_dir) ;
    EEG = gettechnincallycleanEEG40hz(EEG);
    
    
    % Onset of EEG in UTC
    dvec_t = EEG.urevent(1).bvtime;
    dvec =[dvec_t{1,1}(1:4),'\' , dvec_t{1,1}(5:6), '\', dvec_t{1,1}(7:8), '\', dvec_t{1,1}(9:10), '\', dvec_t{1,1}(11:12), '\', dvec_t{1,1}(13:14)];
    dvec_ms = str2num(dvec_t{1,1}(15:end))./1000;
    onset = posixtime(datetime(dvec,'InputFormat','yyyy\MM\dd\HH\mm\ss','TimeZone', 'Europe/Amsterdam'))+(dvec_ms/1000);
    
    %% loop through in 5 minute windows to create an spectral matrix ch x spectral x window, and with the properties noted in the window (indoor, outdoor
    minu = 1; % number of min
    tmp_sample = EEG.srate;% samplerate
    tmp_idx = strcmp({EEG.event.type},'M  1');
    tmp_stimstmap = [EEG.event(tmp_idx).latency];
    dur = [tmp_stimstmap(end)-tmp_stimstmap(1)]./(EEG.srate*60*minu);   % number of 1 minutes
    
    
    % setup EEG indices
    tstamp = tmp_stimstmap(1);
    startidx = floor(tstamp);
    endidx = tmp_stimstmap(1)+[(minu*1000*60*(floor(dur)))*(EEG.srate/1000)];

    
    
     tmdatamatrix = []; Description = [];
    Nlength = ([endidx-startidx]);
    
   
    for ch = 1:64
        tic; 
             px_z = zscore(EEG.data(ch, startidx:endidx-1));
             
             for m = 1:(floor(dur)-1)
                tstamp = 1+[(minu*1000*60*m-1)*(EEG.srate/1000)];
                startx = floor(tstamp);
                endx = floor(tstamp+[(minu*1000*60)*(EEG.srate/1000)]);
                px_z_tmp = px_z(startx:endx); 

                % estimate PSD from Cox & Fell (2020) % https://zzz.bwh.harvard.edu/luna/img/cf/matlab.html
                fftWindowLength=5; %in seconds
                fftWindowOverlap=0.5; %fraction of 1.
                [rawPSD,freq_welch]=pwelch(double(px_z_tmp),EEG.srate*fftWindowLength,fftWindowOverlap*EEG.srate*5,[0.1:0.1:40],EEG.srate);
                
                tmdatamatrix(ch,:,m) = rawPSD;
             clear px_z_tmp rawPSD; 

             end

            display(['completed channel...',num2str(ch)]); 
            toc
    end
    
    % Gather descriptors
    
    
    for m = 1:(floor(dur)-1)
        % List of events
        tstamp = tmp_stimstmap(1)+[(minu*1000*60*m-1)*(EEG.srate/1000)];
        startidx = floor(tstamp);
        endidx = floor(tstamp+[(minu*1000*60)*(EEG.srate/1000)]);
            
        % Number of touches
        scount = and(EEG.Aligned.Phone.Blind{1,1}(:,2)>tstamp, EEG.Aligned.Phone.Blind{1,1}(:,2) < (tstamp+[(minu*1000*60)*(EEG.srate/1000)])) ;
        Description(m,2) = sum(scount);
        
        % Typical gap between touches in ms
        if sum(scount)>0
            Description(m,3) = median(diff(EEG.Aligned.Phone.Blind{1,1}(scount,1) ));
        else
            Description(m,3) = minu*60*1000;
        end
        % In(1) or out(0) determined based on the onset sample
        currentTime = onset+(tstamp*[1000/EEG.srate]./1000); % Current time
        %if isfield(EEG, 'WatchData')
            %Description(m,4) =  and(currentTime>EEG.WatchData{1,1}.Indoor_onset,currentTime<EEG.WatchData{1,1}.Indoor_offset);
        %else
        tmp_nidx = max(find(strcmp(rawsheet(:,1),tmp_name) ==1)); % latest index used
        Indoor_onset = posixtime(datetime(rawsheet{tmp_nidx,27},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
        Indoor_offset = posixtime(datetime(rawsheet{tmp_nidx,28},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
        Dummy_onset = posixtime(datetime(rawsheet{tmp_nidx,24},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
        Dummy_offset = posixtime(datetime(rawsheet{tmp_nidx,25},'ConvertFrom', 'excel', 'TimeZone', 'Europe/Amsterdam'));
        
        Description(m,4) =  and(currentTime>Indoor_onset,currentTime<Indoor_offset);
        Description(m,5) =  and(currentTime>Dummy_onset,currentTime<Dummy_offset);
            
            
            
        %end
    end


    g = 1:m;



%% Outdoor LIMO model

d_id = and([Description(:,4) ==0],[Description(:,5) ==0]);
tmp_EEG_all = tmdatamatrix(:,:,d_id);
tmp_var_all = [sqrt(Description(d_id,2)), g(d_id)'];

if length(tmp_var_all) > 10 && (sum (tmp_var_all(:,1)) > 0)
    tmp_dir_LIMO = [dat_out,'\',tmp_name,'\RegOutdoor\'];
    mkdir(tmp_dir_LIMO);
    cd(tmp_dir_LIMO);
    LIMO.dir = tmp_dir_LIMO;
    LIMO.Analysis = 'Regression';
    LIMO.Type = 'Channels';
    LIMO.Level = 1;
    LIMO.Freq = freq_welch;
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
%     MsgBoxH = findall(0,'Type','figure','Name','TitleOfMsgBox');
%     close(MsgBoxH);






    
