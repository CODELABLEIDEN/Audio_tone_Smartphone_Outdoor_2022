function getdataorganised(sub,inputpath,outputpath,phonepath,overwrite)
% Usage getdataorganised(sub,inputpath,outputpath,overwrite)
% Organizes the data for a given subject with the following outputs, 
% Applicable to all data collected through 2018 - 2020 at CODELAB
% INPUT
% sub : subject name as in 'AG01'
% inputpath: folder where the subject data is stored 
% Phone path: folder where processed phone data is stored 
% Overwrite: If set true then overweite pre-existing data(default true)
%
% OUTPUT 
% For each subject coded in traditonal form like AG01 creates a new EEG
% folder with EEG.set for each session
% Attys or NIH daq data is added into the .set file 
% Phone data is added as well 
% Contains a companion summary file containing the following information per file
% Subject level ##Curfew##Age(at first EEG session)##Gender##FilterSettings
% Session level ##Pre-Curfew##Post-Curfew
% ##{Triggers,Num of Triggers}##ICA status##ONSET UTC##OFFSET UTC##NIH DAQ
% Status##Attys Status##Experiment Type: RT, Passive Tactile, Phone, Active
% Tactile, Audio, Box with FS, Indoor vs. Outdoor 
%
% Performs the following EEG steps 
%   1. Removes bad channels based on poor impedence
%   2. Filter (0.1 to 70 Hz)
%   3. NOT DONE. RunICA > Use Run_ICA.m for that
%   4. NOT DONE. wICA to remove ICA, Blinkremoval based on ICA OR
%      interpolation of missing channels and then average referencing 
% Performs the following Phone steps
%   1. Trim phone data to the session overlap only 
%
% Arko Ghosh, Leiden University, 27th June 2020% 
% Requires xlxs file paths to be manually set(stored in the paper specific one-drive-folder):
% DataDetailList_MvdR
% MASS_Subject_list


%% Set output path 
outputdir = strcat(outputpath,'\',sub);

%% Set input dir 
inputdir = strcat(inputpath,'\',sub);

%% Manually set excel paths 
Mass_sheet = '\\CODELABGAMMA\Users\aghos\OneDrive - fsw.leidenuniv.nl\Leiden_CODELAB\Paper_specific_codes\MassOrganizationEEG_Phone_Audio_Tactile\MASS_Subject_list.xlsx';
RW_sheet = '\\CODELABGAMMA\Users\aghos\OneDrive - fsw.leidenuniv.nl\Leiden_CODELAB\Paper_specific_codes\MassOrganizationEEG_Phone_Audio_Tactile\DataDetailList_MvdR.xlsx';

%% if Status file exists load that or create new 
if and(exist(strcat(outputdir,'\','Status.mat'),'file') == 2,~overwrite)
    load(strcat(outputdir,'\','Status.mat'));
else 
    

% Find the Subject (sub) folder 
% List all files that may contain data 
dir_files = dir(inputpath); 
dir_files=dir_files(~ismember({dir_files.name},{'.','..'}));
% remove all files excep the sub
if exist('sub')
idx_d = ismember({dir_files.name}, sub);
dir_files(~idx_d) = [];
display('Running single subject mode'); 
end

% Decide to continue 
if isempty(dir_files)
    display('No such subject, no processed data created')
    return
end

% Find the corresponding EEG folder
dirs = regexp(genpath(strcat(dir_files.folder,'\',dir_files.name)),['[^;]*'],'match');

% filter out EEG folders in the subject data 

for d = 1:length(dirs)
    s = ([dirs{1,d}]);
    Find_EEG{1,1} = strfind(s, 'EEG');
    Find_EEG{2,1} = strfind(s, 'Eeg');
    Find_EEG{3,1} = strfind(s, 'eeg');

    if [Find_EEG{:,1}] > 5 
    Idx_EEG(d) = true;  
    else
            Idx_EEG(d) = false;  
    end
end

% Decide to continue 
if sum(Idx_EEG)<1
    display('No EEG data, no processed data created')
    return
end

EEG_folder = dirs(Idx_EEG);

% List of VHDR files anf find processed names 
k = 0; 
for f = 1:length(EEG_folder);
eegfiles = dir(strcat(EEG_folder{1,f},'\','*.eeg'));
for e = 1:length(eegfiles)
    k = k+1; 
eeg_name(k).name = strrep(eegfiles(e).name, '.eeg', '.vhdr'); 
eeg_name(k).path = EEG_folder{1,f}; 

tmp = strrep(eegfiles(e).name, '.eeg', '.vmrk'); 

[eeg_name(k).start eeg_name(k).stop] = getEEGonsetoffset(eeg_name(k).name, eeg_name(k).path); % Start stop times in UTC seconds 

eeg_name(k).processed_name  = datestr(datetime(eeg_name(k).start, 'ConvertFrom', 'Posixtime'), 'HH_MM_dd_mm_yy');
% set the process status for the different sensors 
eeg_name(k).EEG = false; 
eeg_name(k).Sensor = false; 
eeg_name(k).Phone = false; 
eeg_name(k).Watch = false; 
eeg_name(k).ICAstatus = false; 
eeg_name(k).CurfewExp = false; 

end
end

% Save intermediate file eeg_name to show and log work in progress
mkdir(outputdir); 
statusname = strcat(outputdir,'\Status.mat');
save(statusname,'eeg_name')
end



%% Gather BS/Sensor data 
    % use standard DAQ import 
    try
    [BS,bs_name,attysflag] = getBSdata(inputdir);
	if attysflag
    [BS,bs_name] = getBSaudioattysdata(inputdir);
    end 
    catch
    BS = [];bs_name =[];
    display('No sensor data found');
    end
    
%% Get QA ID for the identified subject 
    [~,~, Curfew, ~, ~, ParticipantID] = getdemoinfo({sub}, Mass_sheet, [], 2);
    
%% Gather phone data
    % Gather phone data 
    try
    Data = gettapdata(ParticipantID{1,1},phonepath);
    % Trim it to EEG range to save memory 
    overallEEGonset = min([eeg_name.start]).*1000; %convert to UTC ms
    overallEEGoffset = max([eeg_name.stop]).*1000; % convert to UTC ms
    
    Timestamp = [Data{1,1}.SUBJECT.tap.timestamp]; % original timestamps
    IdxKeep = and(Timestamp>overallEEGonset,Timestamp<overallEEGoffset);
    Data{1,1}.SUBJECT.tap(~IdxKeep,:)=[];
    catch
    Data{1,1}.SUBJECT.tap =[];
    display('No phone data found');    
    end
    
%% Gather watch data 
    % This is only worth doing for RW data 
    % It is possible that the Watch Data was not recorded 
    try 
    if strfind(sub,'RW')>0
    WatchData = getGendata(inputdir);
    % Use the LUX values to mark indoor vs. outdoor times 
    Idx = findchangepts(movsum(abs((WatchData{1,1}.light)),100*60));% smooth the data at min resolution
    % Determine Indoor vs. OutDoor
    p1 = std(WatchData{1,1}.light(1:Idx));
    p2 = std(WatchData{1,1}.light(Idx:end));
    if p1>p2 % higher lux outside 
        WatchData{1,1}.Indoor_onset = WatchData{1,1}.time(Idx); WatchData{1,1}.Indoor_offset = WatchData{1,1}.time(end); 
        WatchData{1,1}.Outdoor_onset = WatchData{1,1}.time(1); WatchData{1,1}.Outdoor_offset = WatchData{1,1}.time(Idx); 
    else
        WatchData{1,1}.Indoor_onset = WatchData{1,1}.time(1); WatchData{1,1}.Indoor_offset = WatchData{1,1}.time(Idx); 
        WatchData{1,1}.Outdoor_onset = WatchData{1,1}.time(Idx); WatchData{1,1}.Outdoor_offset = WatchData{1,1}.time(end);
    end
    
    end
    catch
       display(strcat('Perahps no watch data:::::',inputdir)); 
    end

%% Now go through EEG steps [Depends on UTC time conversations all based on ]
% 1. If mobile EEG exists process only that if it meets the criteria incl. of MB > non mobile 
% 2. Find corresponding BendSensor files and integrate into the .set file
% 3. Find the corresponding Phone Data and integrate into the .set file 
% 4. Get a summary of all triggers and store it in eeg_name 
% 5. Insert Channel Map information 
% 5. Delete poor impedence channels (10kOhm, 10.1097/WNP.0000000000000308)
% 6. Filter data between 0.1 Hz to 70 Hz (no notch)
% 7. Set ICA status
% 8. Get Excel data if ti exists and place into EEG.set file 
% 9. Identify CURFEW exp [0/1]; Group[0/1]; PrePost[0/1]
% 10. Save data 

%% 1. Are there overlapping EEG files? If yes, choose the larger one...

for e = 1:length(eeg_name)
interval1 = [eeg_name(e).start:eeg_name(e).stop];% *****
vidx = [1:length(eeg_name)];
pidx = find(vidx==e); % find the file index so it is not compared to self 
vidx(pidx) =[]; k = 0; larger =[];
for f = vidx
k = k+1; 
% are these intervals overlapping with anover? 
interval2 = [eeg_name(f).start:eeg_name(f).stop];% *****
check = length(intersect(floor(interval1), floor(interval2)))>1;
larger =[];
    % if it overlaps, then is this intevallarger than e? 
    if and (check,abs(diff([interval1(1) interval1(end)]))<abs(diff([interval2(1) interval2(end)])))
    larger(k) = true;
    else
    larger(k) = false;
    end
end
% if any larger val. exists then mark it for removal
if sum(larger)>0
    Idx_rej(e) = true;
elseif isempty(larger)
    Idx_rej(e) = false;
else
    Idx_rej(e) = false;
end
end


idexeeg = 1:length(eeg_name); % indices of eegfiles gathered here 
idexeeg(Idx_rej) = []; % remove from consideration 

%% 2-10 Now, perform the corresponding analysis 
for s = idexeeg

% load EEG data 
[EEG, com] = pop_loadbv(eeg_name(s).path,eeg_name(s).name);
eeg_name(s).urevent = EEG.urevent; 

% Apply channel locations 
chanfile = fileparts(which('expected_chanlocs.mat'));
load(strcat(chanfile,'\','expected_chanlocs.mat'));

for i = 1:64
    EEG.chanlocs(i).X = expected_chanlocs(i).X;
    EEG.chanlocs(i).Y = expected_chanlocs(i).Y;
    EEG.chanlocs(i).Z = expected_chanlocs(i).Z;
    EEG.chanlocs(i).labels = expected_chanlocs(i).labels;
    EEG.chanlocs(i).theta = expected_chanlocs(i).theta;
    EEG.chanlocs(i).radius = expected_chanlocs(i).radius;
    EEG.chanlocs(i).sph_theta= expected_chanlocs(i).sph_theta;
    EEG.chanlocs(i).sph_phi= expected_chanlocs(i).sph_phi;
    EEG.chanlocs(i).sph_radius = expected_chanlocs(i).sph_radius;
    EEG.chanlocs(i).type = 'EEG';
    EEG.chanlocs(i).urchan = i; 
end

% Identify and remove bad channels 
    [eeg_name(s).info.impedence, eeg_name(s).info.imptimes, eeg_name(s).info.fs, eeg_name(s).info.chanlabel] =  getvhdrvals (eeg_name(s).name, eeg_name(s).path);
    % Index of bad channels 
    EEG.Orignalchanlocs = EEG.chanlocs;
    if ~isempty(eeg_name(s).info.impedence)
    bad_chan = find(nanmedian(eeg_name(s).info.impedence(1:64),2)>10); % remove channel with impedences higher than 10kOhm
    EEG = pop_select(EEG, 'nochannel', bad_chan); % Trimmed data set 
    EEG.BadChannels = bad_chan; 
    eeg_name(s).info.BadChannels = bad_chan; 
    end
% Filter EEG data 
    
    % Perform filtering

    [EEG] = pop_eegfiltnew(EEG, 0.1, 70);

    
%   Set EEG status 
    eeg_name(s).EEG = true; 
    eeg_name(s).CurfewExp = Curfew;

%%   If BS data exists and overlaps then insert it into EEG
     if ~isempty(BS)
         interval1 = [eeg_name(s).start:eeg_name(s).stop]; j = 0; % set index for storage 
         for b = 1:length(bs_name)
         interval2 = [bs_name(1,b).start:bs_name(1,b).stop];
         overlap = intersect(floor(interval1),floor(interval2));
         if length(overlap)>1
             j = j+1;
             EEG.BS{1,j}.Data = BS{1,b}; 
             EEG.BS{1,j}.info = bs_name(1,b);
             eeg_name(s).Sensor = true; 
         end
         end
     end

%%   If Phone data exists and overlaps then insert into EEG
     if ~isempty(Data{1,1}.SUBJECT.tap)
        EEG.PhoneData = Data;  
        eeg_name(s).Phone = true; 
     end

%%   If WatchData exists and overlaps insert into EEG
     if exist('WatchData','var') == 1
         WatchOnset = WatchData{1,1}.time(1);
         Watchoffset = WatchData{1,1}.time(end);
         
         interval1 = [eeg_name(s).start:eeg_name(s).stop];
         interval2 = [WatchOnset:Watchoffset];
         overlap = intersect(floor(interval1),floor(interval2)); 
         if length(overlap)>1
         checkwatch = true;
         else
         checkwatch = false; 
         end
         else
         checkwatch = false; 
     end

     if checkwatch
             EEG.WatchData = WatchData; % Insrt dats into EEG.set 
             eeg_name(s).Watch = true; 
     end

     
     %% Save the Processed data and the companion file 
    save(strcat(outputdir, '\Status.mat'), 'eeg_name'); 
    EEG = pop_saveset( EEG, 'filepath', outputdir,'filename', eeg_name(s).processed_name);    clear EEG   

     
end

end


