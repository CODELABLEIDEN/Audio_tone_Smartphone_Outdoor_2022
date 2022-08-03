%% Settings 
Dir1 = '\EEGAPPSENSORORGANISED_2018_2019'; 
Dir2 = '\ProcessedData_EEG_Smartphone_WD\EEGAPPSENSORORGANISED_2018_2019'; 

syncfolder(Dir1, Dir2, 0);% Ensure data is in sync before processing
    

ProcessDir = Dir1; % Choose directory to work from 
overwrite = false; % Overwrite ICA files if true
trimdir = false; % To process on distinct machines starting from 100 set to true. 

%% get summary of data collected and then RUN ICA 
Sublist = dir(ProcessDir);
Sublist(ismember( {Sublist.name}, {'.', '..'})) = []; 
Sublist(~[Sublist.isdir]) = [];

% Now go through each dir and reveal: 

for i = 1:length(Sublist)
    try
   load(strcat(Sublist(i).folder,'\',Sublist(i).name,'\Status.mat'))
   EEGstat(i) = sum([eeg_name.EEG])>0;
   EEGPhonestat(i) = sum(and([eeg_name.Phone],[eeg_name.EEG]))>0;
   EEGSensorstat(i) = sum(and([eeg_name.Sensor],[eeg_name.EEG]))>0;
   EEGICA(i) = sum([eeg_name.ICAstatus])>0;
    end
end

Report.NumEEGsub = sum(EEGstat); 
Report.NumEEGPhonesub = sum(EEGPhonestat);
Report.NumEEGSensorsub = sum(EEGSensorstat); 
Report.NimICA = sum(EEGICA); 

disp(Report)

%% Now go through each subject's folder and RUN ICA for the EEG files in Status.mat 
% keep running till two directories have the same size 
check = false; % variable that runs the loop

clear Sublist; 
% Directory listing
Sublist = dir(ProcessDir);
Sublist(ismember( {Sublist.name}, {'.', '..'})) = []; 
Sublist(~[Sublist.isdir]) = [];

if trimdir
    
    Sublist(1:97) = [];
    
end

while ~check 
 
display('%%%%%%%%%%%%%%RUNNING ICA AGAIN AFTER SYNCING DIRECOTRIES%%%%%%%%%%%%');
syncfolder(Dir1, Dir2, 0);% Ensure data is in sync before processing


for i = 1:length(Sublist)
    clear eeg_name; 
    try
   load(strcat(Sublist(i).folder,'\',Sublist(i).name,'\Status.mat'));
    end
   if exist('eeg_name','var') == 1
       for e = 1:length(eeg_name)
           display(strcat('Starting ICA for',Sublist(i).name))
           tic
           % Perform EEG only if EEG indeed exists + overwrite is true
           if and(or(~eeg_name(e).ICAstatus,overwrite),eeg_name(e).EEG)
           EEG = pop_loadset(strcat(eeg_name(e).processed_name,'.set'),strcat(Sublist(i).folder,'\',Sublist(i).name));
           
           % Actual ICA 
           EEG = eeg_checkset( EEG );k = 0; clear chidx; 
           for d = 1:length(EEG.chanlocs)
           if strcmp(EEG.chanlocs(d).type,'EEG');
           k = k+1; 
           chidx(k) = d; 
           end
           end
           EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', chidx);


           % Save the procesed data 
           EEG = pop_saveset( EEG, 'filepath', strcat(Sublist(i).folder,'\',Sublist(i).name),'filename', strcat(eeg_name(e).processed_name,'.set'));    clear EEG; 
           eeg_name(e).ICAstatus = true; 
           save(strcat(Sublist(i).folder,'\',Sublist(i).name, '\Status.mat'), 'eeg_name');
           toc
         
           end
           
       end
       
       
   end
   
   try
       syncfolder(Dir1, Dir2, 0);% Ensure data is in sync before processing

   end
end



check = isequal(DirSize(Dir1), DirSize(Dir2));
end
