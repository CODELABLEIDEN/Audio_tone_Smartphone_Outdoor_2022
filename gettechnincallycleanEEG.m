function OUTEEG = gettechnincallycleanEEG(INEEG)
% OUTEEG = gettechnincallycleanEEG(INEEG)
% assume 1.ICA performed, 2. Bad channels have been marked, 3. channel
% locations already placed 4. Electrode locations present 
% Perfrorms the following step
% 1. Remove blinks according to ICA
% 2. Interpolate missing channels 
% 3. Re-reference data to average channel 

% Remove non_EEG data 
idxr = find(strcmp({INEEG.chanlocs.labels}, 'x_dir') == true);% are blink channels present 
if idxr>0
   INEEG.data(idxr:idxr+2,:) = []; 
   INEEG.chanlocs(idxr:idxr+2) = []; 
   INEEG.nbchan = idxr-1;
end

% Remove blinks based on ICA 

idx = find(strcmp({INEEG.chanlocs.labels}, 'E64')|strcmp({INEEG.chanlocs.labels}, 'E63') == true);% are blink channels present 

rej = []; 
for i = 1:length(idx)
INEEG.icaquant{i} = icablinkmetrics(INEEG, 'ArtifactChannel', INEEG.data(idx(i),:), 'Alpha', 0.001, 'VisualizeData', 'False');
rej = [rej INEEG.icaquant{1,i}.identifiedcomponents]; 
end

Rej_Comp = unique([rej]);
INEEG = pop_subcomp(INEEG,Rej_Comp,0);

% interpolate missing channels 
% remove channel locations which are not EEG
INEEG.Orignalchanlocs(find(cellfun(@isempty,{INEEG.Orignalchanlocs.type})))=[];
INEEG.chanlocs(find(cellfun(@isempty,{INEEG.chanlocs.type})))=[];


INEEG = pop_interp(INEEG,INEEG.Orignalchanlocs,'spherical');

% Final Filter 
[INEEG] = pop_eegfiltnew(INEEG, [], 45);

% Re-reference to average channel 

OUTEEG = pop_reref (INEEG, [1:62], 'keepref', 'on'); OUTEEG = eeg_checkset(OUTEEG); % use this re-ref to avg reference.

