%% Use regression Betas towards one sample LIMO t-test 
% Manually removed RW08 due to odd signals (x10 fold at scale). 

% Directory 
Processed_dir = [BetaSource,'\','LIMO_onesample'];

% General LIMO
load('\Human_electrodes\expected_chanlocs')
load('\Human_electrodes\channeighbstructmat')
nboot = 1000;
tfce = 0;

sublist = dir([BetaSource,'\','RW*']);
% Gather data from all subjects with Indoor folder
k = 0;
for s = 1:length(sublist)
    if exist([BetaSource,'\',sublist(s).name,'\RegOutdoor\Betas.mat']) ==2
        load([BetaSource,'\',sublist(s).name,'\RegOutdoor\Betas.mat']);
        load([BetaSource,'\',sublist(s).name,'\RegOutdoor\LIMO.mat']);
        if std(LIMO.data.Cont(1,:)) > 1;
            k = k+1;
            Par_yr_raw(:,:,k) = squeeze(Betas(:,:,1));
            SampleN(k) = length(LIMO.design.X);
        end
    end
end

    
    % RUN LIMO 
             
    
    
      tmp_dir_LIMO = [Processed_dir,'\Outdoor'];
      mkdir(tmp_dir_LIMO);
      cd(tmp_dir_LIMO);
      LIMO.dir = tmp_dir_LIMO;
      LIMO.Analysis = 'One sample t-test';
      LIMO.Type = 'Channels';
      LIMO.Level = 2;
      LIMO.data.chanlocs = expected_chanlocs;
      LIMO.data.neighbouring_matrix = channeighbstructmat;
      LIMO.bootstrap = 0;
      LIMO.design.bootstrap =0;
      LIMO.design.tfce = tfce;
      LIMO.design.fullfactorial = 0;
      LIMO.SampleN= SampleN; 
      
      save LIMO LIMO
      save Yr Par_yr_raw
      limo_random_robust(1,Par_yr_raw,[1],nboot,tfce);
      
      % H0 mask 
      load([LIMO.dir, '\H0\H0_one_sample_ttest_parameter_1.mat']);
      load([LIMO.dir, '\one_sample_ttest_parameter_1.mat']);

      [mask] = limo_cluster_correction([one_sample(1:60,:,4)].^2,one_sample(1:60,:,5), [squeeze(H0_one_sample(1:60,:,1,:))].^2, squeeze(H0_one_sample(1:60,:,2,:)),LIMO.data.neighbouring_matrix(1:60,1:60),2,0.05);
      
      save([LIMO.dir,'\mask.mat'], 'mask'); 

      close all;
      clear LIMO tmp*

    
