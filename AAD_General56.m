%% Initial parameters

clear all

nofsubj = 4;
lag = [100 400];                                                    
LP = 2;
HP = 9;
%LP = [0, 1, 2, 3, 4];
%HP = [6,7, 8, 9, 10];
%NSR = 20;
NewSamplingRate = [10 15 20 25 30 35 40 45 50 80]; % 200 500]; % New sampling rate; 200 and 500 take to much computing time
trll = [10 20 30 40 50 60 5 15]; % in seconds
wantplots = false; % want plots?

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?

 subband = false; 
 spacing = 1.5; % OPM: How is this defined? Ask Neetha...
 gammatonefreqs = erbspacebw(150,4000,spacing); % gammatone filter centerfrequencies % OPM: What is all off this? How are the inputs defined? Ask Neetha...
 gammatonebetamul = spacing*ones(size(gammatonefreqs)); % multiplier for gammatone filter bandwidths

intermediate = true; % Neetha does this, check what the effect is
intSamp = [20 70 120 160 200 300 400 500]; % Why 120?
cut_off = 4000; %cut-off frequency of headphones. OPM: Find what this is in our experiment.

remove_DC = false; % Neetha does this, check what the effect is
%% Delete if wished

delete(''); % delete all file in current directory - necessary to prevent loop problems

%% Get the covariances
for LP = 2
    for HP = 9
for nsr = 3:3
    NSR = NewSamplingRate(nsr)
for isr = 8:8
    
    intermediateSampleRate = intSamp(isr)
for TRL = 1:1 % Loop for trial length
    
    triallength = trll(TRL)
    covar = []; X=[]; ALall = []; ARall = [];
    EE = cell(1, nofsubj); Attended = cell(1, nofsubj); Unattended = cell(1, nofsubj);
    gentel4 = 0; sessionlength = [];
for S = 1:nofsubj % Subject loop
    gentel5 = 0; gentel6 = 0;
    ses = 0
    
    for LR = 1:2 % Left - Right loop
        if LR == 1
            left = true;
        else
            left = false;
        end
        
        gentel = 1; gentel2 = 1;  gentel3 = 1;
        for session = 1:4 % 4 right ear session - Train on each one separate
            data = [];
            [data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            data = cat(2,data2(:,1),data3(:,1));
            
            ses = ses+1
            
            if S ==1
                if left 
                    if session == 1 % small correction mistake in pilot data. % ??? Was dit alleen bij het linker oor? 
                        data = cat(2,data(:,2),data(:,1)); % Left and Right audio are switched
                    end
                    EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_L.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                else
                    EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_R.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                end
                regions = [];
            elseif S >1
                    if left
                         EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_L.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                    else
                        EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_R.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                    end
                    %[EEG, regions] = pop_rejcont(EEG, 'elecrange',[1:24] ,'freqlimit',[20 40] ,'threshold',8,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');
                    %EEG = pop_select( EEG,'nopoint',regions); % Test to remove artefact parts of the EEG/Audio
                    regions = [];
            end
            
            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
            [pks,locs] = findpeaks(double(EEG.data(28,:))); % Locate Triggers
 
            if session == 4 % Auto add trigger after 6min30s because incongruent length of audio stories so missing trigger.
                locs(2) = locs(1) + ((395)*EEGsampleRate);
            else
                locs(2) = locs(1) + ((size(data,1)/fs)*EEGsampleRate);
            end
            EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            % EEGaad = double(EEGaad); % Dit stond nog in vorig script;
           
            % Remove DC component (Done in Neethas script in 'biopil_filter') 
            % OPM: Dit wel of niet doen???
            % OPM: Matlab zegt: Do not use filtfilt with differentiator and Hilbert FIR filters, because the operation of these filters depends heavily on their phase response.
            if remove_DC
                [b, a] = butter(1, 0.5/(EEGsampleRate/2), 'high'); % Klopt dit???
                for i = 1: 24     
                    EEGaad(i,:) = filtfilt(b,a,double(EEGaad(i,:)));
                end
            end
            
            AudioL = data(:,1); AudioR = data(:,2);
            data = [];
            
            if subband
                % OPM in the script of Neetha in aad_envelopes the audio is
                % intermediately resampled to an intermediatefs of 8000 Hz. Why? Do we
                % need this? Why not immediately to 500 like Rob suggests?
                % OPM: dit wordt niet alleen in geval van subbands
                % gedaan...
                AudioL = resample(AudioL, cut_off*2,fs);
                AudioR = resample(AudioR, cut_off*2, fs);
                fs = cut_off*2;

                % from amtoolbox>joergensen2011.m
                g = gammatonefir(gammatonefreqs,fs,[],gammatonebetamul,'real'); % create real, FIR gammatone filters. SEE aad_plot_gammatone
                
                AudioL = real(ufilterbank(AudioL,g,1)); %audio is filtered along columns. wat is difference with ufilterbankz ? Real or ABS ? See test_gammatone_differences.m
                AudioR = real(ufilterbank(AudioR,g,1)); %audio is filtered along columns. wat is difference with ufilterbankz ? Real or ABS ? See test_gammatone_differences.m

                % AudioL = reshape(AudioL,size(AudioL,1),[]); % make 2D for
                % easier handling % OPM: Niet nodig in dit script??? -> Ask Neetha
            end
            
            % Rob resamples before making the envelope
            if intermediate
            else           
                % Rob's method: Rob downsamples the audio before making the
                % envelope. Neetha does this after making the envelope.
                % Does this make a difference? 
                AudioL = resample(AudioL,EEGsampleRate,fs); % Downsampled audio per ear
                AudioR = resample(AudioR,EEGsampleRate,fs);
                fs = EEGsampleRate;  
            end
            
            % ENVELOPES
            % Powerlaw
            if strcmpi(envelopemethod, 'powerlaw')
                envelopeL = abs(AudioL).^power;
                envelopeR = abs(AudioR).^power;
            % Hilbert    
            elseif strcmpi(envelopemethode, 'hilbert')
                envelopeL = abs(hilbert(AudioL)); % Envelope Hilbert Transform
                envelopeR = abs(hilbert(AudioR));
            end
            
            if intermediate
            % Intermediary Downsample envelope & EEG before applying the more
            % strict bpfilters (to save computation) % OPM: Rob does not do this,
            % but he downsamples the audio to 500... What is the difference? 
                envelopeL = resample(envelopeL,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                envelopeR = resample(envelopeR,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                EEGaad = (resample((double(EEGaad))', intermediateSampleRate, EEGsampleRate))';
                EEGsampleRate = intermediateSampleRate;
                fs = intermediateSampleRate;
            end                 
            
            New = [];
            % FILTERING
            
                        % Neetha does it this way: 
%              %bpFilter the envelope in the same way as the EEG.
%     bpFilter = aad_construct_bpfilter(params);
%     envelope = filtfilt(bpFilter.numerator,1,envelope);
            
            for ch = 1:24 % EEG
                [EEGaad(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                New(ch,:) = resample(double(EEGaad(ch,:)),NSR,EEGsampleRate);
            end
            EEGaad = New;
            EEGsampleRate = NSR;
                        
            if subband
                nofsubbands = length(gammatonefreqs);
                for sb = 1:nofsubbands
                    [envelopeL(:,sb)] = Rbp(LP, HP, 6, fs, envelopeL(:,sb)); % Band pass % OPM: Hoe ken je de filter order?
                    [envelopeR(:,sb)] = Rbp(LP, HP, 6, fs, envelopeR(:,sb)); %
                end
                % Sum of subbands so that there only is one envelope again
                envelopeL = (sum(envelopeL,2))';
                envelopeR = (sum(envelopeR,2))';
            else
            %[envelopeL(:)] = Rbp(LP, HP, 6, fs, envelopeL(:,1)); % Band pass
            %[envelopeR(:)] = Rbp(LP, HP, 6, fs, envelopeR(:,1)); % Band pass
            
            [envelopeL(:)] = Rbp([], HP, 6, fs, envelopeL(:,1)); % Band pass
            [envelopeR(:)] = Rbp([], HP, 6, fs, envelopeR(:,1)); % Band pass
            end
            
            envelopeL =  resample(double(envelopeL),NSR,fs);
            envelopeR =  resample(double(envelopeR),NSR,fs);
            
            fs = NSR;
            

            
            % CREATE TRIALS
            window = triallength*NSR; % seconds * Sampling-rate
            [trs] = mod(size(EEGaad,2),window);
            trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
            sessionlength(S, ses) = trn; 
           

            % EEG
            s1 = 0; E = [];
            for i = 1:trn
                for ch = 1:24
                    E(i,ch,:) = EEGaad(ch,s1+1:s1+window);
                end
                s1 = s1+window;
            end
            
            % AUDIO
            s1 = 0; AL = []; AR = [];
            for i = 1:trn
                gentel5 = gentel5+1;
                AL(i,:) = envelopeL(s1+1:s1+window);
                AR(i,:) = envelopeR(s1+1:s1+window);
                if left
                    Attended{S}(gentel5,:) = AL(i,:);
                    Unattended{S}(gentel5,:) = AR(i,:);
                else
                    Attended{S}(gentel5,:) = AR(i,:);
                    Unattended{S}(gentel5,:) = AL(i,:);
                end
                s1 =s1+window;
                gentel4 = gentel4+1;             
                
            end

            % Lags
            test = {}; Rxx = {}; Rxy_left = {}; Rxy_right = {}; x = [];
            for trials = 1:trn
                
                gentel6 = gentel6 + 1
                x=squeeze(E(trials,:,:))';
                start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                nofsamples = size(x,1);
                [start,fin] = deal(-fin,-start);
                [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                
                EE{S} = cat(3,EE{S},X);
                
                % Covariance Matrices
                Rxx{trials} = (X'*X)/nofsamples;
                
                Rxy_left{trials} = (X'*squeeze(AL(trials,:))')/nofsamples;
                Rxy_right{trials} = (X'*squeeze(AR(trials,:))')/nofsamples;
                
                % Ryy_left{trials} = ((AL(1,:))*squeeze(AL(trials,:))')/nofsamples;
                % Ryy_right{trials} = ((AR(1,:))*squeeze(AR(trials,:))')/nofsamples;
                
                
                covar.Rxx_all(S, gentel6,:,:) = Rxx{trials};
                if left
                    covar.Rxy_at_all(S, gentel6,:) = Rxy_left{trials};
                else
                    covar.Rxy_at_all(S, gentel6,:) = Rxy_right{trials};
                end
            end
           
        end

    end
 
end
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')

%save(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'], 'EE', 'covar', 'Attended', 'Unattended', 'sessionlength') % , '-v7.3') Vanaf NSR > 50 nodig
save(['TRL_' num2str(TRL) '_ISR' num2str(isr) '.mat'], 'EE', 'covar', 'Attended', 'Unattended', 'sessionlength') % , '-v7.3') Vanaf NSR > 50 nodig

clear EE Attended Unattended

end
end
end
    end
end


%
%% Subject accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
Subject_Accuracy = [];
nofsubj = 4;
%for LP = [0 1 2 3 4]
    %for HP = [6 7 8 9 10]
for nsr = 3:3
for isr = 8:8
for TRL = 1:1
    
    %load(['LP' num2str(LP) '_HP' num2str(HP) '_TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    %load(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    load(['TRL_' num2str(TRL) '_ISR' num2str(isr) '.mat'])

    for S = 1:nofsubj
        
        noftrials = size(Attended{S},1);
        train = (1:noftrials);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            trials2 = train(1:end-1);
        
            RXX = squeeze(mean(covar.Rxx_all(S,trials2,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,trials2,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE{S}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{S}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{S}(train(end),:)'); % Unattended
            train = [train(end) train(1:end-1)];
        end
        difference(S,:) = CA - CUA;
        
        Subject_Accuracy(TRL,isr, nsr,S) = (sum(CA-CUA> 0)/(noftrials))*100; % Final accuracy on all trials
        CA = []; CUA = [];
        
        noftrlses = noftrials/8;
        for i = 1:8
        Session_Accuracy(S,i) = sum(difference(S,noftrlses*(i-1)+1:noftrlses*i)>0)/noftrlses*100;
        end
        

    end

end
end
end
Subject_Accuracy
Session_Accuracy
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
        save(['Result_sessions_again.mat'], 'Subject_Accuracy', 'difference', 'Session_Accuracy')

        %save(['Result_LP' num2str(LP) '_HP' num2str(HP) '.mat'], 'Subject_Accuracy')
%    end
%end

%%
%% Session accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
Session_Accuracy = [];
nofsubj = 4;
%for LP = [0 1 2 3 4]
    %for HP = [6 7 8 9 10]
for nsr = 3:3
for isr = 8:8
for TRL = 1:1
    
    %load(['LP' num2str(LP) '_HP' num2str(HP) '_TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    %load(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    load(['TRL_' num2str(TRL) '_ISR' num2str(isr) '.mat'])

    for S = 1:nofsubj
        
        noftrials = size(Attended{S},1)/8;
        for ses = 1:8
        train = (noftrials*(ses-1)+1:noftrials*ses);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            trials2 = train(1:end-1);
        
            RXX = squeeze(mean(covar.Rxx_all(S,trials2,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,trials2,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE{S}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{S}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{S}(train(end),:)'); % Unattended
            train = [train(end) train(1:end-1)];
        end
        difference(S,ses,:) = CA - CUA;
        Session_Accuracy(TRL,S,ses) = (sum(CA-CUA> 0)/(noftrials))*100; % Final accuracy on all trials

        end
               
        CA = []; CUA = [];
                

    end

end
end
end
Session_Accuracy
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
        save(['Session_Accuracy.mat'],'difference', 'Session_Accuracy')

        %save(['Result_LP' num2str(LP) '_HP' num2str(HP) '.mat'], 'Subject_Accuracy')
%    end
%end

%% Figure of the different LP and HP combinations
% All
figure
for LP = [0 1 2 3 4]
    for HP = [6 7 8 9 10]
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\investigate HL and LP values')
        load(['Result_LP' num2str(LP) '_HP' num2str(HP) '.mat'])
        
        Accuracy = squeeze(Subject_Accuracy(:, 8,3,:));
        plot(10:10:60, Accuracy)
        hold on
        Accuracy = []; Subject_Accuracy =[];
    end
end
title('Investigate LP and HP; all')
xlabel('Trial length')
ylabel('Accuracy')
legend('06', '07', '08', '09', '010','16', '17', '18', '19', '110','26', '27', '28', '29', '210','36', '37', '38', '39', '310','46', '47', '48', '49', '410')
hold off
% a figure for each LP value (because the above figure might be too
% overwhelming...
for LP = [0 1 2 3 4]
    figure
    for HP = [6 7 8 9 10]
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\investigate HL and LP values')
        load(['Result_LP' num2str(LP) '_HP' num2str(HP) '.mat'])
        
        Accuracy = squeeze(Subject_Accuracy(:, 8,3,:));
        plot(10:10:60, Accuracy)
        hold on
        Accuracy = []; Subject_Accuracy =[];
    end
    title(['Investigate LP and HP; LP = ' num2str(LP)])
xlabel('Trial length')
ylabel('Accuracy')
legend('6', '7', '8', '9', '10')
hold off
end
    
%% All other subjects and all other trials
cd('C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\All_data\New bp filter\intermediate')

isr = 8; nsr = 3;

% OPM: changed to not include the data of subject 2
for TRL = 1: length(trll)
    
    load(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    subj = [1, 3, 4];
    trials3 = []; trials2 = [];
   
    for S = 1:3
 
        noftrials = size(Attended{subj(end)}, 1);
        train = (1:noftrials);
        trials3 = subj(1:end-1);
        CA = []; CUA = [];
        
        for HH = 1:noftrials % repeat for all trials - leave one out
        
            trials2 = train(1:end-1);
            
            sxx = squeeze(mean(covar.Rxx_all(trials3,trials2,:,:),1)); % plain averaging of cov matrices
            RXX = squeeze(mean(sxx,1));
           
            
            sxy = squeeze(mean(covar.Rxy_at_all(trials3,trials2,:),1));
            RXY_ATT = squeeze(mean(sxy,1))';
            
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE{subj(end)}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{subj(end)}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{subj(end)}(train(end),:)'); % Unattended
            
            train = [train(end) train(1:end-1)];
        end
        
        General_Accuracy(TRL,subj(end)) = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials
        
%         if wantplots
%             figure
%             subplot(1,3,1)
%             WW = reshape(avg_dec,abs(start - fin)+1,24);
%             WWm = mean(WW,1);
%             topoplot(WWm,'C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\elec_24ch.elp','style','fill','electrodes','on');
%             laggs = mean(WW,2);
%             subplot(1,3,2)
%             plot(laggs)
%             set(gca,'XTick',[1 2.33 4.67 7])
%             set(gca,'XTickLabel',{'1' '100' '200' '300'})
%             xlabel('Time lagg (ms)')
%             [q8 w8] = max(laggs)
%             WWm = mean(WW(w8,:),1);
%             subplot(1,3,3)
%             topoplot(WWm,'C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\elec_24ch.elp','style','fill','electrodes','on');
%         end
        
        subj = [subj(end) subj(1:end-1)];


    end

end


%% One other subjects and all other trials
cd('C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\All_data\New bp filter\intermediate')

isr = 8; nsr = 3;

% OPM: changed to not include the data of subject 2
for TRL = 1: length(trll)
    
    load(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    subj = [1, 2, 3, 4];
    trials3 = []; trials2 = [];
   
    for S = 1:4
 
        noftrials = size(Attended{subj(end)}, 1);
        train = (1:noftrials);
        trials4 = subj(1:end-1);
        
        for SO = 1:3;
            trials3 = trials4(SO);
            
        CA = []; CUA = [];
        
        for HH = 1:noftrials % repeat for all trials - leave one out
        
            trials2 = train(1:end-1);
            
            sxx = squeeze(mean(covar.Rxx_all(trials3,trials2,:,:),1)); % plain averaging of cov matrices
            RXX = squeeze(mean(sxx,1));
           
            
            sxy = squeeze(mean(covar.Rxy_at_all(trials3,trials2,:),1));
            RXY_ATT = squeeze(mean(sxy,1))';
            
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE{subj(end)}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{subj(end)}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{subj(end)}(train(end),:)'); % Unattended
            
            train = [train(end) train(1:end-1)];
        end
        
        General_Accuracy(TRL,subj(end),trials4(SO)) = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials
        
%         if wantplots
%             figure
%             subplot(1,3,1)
%             WW = reshape(avg_dec,abs(start - fin)+1,24);
%             WWm = mean(WW,1);
%             topoplot(WWm,'C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\elec_24ch.elp','style','fill','electrodes','on');
%             laggs = mean(WW,2);
%             subplot(1,3,2)
%             plot(laggs)
%             set(gca,'XTick',[1 2.33 4.67 7])
%             set(gca,'XTickLabel',{'1' '100' '200' '300'})
%             xlabel('Time lagg (ms)')
%             [q8 w8] = max(laggs)
%             WWm = mean(WW(w8,:),1);
%             subplot(1,3,3)
%             topoplot(WWm,'C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\elec_24ch.elp','style','fill','electrodes','on');
%         end
        end
        
        subj = [subj(end) subj(1:end-1)];


    end

end

