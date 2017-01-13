%% Parameters
close all 
clear all 

nofsubj = 4;
LP = 2;
HP = 9;
NSR  = 20
trll = [10 20 30 40 50 60];

power = 0.6; % Why? See literature, or ask Neetha...

intSamp = [20 70 120 160 200 300 400 500]; % Why 120? -> See literature + ask Neetha 

%% All the intermediates

for S = 1:2
    for isr = 3:3
        intermediateSampleRate = intSamp(isr)
    ses = 0;
    for LR = 1:2
        for session = 1:4
            ses = ses +1;
            if S == 1 && LR == 1 % Error corretction
                [AudioL,fs] = audioread(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\stimuli\part' num2str(session) '_track2_dry.wav']);
                [AudioR,fs] = audioread(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\stimuli\part' num2str(session) '_track1_dry.wav']);
            else
                [AudioL,fs] = audioread(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\stimuli\part' num2str(session) '_track1_dry.wav']);
                [AudioR,fs] = audioread(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\stimuli\part' num2str(session) '_track2_dry.wav']);
            end
            if LR ==1
                EEG = pop_loadxdf(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\raw_data\S' num2str(S) '_Pilot_' num2str(session) '_cheap_L.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
            elseif LR == 2
                EEG = pop_loadxdf(['C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\raw_data\S' num2str(S) '_Pilot_' num2str(session) '_cheap_R.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
            end
            
            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
            [pks,locs] = findpeaks(double(EEG.data(28,:))); % Locate Triggers
 
            if session == 4 % Auto add trigger after 6min30s because incongruent length of audio stories so missing trigger.
                locs(2) = locs(1) + ((395)*EEGsampleRate);
            else
                locs(2) = locs(1) + ((size(AudioL,1)/fs)*EEGsampleRate);
            end
            EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            
%             % Remove DC component (Done in Neethas script in 'biopil_filter') 
%             % OPM: Dit wel of niet doen???
%             % OPM: Matlab zegt: Do not use filtfilt with differentiator and Hilbert FIR filters, because the operation of these filters depends heavily on their phase response.
%             if remove_DC
%                 [b, a] = butter(1, 0.5/(EEGsampleRate/2), 'high'); % Klopt dit???
%                 for i = 1: 24     
%                     EEGaad(i,:) = filtfilt(b,a,double(EEGaad(i,:)));
%                 end
%             end
            
%             if subband
%                 % OPM in the script of Neetha in aad_envelopes the audio is
%                 % intermediately resampled to an intermediatefs of 8000 Hz. Why? Do we
%                 % need this? Why not immediately to 500 like Rob suggests?
%                 % OPM: dit wordt niet alleen in geval van subbands
%                 % gedaan...
%                 AudioL = resample(AudioL, cut_off*2,fs);
%                 AudioR = resample(AudioR, cut_off*2, fs);
%                 fs = cut_off*2;
% 
%                 % from amtoolbox>joergensen2011.m
%                 g = gammatonefir(gammatonefreqs,fs,[],gammatonebetamul,'real'); % create real, FIR gammatone filters. SEE aad_plot_gammatone
%                 
%                 AudioL = real(ufilterbank(AudioL,g,1)); %audio is filtered along columns. wat is difference with ufilterbankz ? Real or ABS ? See test_gammatone_differences.m
%                 AudioR = real(ufilterbank(AudioR,g,1)); %audio is filtered along columns. wat is difference with ufilterbankz ? Real or ABS ? See test_gammatone_differences.m
% 
%                 % AudioL = reshape(AudioL,size(AudioL,1),[]); % make 2D for
%                 % easier handling % OPM: Niet nodig in dit script??? -> Ask Neetha
%             end
                        % Rob resamples before making the envelope
%             if intermediate
%             else           
%                 % Rob's method: Rob downsamples the audio before making the
%                 % envelope. Neetha does this after making the envelope.
%                 % Does this make a difference? 
%                 AudioL = resample(AudioL,EEGsampleRate,fs); % Downsampled audio per ear
%                 AudioR = resample(AudioR,EEGsampleRate,fs);
%                 fs = EEGsampleRate;  
%             end
            
            % ENVELOPES
            % Powerlaw
         
                envelopeL = abs(AudioL).^power;
                envelopeR = abs(AudioR).^power;
           
            
            % Intermediary Downsample envelope & EEG before applying the more
            % strict bpfilters (to save computation) % OPM: Rob does not do this,
            % but he downsamples the audio to 500... What is the difference? 
                envelopeL = resample(envelopeL,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                envelopeR = resample(envelopeR,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                EEGaad = (resample((double(EEGaad))', intermediateSampleRate, EEGsampleRate))';
                EEGsampleRate = intermediateSampleRate;
                fs = intermediateSampleRate;
            
            New = [];
            % FILTERING
            
            for ch = 1:24 % EEG
                [EEGaad(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad(ch,:)); % Band pass
                New(ch,:) = resample(double(EEGaad(ch,:)),NSR,EEGsampleRate);
            end
            EEGaad = New;
            EEGsampleRate = NSR;
            
            [envelopeL(:)] = Rbp(LP, HP, 6, fs, envelopeL(:,1)); % Band pass
            [envelopeR(:)] = Rbp(LP, HP, 6, fs, envelopeR(:,1)); % Band pass
            envelopeL =  resample(double(envelopeL),NSR,fs);
            envelopeR =  resample(double(envelopeR),NSR,fs);
            fs = NSR;
            
            EEGses{ses} = EEGaad;
            EEGaad = [];
            
            envelopeLses{ses} = envelopeL;
            envelopeRses{ses} = envelopeR;

        end
       
    end

    
    cd('C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\All_data\fixed\intermediate varied')
    save(['S_' num2str(S) '_isr_' num2str(isr) '.mat'], 'EEGses', 'envelopeLses', 'envelopeRses')
    EEGses = {}; envelopeLses = {}; envelopeRses = {};

    end
end

%% Create trials

for S = 1:2
    for isr = 3:3
        
        cd('C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\All_data\fixed\intermediate varied')
        load(['S_' num2str(S) '_isr_' num2str(isr) '.mat'])
        
        for TRL = 5:6
            triallength = trll(TRL);
            
            gentel = 0; gentel6 = 0;
            covar = []; X=[]; ALall = []; ARall = [];
            EE = []; Attended = []; Unattended = [];
            
            for ses = 1:8
                EEGaad = []; envelopeL = []; envelopeR = [];
                
                EEGaad = EEGses{ses};
                envelopeL = envelopeLses{ses};
                envelopeR = envelopeRses{ses};
                
                % CREATE TRIALS
                window = triallength*NSR; % seconds * Sampling-rate
                [trs] = mod(size(EEGaad,2),window);
                trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
                trs = [];
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
                    gentel = gentel+1;
                    AL(i,:) = envelopeL(s1+1:s1+window);
                    AR(i,:) = envelopeR(s1+1:s1+window);
                    if ses < 5 
                        Attended(gentel,:) = AL(i,:);
                        Unattended(gentel,:) = AR(i,:);
                    elseif ses > 4
                        Attended(gentel,:) = AR(i,:);
                        Unattended(gentel,:) = AL(i,:);
                    end
                    s1 =s1+window;           
                
                end
                
                % Lags
                test = {}; Rxx = {}; Rxy_left = {}; Rxy_right = {}; x = [];
                for trials = 1:trn
                    gentel6 = gentel6 + 1;
                    x=squeeze(E(trials,:,:))';
                    start = floor(100/1e3*NSR); %convert milliseconds to samples 100
                    fin = ceil(400/1e3*NSR); %convert milliseconds to samples 400
                    nofsamples = size(x,1);
                    [start,fin] = deal(-fin,-start);
                    [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                
                    EE = cat(3,EE,X);
                
                    % Covariance Matrices
                    Rxx{trials} = (X'*X)/nofsamples;
                
                    Rxy_left{trials} = (X'*squeeze(AL(trials,:))')/nofsamples;
                    Rxy_right{trials} = (X'*squeeze(AR(trials,:))')/nofsamples;
                
                    % Ryy_left{trials} = ((AL(1,:))*squeeze(AL(trials,:))')/nofsamples;
                    % Ryy_right{trials} = ((AR(1,:))*squeeze(AR(trials,:))')/nofsamples;
                
                
                    covar.Rxx_all(gentel6,:,:) = Rxx{trials};
                    if ses < 5
                        covar.Rxy_at_all(gentel6,:) = Rxy_left{trials};
                    elseif ses > 4
                        covar.Rxy_at_all(gentel6,:) = Rxy_right{trials};
                    end
                end
            
            
            end
            
        noftrials = size(Attended,1);
        train = (1:noftrials);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            trials2 = train(1:end-1);
        
            RXX = squeeze(mean(covar.Rxx_all(trials2,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(trials2,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended(train(end),:)'); % Unattended
            train = [train(end) train(1:end-1)];
        end
        
        Accuracy(isr,TRL) = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials
        CA = []; CUA = [];
            
        end
        
    end
    
    cd('C:\Users\Annelies\OneDrive\Studie\Master\thesis\script\All_data\fixed\intermediate varied')
    save(['Result_S' num2str(S) '.mat'], 'Accuracy')
    Accuracy = [];
end

