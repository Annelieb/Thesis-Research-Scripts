%% Cut first

clear all

nofsubj = 4;
lag = [100 400];                                                    
LP = 2;
HP = 9;
%LP = [0, 1, 2, 3, 4];
%HP = [6,7, 8, 9, 10];
NSR = 20;
%NewSamplingRate = [10 15 20 25 30 35 40 45 50 80]; % 200 500]; % New sampling rate; 200 and 500 take to much computing time
trll = [10 20 30 40 50 60]; % in seconds

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?
intermediateSampleRate = 500
 subband = false; 
% spacing = 1.5; % OPM: How is this defined? Ask Neetha...
% gammatonefreqs = erbspacebw(150,4000,spacing); % gammatone filter centerfrequencies % OPM: What is all off this? How are the inputs defined? Ask Neetha...
% gammatonebetamul = spacing*ones(size(gammatonefreqs)); % multiplier for gammatone filter bandwidths

intermediate = true; % Neetha does this, check what the effect is
%intSamp = [20 70 120 160 200 300 400 500]; % Why 120?
%cut_off = 4000; %cut-off frequency of headphones. OPM: Find what this is in our experiment.

remove_DC = false; % Neetha does this, check what the effect is

%% Get the covariances

for TRL = 1:3 % Loop for trial length
    
    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
    load(['TRL' num2str(TRL) '_ISR8_NSR3.mat'], 'Attended', 'Unattended')
    Attended_nc = Attended; Unattended_nc = Unattended;
    clear Attended Unattended
    triallength = trll(TRL)
    covar = []; X=[]; ALall = []; ARall = [];
    EE = cell(1, nofsubj); Attended = cell(1, nofsubj); Unattended = cell(1, nofsubj);
    Attended_cut = cell(1, nofsubj); Unattended_cut = cell(1, nofsubj);
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
            
            AudioL = data(:,1); AudioR = data(:,2);
            data = [];
            
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
            
               envelopeL = resample(envelopeL,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
               envelopeR = resample(envelopeR,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                fs = intermediateSampleRate;
                
            % CREATE TRIALS
            window = triallength*intermediateSampleRate; % seconds * Sampling-rate
            [trs] = mod(size(EEGaad,2),window);
            trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
            sessionlength(S, ses) = trn; 
           
            % EEG
            s1 = 0; E = [];
            s1 = 0; AL = []; AR = [];
            for i = 1:trn
                for ch = 1:24
                    E(i,ch,:) = EEGaad(ch,s1+1:s1+window);
                    [E(i,ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, squeeze(E(i,ch,:))); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                    New(i,ch,:) = resample(double(E(i,ch,:)),NSR,EEGsampleRate);
                end
                ok = 1
                      
            % AUDIO
                gentel5 = gentel5+1;
                AL(i,:) = envelopeL(s1+1:s1+window);
                AR(i,:) = envelopeR(s1+1:s1+window);
                AL(i,:) = Rbp([], HP, 6, fs, AL(i,:)); % Band pass
                AR(i,:) = Rbp([], HP, 6, fs, AR(i,:)); % Band pass
                NewAL(i,:) =  resample(double(AL(i,:)),NSR,fs);
                NewAR(i,:) =  resample(double(AR(i,:)),NSR,fs);
            
                if left
                    Attended_cut{S}(gentel5,:) = NewAL(i,:);
                    Unattended_cut{S}(gentel5,:) = NewAR(i,:);
                else
                    Attended_cut{S}(gentel5,:) = NewAR(i,:);
                    Unattended_cut{S}(gentel5,:) = NewAL(i,:);
                end
                s1 =s1+window;         
                
            end
            
            E = New; AL = NewAL; AR = NewAR;          
            New = []; NewAL = []; NewAR = [];

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
                
                covar.Rxy_at_all_cut(S, gentel6,:) = (X'*squeeze(Attended_cut{S}(gentel6,:))')/nofsamples;
                covar.Rxy_at_all_nc(S, gentel6,:) = (X'*squeeze(Attended_nc{S}(gentel6,:))')/nofsamples;              
                
                covar.Rxx_all(S, gentel6,:,:) = Rxx{trials};

            end
           
        end

    end
 
end
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\cut first')

%save(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'], 'EE', 'covar', 'Attended', 'Unattended', 'sessionlength') % , '-v7.3') Vanaf NSR > 50 nodig
save(['TRL' num2str(TRL) '.mat'], 'EE', 'covar', 'Attended_cut', 'Unattended_cut', 'Attended_nc', 'Unattended_nc', 'sessionlength') % , '-v7.3') Vanaf NSR > 50 nodig

clear EE Attended Unattended

end

% Subject accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\cut first')
Subject_Accuracy = [];

for TRL = 1:3
    
    load(['TRL' num2str(TRL) '.mat'])

    for S = 1:nofsubj
        
        noftrials = size(Attended_cut{S},1);
        train = (1:noftrials);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            trials2 = train(1:end-1);
        
            RXX = squeeze(mean(covar.Rxx_all(S,trials2,:,:),2)); % plain averaging of cov matrices
            RXY_ATT_cut = squeeze(mean(covar.Rxy_at_all_cut(S,trials2,:),2));
            RXY_ATT_nc = squeeze(mean(covar.Rxy_at_all_nc(S,trials2,:),2));

            avg_dec_cut = RXX \ RXY_ATT_cut; % LS solution.
            avg_dec_nc = RXX \ RXY_ATT_nc; % LS solution.
            
            recenv_cut = EE{S}(:,:,train(end)) * avg_dec_cut;
            recenv_nc = EE{S}(:,:,train(end)) * avg_dec_nc;
        
            CA_cut(HH) = corr2(recenv_cut,Attended_cut{S}(train(end),:)'); % Attended
            CUA_cut(HH) = corr2(recenv_cut,Unattended_cut{S}(train(end),:)'); % Unattended
            
            CA_nc(HH) = corr2(recenv_nc,Attended_nc{S}(train(end),:)'); % Attended
            CUA_nc(HH) = corr2(recenv_nc,Unattended_nc{S}(train(end),:)'); % Unattended

            train = [train(end) train(1:end-1)];
        end
        
        Subject_Accuracy_cut(TRL,S) = (sum(CA_cut-CUA_cut> 0)/(noftrials))*100 % Final accuracy on all trials
        Subject_Accuracy_nc(TRL,S) = (sum(CA_nc-CUA_nc> 0)/(noftrials))*100 % Final accuracy on all trials

        CA_cut = []; CUA_cut = [];CA_nc = []; CUA_nc = [];

    end

end

cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\cut first')
        save(['Result.mat'], 'Subject_Accuracy_cut', 'Subject_Accuracy_nc')
