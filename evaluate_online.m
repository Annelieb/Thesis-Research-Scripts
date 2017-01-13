%% Initial parameters

clear all

S = 1; % Subjects that need to be evaluated
sessions = 1:4; 

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;
TRL = [10 30 60]


envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?

%% Get covariances    
for triallength = TRL
    for TF = 1:2 % Left - Right loop
        if startleft
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end

        for session = 1:4 % 4 right ear session - Train on each one separate
            data = [];
            covar = []; X=[];
            EE = []; Attended = []; Unattended = [];
            sessionlength = [];

            
            [data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            data = cat(2,data2(:,1),data3(:,1));
            
            if session < 3
                AudioL = audiodata(:,1); AudioR = audiodata(:,2);
            elseif session > 2
                AudioR = audiodata(:,1); AudioL = audiodata(:,2);
            end
  
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
           
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            
            if TF == 1
                EEG = pop_biosig(['Training_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
            elseif TF == 2
                EEG = pop_biosig(['Feedback_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
            end
            
            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
            locs(1) = EEG.event(1).latency;
 
            if session == 4 % Auto add trigger after 6min30s because incongruent length of audio stories so missing trigger.
                locs(2) = locs(1) + ((395)*EEGsampleRate);
            else
                locs(2) = locs(1) + ((size(data,1)/fs)*EEGsampleRate);
                % locs(2) = EEG.event(2).latency; % OPM: Normaal zouden
                % deze twee hetzelfde resultaat moeten geven. Toch
                % verschilt het lichtjes... 
            end
             
            
            EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            
            
            
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
            
                envelopeL = resample(envelopeL,500,fs); % OPM: Is dit nodig in ons script?
                envelopeR = resample(envelopeR,500,fs); % OPM: Is dit nodig in ons script?
                fs = 500;              
            
            New = [];
            % FILTERING
                        
            for ch = 1:24 % EEG
                [EEGaad(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                New(ch,:) = resample(double(EEGaad(ch,:)),NSR,EEGsampleRate);
            end
            EEGaad = New;
            EEGsampleRate = NSR;
                                   
            [envelopeL(:)] = Rbp([], HP, 6, fs, envelopeL(:,1)); % Band pass
            [envelopeR(:)] = Rbp([], HP, 6, fs, envelopeR(:,1)); % Band pass
            
            envelopeL =  resample(double(envelopeL),NSR,fs);
            envelopeR =  resample(double(envelopeR),NSR,fs);
            
            fs = NSR;
                        
            % CREATE TRIALS
            window = triallength*NSR; % seconds * Sampling-rate
            [trs] = mod(size(EEGaad,2),window);
            trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
            sessionlength = trn; 
           

            % EEG
            s1 = 0; E = [];
            for i = 1:trn
                for ch = 1:24
                    E(i,ch,:) = EEGaad(ch,s1+1:s1+window);
                end
                s1 = s1+window;
            end
            
            % AUDIO
            s1 = 0; AL = []; AR = []; % for some reason index is out of bounds in this way... 
            for i = 1:trn

                AL(i,:) = envelopeL(s1+1:s1+window);
                AR(i,:) = envelopeR(s1+1:s1+window);
               s1 =s1+window;          
            end
            
            if ear(TF, session) == 1
                left = true; 
            else
                left = false;
            end
                if left
                    Attended = AL;
                    Unattended = AR;
                else
                    Attended = AR;
                    Unattended = AL;
                end
                

            % Lags
            test = {}; Rxx = {}; Rxy_left = {}; Rxy_right = {}; x = [];
            for trials = 1:trn
                

                x=squeeze(E(trials,:,:))';
                start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                nofsamples = size(x,1);
                [start,fin] = deal(-fin,-start);
                [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                
                EE = cat(3,EE,X);
                
                % Covariance Matrices
                Rxx{trials} = (X'*X)/nofsamples;
                
                Rxy_left{trials} = (X'*squeeze(AL(trials,:))')/nofsamples;
                Rxy_right{trials} = (X'*squeeze(AR(trials,:))')/nofsamples;
                
                covar.Rxx_ses(trials,:,:) = Rxx{trials};
                if left
                    covar.Rxy_at_ses(trials,:) = Rxy_left{trials};
                else
                    covar.Rxy_at_ses(trials,:) = Rxy_right{trials};
                end
            end
            
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Off line analysis'])

            if TF == 1
                save(['OFFCUT103060_Training_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
            elseif TF == 2
                save(['OFFCUT103060_Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
            end
            clear EE Attended Unattended
           
        end

    end
end



%% Session accuracy
Session_Accuracy = [];

% cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\General decoders')
% load(['TRL1_Decoder_S3.mat'])
% avg_dec3 = avg_dec;
% avg_dec =[];
% 
% cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\General decoders')
% load(['TRL1_Decoder_S4.mat'])
% avg_dec4 = avg_dec;
% avg_dec = [];

cd('C:\Users\Annelies\Documents\Studie\Thesis\script')
load('TEST_S4.mat')
avg_dec4 = avg_dec;
avg_dec = [];



    for LR = 1:2 % Left - Right loop
        if LR == 1
            left = true;
        else
            left = false;
        end
    for session = sessions
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\Online Rob\pilot_2\Analysis')
        if left
            load(['Left_Session_' num2str(session) '.mat'], 'EE', 'covar', 'Attended', 'Unattended', 'sessionlength') 
            else
            load(['Right_Session_' num2str(session) '.mat'], 'EE', 'covar', 'Attended', 'Unattended', 'sessionlength') 
            end
        
        noftrials = sessionlength;
        train = 1:noftrials;
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            
            trials2 = train(1:end-1);
            
            % Accuracy, decoder trained on self
            RXX = squeeze(mean(covar.Rxx_ses(trials2,:,:),1)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_ses(trials2,:),1));
        
            avg_dec_self = RXX \ RXY_ATT'; % LS solution.
            recenv_self = EE(:,:,train(end)) * avg_dec_self;
            CA_self(HH) = corr2(recenv_self,Attended(train(end),:)'); % Attended
            CUA_self(HH) = corr2(recenv_self,Unattended(train(end),:)'); % Unattended
            
%             recenv_dec3 = EE(:,:,train(end)) * avg_dec3;
%             CA_dec3(HH) = corr2(recenv_dec3,Attended(train(end),:)'); % Attended
%             CUA_dec3(HH) = corr2(recenv_dec3,Unattended(train(end),:)'); % Unattended
            
            recenv_dec4 = EE(:,:,train(end)) * avg_dec4;
            CA_dec4(HH) = corr2(recenv_dec4,Attended(train(end),:)'); % Attended
            CUA_dec4(HH) = corr2(recenv_dec4,Unattended(train(end),:)'); % Unattended
            
            train = [train(end) train(1:end-1)];
        end
        
        Session_Accuracy_self(session, LR) = (sum(CA_self-CUA_self> 0)/(noftrials))*100 % Final accuracy on all trials
        %Session_Accuracy_dec3(session, LR) = (sum(CA_dec3-CUA_dec3> 0)/(noftrials))*100 % Final accuracy on all trials
        Session_Accuracy_dec4(session, LR) = (sum(CA_dec4-CUA_dec4> 0)/(noftrials))*100 % Final accuracy on all trials

        CA_self = []; CUA_self = [];
        %CA_dec3 = []; CUA_dec3 = [];
        CA_dec4 = []; CUA_dec4 = [];
    end
    end

        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\Online Rob\pilot_2\Analysis')
        save(['Result_session_accuracies.mat'], 'Session_Accuracy_self', 'Session_Accuracy_dec3', 'Session_Accuracy_dec4')
%% Subject accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\Online Rob\Analysis')
load(['Covars_etc.mat'])
Subject_Accuracy = [];

RXX = squeeze(mean(covar.Rxx_all(1,:,:,:),2)); % plain averaging of cov matrices
RXY_ATT = squeeze(mean(covar.Rxy_at_all(1,:,:),2));
        
avg_dec = RXX \ RXY_ATT; % LS solution.

    for S = 1:nofsubj
        
        noftrials = size(Attended{S},1);
        train = (1:noftrials);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
%             trials2 = train(1:end-1);
%         
%             RXX = squeeze(mean(covar.Rxx_all(S,trials2,:,:),2)); % plain averaging of cov matrices
%             RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,trials2,:),2));
        
%             avg_dec = RXX \ RXY_ATT; % LS solution.
%             cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\General decoders')
%             load(['TRL1_Decoder_S3.mat'])
            recenv = EE{S}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{S}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{S}(train(end),:)'); % Unattended
            train = [train(end) train(1:end-1)];
        end
        
        Subject_Accuracy(S) = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials
        CA = []; CUA = [];


    end

        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\Online Rob\Analysis')
        save(['Result_dec_van_eerste_online.mat'], 'Subject_Accuracy')
