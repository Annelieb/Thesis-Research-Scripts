%% Evaluate online; Subject accuracy

% Initial parameters

clear all

Subjects = 7; % Subject that needs to be evaluated
startleft = [1 0 1 0 1 0 0];
sessions = 4; 

TRL = [1 3 6]; 

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;
%triallength = 10

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; 

%%
for S = Subjects    
    for TF = 1:2
        EEGaad = {}; Att = {}; Unatt ={};
        for session = 1:4
            
            audiodata = []; X=[]; EE = [];
            EEGaad_data2 ={}; AL = {}; AR = {};
            
            % Load & process the audio
            [audio_data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [audio_data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            audiodata = cat(2,audio_data2(:,1),audio_data3(:,1));
            
            if session < 3
                AudioL = audiodata(:,1); AudioR = audiodata(:,2);
            elseif session > 2
                AudioR = audiodata(:,1); AudioL = audiodata(:,2);
            end
            
            envelopeL = abs(AudioL).^power;
            envelopeR = abs(AudioR).^power;
            
            envelopeL = resample(envelopeL,500,fs);
            envelopeR = resample(envelopeR,500,fs);
            fs = 500;
            
            % Load the cut EEG
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if startleft(S) == 1
                ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
            else
                ear = [0 0 1 1; 1 1 0 0];
            end
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'data2', 'WIN')
            elseif TF == 2
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'data2', 'WIN')
            end
            
            EEGsampleRate = 500;
            for i = 1: size(data2,2)
                k = find(data2{i}(25,:));
                if isempty(k)
                    EEGaad_data2{i} = data2{i}(1:24,:);
                else
                    EEGaad_data2{i} = [];
                end
            end
            EEGaad_data2(cellfun('isempty',EEGaad_data2)) = []; %Remove empty cells
            
            if size(EEGaad_data2,2) > 39 % Get rid of extra data when longer than 6.30 min
                EEGaad_data2(40:end) = [];
            end
            
            % Cut the audio
            WIN(cellfun('isempty', WIN)) = [];
            for trl = 1:39
                strt = WIN{trl}+1; % for some reason +2 is better than just +1. Don't know why?
                stp = WIN{trl+1};
                if strt < stp
                    AL{trl} = envelopeL(strt:stp);
                    AR{trl} = envelopeR(strt:stp);
                    AL{trl} = double(AL{trl});
                    AR{trl} = double(AR{trl});                
                end
            end
            
            % Put all sessions together
            EEGaad = [EEGaad, EEGaad_data2];
            if strcmp(EAR, 'left')
                Att = [Att, AL];
                Unatt = [Unatt, AR];
            elseif strcmp(EAR, 'right')
                Att = [Att, AR];
                Unatt = [Unatt, AL];
            end
        end
        
        check = 1
        
        for trll = 1:length(TRL)
            clear E EE Rxx Rxy_at New
            
            noftrials = floor(size(Att,2)/TRL(trll));
            Att_recut = cell(1, noftrials);
            Unatt_recut = cell(1, noftrials);
            EEGaad_recut = cell(1, noftrials);
            
            for i = 1: noftrials
                for k = 1:TRL(trll)
                    Att_recut{i} = [Att_recut{i}; Att{i+k-1}];
                    Unatt_recut{i} = [Unatt_recut{i}; Unatt{i+k-1}];
                    EEGaad_recut{i} = cat(2, EEGaad_recut{i}, EEGaad{i+k-1});
                end
            end
            
            check = 3
            
            %FILTERING and RESAMPLING
            %EEG
            for trl = 1:noftrials
                for ch = 1:24                  
                    E{trl}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, EEGaad_recut{trl}(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                    New(ch,:) = resample(double(E{trl}(ch,:)),NSR,EEGsampleRate);
                end
                E{trl} = New;                
                New = [];
            end 
            
            %Audio
            for trl = 1:noftrials
                Att_recut{trl} = Rbp([], HP, 6, fs, Att_recut{trl});
                Unatt_recut{trl} = Rbp([], HP, 6, fs, Unatt_recut{trl});
                
                Att_recut{trl} = resample(double(Att_recut{trl}), NSR, fs); 
                Unatt_recut{trl} = resample(double(Unatt_recut{trl}), NSR, fs);
            end
            
            check = 4
            
            for trl = 1:noftrials               
                % Lags                          
                x = E{trl}';
                start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                [start,fin] = deal(-fin,-start);
                nofsamples = size(x,1);
                [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                
                EE{trl} = X;
                
                % Covariance Matrices
                Rxx(trl,:,:) = (X'*X)/nofsamples;
                                
                Rxy_at(trl,:) = (X'*Att_recut{trl}')/nofsamples;
            end
            
            check = 5
            
            train = 1:noftrials;
            
            for HH = 1:noftrials % repeat for all trials - leave one out
                trials2 = train(1:end-1);

                RXX = squeeze(mean(Rxx(trials2,:,:),1)); % plain averaging of cov matrices
                RXY_ATT = squeeze(mean(Rxy_at(trials2,:),1));

                avg_dec = RXX \ RXY_ATT'; % LS solution.
                recenv = EE{train(end)} * avg_dec;

                CA{trll}(S,TF,train(end)) = corr2(recenv,Att_recut{train(end)}'); % Attended
                CUA{trll}(S,TF,train(end)) = corr2(recenv,Unatt_recut{train(end)}'); % Unattended
                train = [train(end) train(1:end-1)];
            end
            
            check = 6
            
            if TF == 1
                difference_T{trll}(S,:) = CA{trll}(S,TF,:) - CUA{trll}(S,TF,:);
                Accuracy_T(S,trll) = sum(difference_T{trll}(S,:)>0)/noftrials*100
            elseif TF == 2
                difference_F{trll}(S,:) = CA{trll}(S,TF,:) - CUA{trll}(S,TF,:);
                Accuracy_F(S,trll) = sum(difference_F{trll}(S,:)>0)/noftrials*100
            end
        end
    end
end
