%% Evaluate online; correct cutting

% Initial parameters

clear all

startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;
%triallength = 10

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; 

%% Get covariances    
for S = 1:12
    for TF = 1:2 % Training Feedback loop

        if startleft(S)
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end
        
        for session = 1:4 % 4 right ear session - Train on each one separate
            audiodata = []; covar = []; X=[]; EE = []; Attended = []; Unattended = [];
            EEGaad_data2 ={}; EEGaad = [];
            
            [audio_data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [audio_data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            audiodata = cat(2,audio_data2(:,1),audio_data3(:,1));
            
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
                load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'data2', 'WIN', 'Attended', 'Unattended')
            elseif TF == 2
                EEG = pop_biosig(['Feedback_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'data2', 'WIN', 'Attended', 'Unattended')
            end
            
            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
            locs(1) = EEG.event(1).latency;
            if session == 4
                locs(2) = locs(1) + ((395)*EEGsampleRate);
            else
                locs(2) = EEG.event(2).latency;
            end
            
            EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            
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
                Attended(40:end) = []; 
                Unattended(40:end) = [];
            end
            
            audiodata = [];
            
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
            
                envelopeL = resample(envelopeL,500,fs);
                envelopeR = resample(envelopeR,500,fs);
                fs = 500;              
            
            New = [];
            
            % CUT BEFORE FILTERING
            E = {};Rxx = {}; x = []; Rxy_at = {}; AL = {}; AR = {};
            WIN(cellfun('isempty',WIN)) = []; %Remove empty cells
            size_WIN = size(WIN,2)-1

            for trl = 1:size_WIN
                strt = WIN{trl}+1; % for some reason +2 is better than just +1. Don't know why?
                stp = WIN{trl+1};
                if stp <= size(EEGaad,2)
                    AL{trl} = envelopeL(strt:stp);
                    AR{trl} = envelopeR(strt:stp);
                    AL{trl} = double(AL{trl});
                    AR{trl} = double(AR{trl});                
                for ch = 1:24
                    E{trl}(ch,:) = EEGaad(ch, strt:stp);
                end
                end
            end
            
             
            E(cellfun('isempty',E)) = []; %Remove empty cells
            AL(cellfun('isempty',AL)) = []; %Remove empty cells
            AR(cellfun('isempty',AR)) = []; %Remove empty cells

            noftrials = size(E,2);
            noftrials_data2 = size(EEGaad_data2,2);
            
            if noftrials < noftrials_data2
                EEGaad_data2(noftrials+1:end) = [];
            elseif noftrials > noftrials_data2
                error('noftrials > noftrials_data2')
                break
            end
            

            
            if noftrials ~= noftrials_data2
                error('Still a mistake')
            elseif noftrials == noftrials_data2
                for i = 1: noftrials
                    sizes_E(i) = size(E{i},2);
                    sizes_data2(i) = size(EEGaad_data2{i},2);
                    if isequal(sizes_E, sizes_data2)
                    else
                        error('sizes_E and sizes_data2 are not equal')
                    end
                end
            end
           
             E_data2 = {};
            %FILTERING and RESAMPLING
            for trl = 1:noftrials_data2
                for ch = 1:24                  
                    E{trl}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, E{trl}(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                    E_data2{trl}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, EEGaad_data2{trl}(ch,:));
                    New(ch,:) = resample(double(E{trl}(ch,:)),NSR,EEGsampleRate);
                    New_data2(ch,:) = resample(double(E_data2{trl}(ch,:)),NSR,EEGsampleRate);
                end
                
                E{trl} = New;
                E_data2{trl} = New_data2;
                
                New = []; 
                New_data2 = [];                
                                   
                AL{trl} = Rbp([], HP, 6, fs, AL{trl}); % Band pass
                AR{trl} = Rbp([], HP, 6, fs, AR{trl}); % Band pass
            
                AL{trl} =  resample(double(AL{trl}),NSR,fs);
                AR{trl} =  resample(double(AR{trl}),NSR,fs);           
                
            end

            
            for trl = 1:noftrials_data2               
                % Lags                          
                x = E{trl}';
                x_data2 = E_data2{trl}';
                start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                [start,fin] = deal(-fin,-start);
                nofsamples = size(x_data2,1);
                [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                [X_data2] = aad_LagGenerator(x_data2,start:fin); % Time lags of 2-500ms in samples
                
                EE{trl} = X;
                EE_data2{trl} = X_data2;
                
                % Covariance Matrices
                Rxx{trl} = (X'*X)/nofsamples;
                Rxx_data2{trl} = (X_data2'*X_data2)/nofsamples;
                
                if ear(TF, session) == 1
                Rxy_at_made{trl} = (X'*squeeze(AL{trl})')/nofsamples;
                Rxy_at_data2_made{trl} = (X_data2'*squeeze(AL{trl})')/nofsamples;
                
                else
                Rxy_at_made{trl} = (X'*squeeze(AR{trl})')/nofsamples;
                Rxy_at_data2_made{trl} = (X_data2'*squeeze(AR{trl})')/nofsamples;
                end
                
                Rxy_at{trl} = (X'*squeeze(Attended{trl}))/nofsamples;
                Rxy_at_data2{trl} = (X_data2'*squeeze(Attended{trl}))/nofsamples;

            end
            
            if ear(TF, session) == 1
                Attended_made = AL;
                Unattended_made = AR;
                
            else
                Attended_made = AR;
                Unattended_made = AL;
                
            end
            
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])

            if exist('Off line analysis','dir') == 7
            else
                mkdir('Off line analysis')
            end
            
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Off line analysis'])
            
            if TF == 1
                save(['Training_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
            elseif TF == 2
                save(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
            end            
            
            clear EE Attended Unattended EE_data2 Attended_made Unattended_made
           
            end
    end
end
    
    
    
    
    
    %% Session accuracy

Session_Accuracy = [];

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
load(['Best_decoder_Subject' num2str(S) '.mat'])
avg_best = avg_dec;
avg_dec = [];

        difference_self_data2 = [];
        difference_self = [];
%         difference_dec4_data2 = [];
%         difference_dec4 = [];

for TF = 1:2 % Training Feedback loop

        if startleft
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end

     for session = 1:4

            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Off line analysis'])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Training_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
            
                noftrials = size(Attended_made,2);
                train = 1:noftrials;
    
                for HH = 1:noftrials % repeat for all trials - leave one out
            
                    trials2 = train(1:end-1);

                    % Accuracy, decoder trained on self
                    RXX = zeros(size(Rxx{1})); RXY_ATT = zeros(size(Rxy_at{1}));RXY_ATT_made = zeros(size(Rxy_at{1}));
                    RXX_data2 = zeros(size(Rxx_data2{1})); RXY_ATT_data2 = zeros(size(Rxy_at_data2{1}));RXY_ATT_data2_made = zeros(size(Rxy_at_data2{1}));
                    
                    for i = trials2
                        RXX = RXX + Rxx{i}; % plain averaging of cov matrices
                        RXY_ATT = RXY_ATT + Rxy_at{i};
                        RXY_ATT_made = RXY_ATT_made + Rxy_at_made{i};

                        RXX_data2 = RXX_data2 + Rxx_data2{i}; % plain averaging of cov matrices
                        RXY_ATT_data2 = RXY_ATT_data2 + Rxy_at_data2{i};
                        RXY_ATT_data2_made = RXY_ATT_data2_made + Rxy_at_data2_made{i};
                    end
                    
                    % gdf
                    RXX = RXX/(noftrials-1);
                    RXY_ATT = RXY_ATT/(noftrials-1);
                    RXY_ATT_made = RXY_ATT_made/(noftrials-1);
                    
                    % data2
                    RXX_data2 = RXX_data2/(noftrials-1);
                    RXY_ATT_data2 = RXY_ATT_data2/(noftrials-1);
                    RXY_ATT_data2_made = RXY_ATT_data2_made/(noftrials-1);
                    
                    % gdf
                    avg_dec_self = RXX \ RXY_ATT; % LS solution.
                    recenv_self = EE{train(end)} * avg_dec_self;
                    
                    CA_self(HH) = corr2(recenv_self,Attended{train(end)}); % Attended
                    CUA_self(HH) = corr2(recenv_self,Unattended{train(end)}); % Unattended
%             
                    avg_dec_made = RXX \ RXY_ATT_made; % LS solution.
                    recenv_self_made = EE{train(end)} * avg_dec_made;
                    
                    CA_self_made(HH) = corr2(recenv_self_made,Attended_made{train(end)}'); % Attended
                    CUA_self_made(HH) = corr2(recenv_self_made,Unattended_made{train(end)}'); % Unattended
                    
                    % data2
                    avg_dec_self_data2 = RXX_data2 \ RXY_ATT_data2; % LS solution.
                    recenv_self_data2 = EE_data2{train(end)} * avg_dec_self_data2;
                    
                    CA_self_data2(HH) = corr2(recenv_self_data2,Attended{train(end)}); % Attended
                    CUA_self_data2(HH) = corr2(recenv_self_data2,Unattended{train(end)}); % Unattended
                                 
                    avg_dec_self_data2_made = RXX_data2 \ RXY_ATT_data2_made; % LS solution.
                    recenv_self_data2_made = EE_data2{train(end)} * avg_dec_self_data2_made;
                    
                    CA_self_data2_made(HH) = corr2(recenv_self_data2_made,Attended_made{train(end)}'); % Attended
                    CUA_self_data2_made(HH) = corr2(recenv_self_data2_made,Unattended_made{train(end)}'); % Unattended
            
                    train = [train(end) train(1:end-1)];
                end
                
                Session_Accuracy_Training(session, :) = [(sum(CA_self-CUA_self> 0)/(noftrials))*100, (sum(CA_self_made-CUA_self_made> 0)/(noftrials))*100, (sum(CA_self_data2-CUA_self_data2> 0)/(noftrials))*100, (sum(CA_self_data2_made-CUA_self_data2_made> 0)/(noftrials))*100]

            elseif TF == 2
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2','Rxy_at_made', 'Rxy_at_data2_made', 'Attended', 'Unattended', 'Attended_made', 'Unattended_made')
                
                noftrials = size(EE,2);
                for i = 1: noftrials
                    recenv = EE{i} * avg_best;
                    CA_best(i) = corr2(recenv,Attended{i}); % Attended
                    CUA_best(i) = corr2(recenv,Unattended{i}); % Unattended

                    CA_best_made(i) = corr2(recenv,Attended_made{i}'); % Attended
                    CUA_best_made(i) = corr2(recenv,Unattended_made{i}'); % Unattended

                    recenv_data2 = EE_data2{i} * avg_best;
                    CA_best_data2(i) = corr2(recenv_data2,Attended{i}); % Attended
                    CUA_best_data2(i) = corr2(recenv_data2,Unattended{i}); % Unattended            

                    CA_best_data2_made(i) = corr2(recenv_data2,Attended_made{i}'); % Attended
                    CUA_best_data2_made(i) = corr2(recenv_data2,Unattended_made{i}'); % Unattended            
                end
                
                Session_Accuracy_Feedback(session, :) = [(sum(CA_best-CUA_best> 0)/(noftrials))*100, (sum(CA_best_made-CUA_best_made> 0)/(noftrials))*100, (sum(CA_best_data2-CUA_best_data2> 0)/(noftrials))*100, (sum(CA_best_data2_made-CUA_best_data2_made> 0)/(noftrials))*100]
            end
     end
end

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Off line analysis'])
save(['Accuracy'], 'Session_Accuracy_Training', 'Session_Accuracy_Feedback')

%% Online session accuracies samenvoegen voor duidelijkheid
S = 
startleft = false

for TF = 1:2 % Training Feedback loop

        if startleft
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end

     for session = 1:4

            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'Session_Accuracy')
                Session_Accuracy_online_Training(session,1) = Session_Accuracy
            else
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'Session_Accuracy') 
                Session_Accuracy_online_Feedback(session,1) = Session_Accuracy
            end
     end
end

%% Online session accuracies samenvoegen voor duidelijkheid
S = 7
startleft = false

for TF = 1:2 % Training Feedback loop

        if startleft
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end

     for session = 1:4

            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'Session_Accuracy')
                Session_Accuracy_online_Training(session,1) = Session_Accuracy
            else
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'Session_Accuracy') 
                Session_Accuracy_online_Feedback(session,1) = Session_Accuracy
            end
     end
end

%% Online session accuracies based on 'differences'

for TF = 1:2 % Training Feedback loop

        if startleft
            ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
        else
            ear = [0 0 1 1; 1 1 0 0];
        end

     for session = 1:4

            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                Session_Accuracy_online_Training(session,1) = sum(difference(1:39) > 0)/39*100           
            else
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference') 
                Session_Accuracy_online_Feedback(session,1) = sum(difference(1:39) > 0)/39*100
            end
     end
end