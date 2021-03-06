%% Evaluate online; correct cutting

% Initial parameters

clear all

subjects = 1; % Subjects that need to be evaluated
sessions = 2; 

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;
%triallength = 10

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?

%% Get covariances    
    for LR = 1 % Left - Right loop
        if LR == 1
            left = true;
        else
            left = false;
        end

        for session = 2 % 4 right ear session - Train on each one separate
            if LR == 1 && session == 1
            else
            audiodata = []; covar = []; X=[]; EE = []; Attended = []; Unattended = [];
            EEGaad_data2 ={}; EEGaad = [];
            
            [audio_data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [audio_data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            audiodata = cat(2,audio_data2(:,1),audio_data3(:,1));
                        
            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data')
                    if left
                        %EEG = pop_biosig(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\pilot_2_sessie' num2str(session) '_L.gdf']);
                        load(['pilot_2_sessie2_left.mat'])
                    else
                        %EEG = pop_biosig(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\pilot_2_sessie' num2str(session) '_R.gdf']);
                        %load(['pilot_2_sessie' num2str(session) '_righ.mat'])
                    end

            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
%             locs(1) = EEG.event(1).latency;
%             if session == 4
%                 locs(2) = locs(1) + ((395)*EEGsampleRate);
%             else
%                 locs(2) = EEG.event(2).latency;
%             end
            
%             EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            
            for i = 1: size(data2,2)
                k = find(data2{i}(25,:));
                if isempty(k)
                    EEGaad_data2{i} = data2{i}(1:24,:);
                else
                    EEGaad_data2{i} = [];
                end
            end
            
            EEGaad_data2(cellfun('isempty',EEGaad_data2)) = []; %Remove empty cells
 
            AudioL = audiodata(:,1); AudioR = audiodata(:,2);
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
%                 %if stp <= size(EEGaad,2)
                AL{trl} = envelopeL(strt:stp);
                AR{trl} = envelopeR(strt:stp);
                AL{trl} = double(AL{trl});
                AR{trl} = double(AR{trl});
                 
%                 for ch = 1:24
%                     E{trl}(ch,:) = EEGaad(ch, strt:stp);
%                 end
                %end
            end
            
             
            %E(cellfun('isempty',E)) = []; %Remove empty cells
            AL(cellfun('isempty',AL)) = []; %Remove empty cells
            AR(cellfun('isempty',AR)) = []; %Remove empty cells

            %noftrials = size(E,2);
            noftrials_data2 = size(EEGaad_data2,2);
            
%             if noftrials < noftrials_data2
%                 EEGaad_data2(noftrials+1:end) = [];
%             elseif noftrials > noftrials_data2
%                 error('noftrials > noftrials_data2')
%                 break
%             end
            
            %noftrials = size(E,2);
            %noftrials_data2 = size(EEGaad_data2,2);
            
%             if noftrials ~= noftrials_data2
%                 error('Still a mistake')
%             elseif noftrials == noftrials_data2
%                 for i = 1: noftrials
%                     %sizes_E(i) = size(E{i},2);
%                     sizes_data2(i) = size(EEGaad_data2{i},2);
%                     if isequal(sizes_E, sizes_data2)
%                     else
%                         error('sizes_E and sizes_data2 are not equal')
%                     end
%                 end
%             end
% %             
             E_data2 = {};
            %FILTERING and RESAMPLING
            for trl = 1:noftrials_data2
                for ch = 1:24                  
                    %E{trl}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, E{trl}(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                    E_data2{trl}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, EEGaad_data2{trl}(ch,:));
                    %New(ch,:) = resample(double(E{trl}(ch,:)),NSR,EEGsampleRate);
                    New_data2(ch,:) = resample(double(E_data2{trl}(ch,:)),NSR,EEGsampleRate);
                end
                
                %E{trl} = New;
                E_data2{trl} = New_data2;
                
                %New = []; 
                New_data2 = [];                
                                   
                AL{trl} = Rbp([], HP, 6, fs, AL{trl}); % Band pass
                AR{trl} = Rbp([], HP, 6, fs, AR{trl}); % Band pass
            
                AL{trl} =  resample(double(AL{trl}),NSR,fs);
                AR{trl} =  resample(double(AR{trl}),NSR,fs);           
                
            end

            
            for trl = 1:noftrials_data2               
                % Lags                          
                %x = E{trl}';
                x_data2 = E_data2{trl}';
                start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                [start,fin] = deal(-fin,-start);
                nofsamples = size(x_data2,1);
                %[X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples
                [X_data2] = aad_LagGenerator(x_data2,start:fin); % Time lags of 2-500ms in samples
                
                %EE{trl} = X;
                EE_data2{trl} = X_data2;
                
                % Covariance Matrices
                %Rxx{trl} = (X'*X)/nofsamples;
                Rxx_data2{trl} = (X_data2'*X_data2)/nofsamples;
                
                if left
                %Rxy_at{trl} = (X'*squeeze(AL{trl})')/nofsamples;
                Rxy_at_data2{trl} = (X_data2'*squeeze(AL{trl})')/nofsamples;
                else
                %Rxy_at{trl} = (X'*squeeze(AR{trl})')/nofsamples;
                Rxy_at_data2{trl} = (X_data2'*squeeze(AR{trl})')/nofsamples;
                end
                
            end
            
            if left
                Attended = AL;
                Unattended = AR;
                
            else
                Attended = AR;
                Unattended = AL;
                
            end
                                       
            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data')
            if left
            %save(['Left_Session_' num2str(session) '_cut_AL.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2', 'Attended', 'Unattended') 
            save(['test_ses2_left2'], 'EE_data2', 'Rxx_data2', 'Rxy_at_data2', 'Attended', 'Unattended') 

            else
            save(['Right_Session_' num2str(session) '_cut_AL.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2', 'Attended', 'Unattended') 
            end
            
            clear EE Attended Unattended EE_data2 Attended_WIN
           
            end
        end

    end



% Session accuracy
Session_Accuracy = [];

            %cd('')
%load(['TEST_S4'])
%avg_dec4 = avg_dec;
%avg_dec = [];

        difference_self_data2 = [];
        %difference_self = [];
        %difference_dec4_data2 = [];
        %difference_dec4 = [];

%     for LR = 1 % Left - Right loop
%         if LR == 1
%             left = true;
%         else
%             left = false;
%         end
%     for session = 2
%                     if LR == 1 && session == 1
%             else

            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data')

            %if left
            load(['test_ses2_left2']) 
%             else
%             load(['Right_Session_' num2str(session) '_cut_AL.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2', 'Attended', 'Unattended') 
%             end
        
        noftrials = size(Attended,2);
        train = 1:noftrials;
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            
            trials2 = train(1:end-1);
            
            % Accuracy, decoder trained on self
            %RXX = zeros(size(Rxx{1})); RXY_ATT = zeros(size(Rxy_at{1}));
            RXX_data2 = zeros(size(Rxx_data2{1})); RXY_ATT_data2 = zeros(size(Rxy_at_data2{1}));
            for i = trials2
                %RXX = RXX + Rxx{i}; % plain averaging of cov matrices
                %RXY_ATT = RXY_ATT + Rxy_at{i};
                
                RXX_data2 = RXX_data2 + Rxx_data2{i}; % plain averaging of cov matrices
                RXY_ATT_data2 = RXY_ATT_data2 + Rxy_at_data2{i};
            end
            
%             RXX = RXX/(noftrials-1);
%             RXY_ATT = RXY_ATT/(noftrials-1);
            RXX_data2 = RXX_data2/(noftrials-1);
            RXY_ATT_data2 = RXY_ATT_data2/(noftrials-1);
        
%             avg_dec_self = RXX \ RXY_ATT; % LS solution.
%             recenv_self = EE{train(end)} * avg_dec_self;
%             CA_self(HH) = corr2(recenv_self,Attended{train(end)}'); % Attended
%             CUA_self(HH) = corr2(recenv_self,Unattended{train(end)}'); % Unattended
%             
            avg_dec_self_data2 = RXX_data2 \ RXY_ATT_data2; % LS solution.
            recenv_self_data2 = EE_data2{train(end)} * avg_dec_self_data2;
            CA_self_data2(HH) = corr2(recenv_self_data2,Attended{train(end)}'); % Attended
            CUA_self_data2(HH) = corr2(recenv_self_data2,Unattended{train(end)}'); % Unattended
            % OPM: Unattended_WIN bestaat niet... 
            %CA_self_data2(HH) = corr2(recenv_self_data2,Attended_WIN{train(end)}'); % Attended
            %CUA_self_data2(HH) = corr2(recenv_self_data2,Unattended_WIN{train(end)}'); % Unattended

            
%             recenv_dec4 = EE{train(end)} * avg_dec4;
%             CA_dec4(HH) = corr2(recenv_dec4,Attended{train(end)}'); % Attended
%             CUA_dec4(HH) = corr2(recenv_dec4,Unattended{train(end)}'); % Unattended
%             
%             recenv_dec4_data2 = EE_data2{train(end)} * avg_dec4;
%             CA_dec4_data2(HH) = corr2(recenv_dec4_data2,Attended{train(end)}'); % Attended
%             CUA_dec4_data2(HH) = corr2(recenv_dec4_data2,Unattended{train(end)}'); % Unattended            
%             
            train = [train(end) train(1:end-1)];
        end
        
%         Session_Accuracy_self(session, LR) = (sum(CA_self-CUA_self> 0)/(noftrials))*100 % Final accuracy on all trials
%         Session_Accuracy_dec4(session, LR) = (sum(CA_dec4-CUA_dec4> 0)/(noftrials))*100 % Final accuracy on all trials

        Session_Accuracy_self_data2(session, LR) = (sum(CA_self_data2-CUA_self_data2> 0)/(noftrials))*100 % Final accuracy on all trials
%         Session_Accuracy_dec4_data2(session, LR) = (sum(CA_dec4_data2-CUA_dec4_data2> 0)/(noftrials))*100 % Final accuracy on all trials
        
       
        difference_self_data2 = cat(2,difference_self_data2,CA_self_data2-CUA_self_data2);
%         difference_self = cat(2,difference_self,CA_self-CUA_self);
%         difference_dec4_data2 = cat(2,difference_dec4_data2,CA_dec4_data2-CUA_dec4_data2);
%         difference_dec4 = cat(2,difference_dec4,CA_dec4-CUA_dec4);
        
%         CA_self = []; CUA_self = [];CA_dec4 = []; CUA_dec4 = [];
        CA_self_data2 = []; CUA_self_data2 = [];
        %CA_dec4_data2 = []; CUA_dec4_data2 = [];

%                     end
%                     end
%     end

            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data')
        save(['test2.mat'], 'Session_Accuracy_self_data2', 'difference_self_data2')
%% Decoder over all sessions
EE_all = {}; Attended_all ={}; Unattended_all = {}; Rxx_all = {}; Rxy_all = {};
    for LR = 1:2 % Left - Right loop
        if LR == 1
            left = true;
        else
            left = false;
        end
    for session = 1:2
                    if LR == 1 && session == 1
            else

            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data')

            if left
            load(['Left_Session_' num2str(session) '_cut_AL.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2', 'Attended', 'Unattended') 
            else
            load(['Right_Session_' num2str(session) '_cut_AL.mat'], 'EE', 'EE_data2', 'Rxx', 'Rxx_data2', 'Rxy_at', 'Rxy_at_data2', 'Attended', 'Unattended') 
            end
            
            EE_all = [EE_all, EE];
            Attended_all = [Attended_all, Attended];
            Unattended_all = [Unattended_all, Unattended];
            Rxx_all = [Rxx_all, Rxx];
            Rxy_all = [Rxy_all, Rxy_at];
                    end
    end
    end
    
    % Decoder
            RXX = zeros(size(Rxx_all{1})); RXY_ATT = zeros(size(Rxy_all{1}));
            noftrials = size(Rxx_all,2);
            for i = noftrials
                RXX = RXX + Rxx_all{i}; % plain averaging of cov matrices
                RXY_ATT = RXY_ATT + Rxy_all{i};
            end
            
            RXX = RXX/(noftrials);
            RXY_ATT = RXY_ATT/(noftrials);
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            
            save(['decoder'], 'avg_dec')