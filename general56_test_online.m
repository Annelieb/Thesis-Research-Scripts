%% Initial parameters
close all 
clear all

subjects = 1:2;
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

lag = [100 400];                                                    
LP = 1; % in pilot 2 tot 9 Hz bekeken
HP = 8;
%LP = [0, 1, 2, 3, 4];
%HP = [6,7, 8, 9, 10];
%NSR = 20;
NewSamplingRate = 20 %[10 15 20 25 30 35 40 45 50 80]; % 200 500]; % New sampling rate; 200 and 500 take to much computing time
trll = [10 30 60]; %[10 20 30 40 50 60 5 15]; % in seconds

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?

 subband = false; 
 if subband
    spacing = 1.5; % OPM: How is this defined? Ask Neetha...
    gammatonefreqs = erbspacebw(150,4000,spacing); % gammatone filter centerfrequencies % OPM: What is all off this? How are the inputs defined? Ask Neetha...
    gammatonebetamul = spacing*ones(size(gammatonefreqs)); % multiplier for gammatone filter bandwidths
 end
 
intSamp = 500; %[20 70 120 160 200 300 400 500]; % Why 120?
cut_off = 4000; %cut-off frequency of headphones. OPM: Find what this is in our experiment.

cutfirst = false;
audioLP = false;

start_online = false;

mimmick_online = false;
if mimmick_online
    trll = 10
    lag = [100 400];                                                    
    LP = 1; % in pilot 2 tot 9 Hz bekeken
    HP = 8;
    NewSamplingRate = 20
    subband = false;
    intSamp = 500;
    cutfirst = false;
    audioLP = false;
    start_online = false;
end

%% Get the covariances
% for LP = 2
% for HP = 9
% for nsr = 3:3
% NSR = NewSamplingRate(nsr)
NSR = NewSamplingRate
% for isr = 8:8

%intermediateSampleRate = intSamp(isr)
intermediateSampleRate = intSamp

for S = subjects % Subject loop
        EE = {}; Attended = {}; Unattended = {}; Rxx_all ={}; Rxy_at_all = {};
    
    for TRL = 1:length(trll) % Loop for trial length

        triallength = trll(TRL)
        X=[];
        gentel4 = 0;
        gentel5 = 0; gentel6 = 0;
        ses = 0

        for TF = 1:2 % Left - Right loop

            for session = 1:4 % 4 right ear session - Train on each one separate
                data = [];
                [data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
                [data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
                data = cat(2,data2(:,1),data3(:,1));

                ses = ses+1

                if startleft(S) == 1
                    ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
                else
                    ear = [0 0 1 1; 1 1 0 0];
                end
                if ear(TF, session) == 1
                    EAR = 'left'; 
                    left = true;
                else
                    EAR = 'right';
                    left = false;
                end
                cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])

                if TF == 1
                    EEG = pop_biosig(['Training_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
                    if start_online || mimmick_online
                        load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'WIN')
                    end
                elseif TF == 2
                    EEG = pop_biosig(['Feedback_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
                    if start_online || mimmick_online
                       load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'WIN')
                    end
                end

                EEGsampleRate = 500;
                              
                locs(1) = EEG.event(1).latency;
                if session == 4
                    locs(2) = locs(1) + ((395)*EEGsampleRate);
                else
                    locs(2) = EEG.event(2).latency; %Rob doet dit iets anders... => maakt geen verschil
                end
                
                if start_online
                    WIN(cellfun('isempty',WIN)) = []; %Remove empty cells
                    locs(1) = locs(1) + WIN{1};
                    data(1:round(WIN{1}/500*fs),:) = [];
                end

                EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG

                if session < 3
                    AudioL = data(:,1); AudioR = data(:,2);
                elseif session > 2
                    AudioR = data(:,1); AudioL = data(:,2);
                end
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

                AudioL = []; AudioR = []; % save RAM

                envelopeL = resample(envelopeL,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                envelopeR = resample(envelopeR,intermediateSampleRate,fs); % OPM: Is dit nodig in ons script?
                EEGaad = (resample((double(EEGaad))', intermediateSampleRate, EEGsampleRate))';
                EEGsampleRate = intermediateSampleRate;
                fs = intermediateSampleRate; 

                % FILTERING

                if cutfirst || mimmick_online
                    % CREATE TRIALS
                    if mimmick_online
                        WIN(cellfun('isempty',WIN)) = []; %Remove empty cells
                        trn = 39;
                        E = {}; AL = {}; AR = {};
                        for i = 1:trn
                            for ch = 1:24
                                E{i}(ch,:) = EEGaad(ch,WIN{i}+1:WIN{i+1});
                            end
                            AL{i} = envelopeL(WIN{i}+1:WIN{i+1});
                            AR{i} = envelopeR(WIN{i}+1:WIN{i+1});
                        end

                    else
                        window = triallength*EEGsampleRate; % seconds * Sampling-rate
                        [trs] = mod(size(EEGaad,2),window);
                        trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
                        
                        s1 = 0; E = {}; AL = {}; AR = {};
                        for i = 1:trn
                            for ch = 1:24
                                E{i}(ch,:) = EEGaad(ch,s1+1:s1+window);
                            end
                            if subband
                                nofsubbands = length(gammatonefreqs);
                                for sb = 1:nofsubbands
                                    AL{i}(sb,:) = envelopeL(s1+1:s1+window,sb);
                                    AR{i}(sb,:) = envelopeR(s1+1:s1+window,sb);
                                end
                            else
                                AL{i} = envelopeL(s1+1:s1+window);
                                AR{i} = envelopeR(s1+1:s1+window);
                            end
                            s1 =s1+window; 
                        end
                    end

                    % PROCESSING
                    New = {};
                    for i = 1:trn
                        for ch = 1:24
                            E{i}(ch,:) = Rbp(LP, HP, 6, EEGsampleRate, squeeze(E{i}(ch,:))); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                            New{i}(ch,:) = resample(double(E{i}(ch,:)),NSR,EEGsampleRate);
                        end

                        % AUDIO
                        gentel5 = gentel5+1;
                        if subband
                            nofsubbands = length(gammatonefreqs);
                            for sb = 1:nofsubbands
                                if audioLP
                                    AL{i}(sb,:) = Rbp(LP, HP, 6, fs, AL{i}(sb,:));
                                    AR{i}(sb,:) = Rbp(LP, HP, 6, fs, AR{i}(sb,:));
                                else
                                    AL{i}(sb,:) = Rbp([], HP, 6, fs, AL{i}(sb,:));
                                    AR{i}(sb,:) = Rbp([], HP, 6, fs, AR{i}(sb,:));
                                end
                            end
                            % Sum of subbands so that there only is one envelope again
                            AL_new{i} = squeeze(sum(AL{i},1));
                            AR_new{i} = squeeze(sum(AR{i},2));
                            AL = AL_new;
                            AR = AR_new;
                        else
                            if audioLP
                                AL{i} = Rbp(LP, HP, 6, fs, AL{i}); % Band pass
                                AR{i} = Rbp(LP, HP, 6, fs, AR{i}); % Band pass
                            else
                                AL{i} = Rbp([], HP, 6, fs, AL{i}); % Band pass
                                AR{i} = Rbp([], HP, 6, fs, AR{i}); % Band pass
                            end
                        end

                        NewAL{i} =  (resample(double(AL{i}),NSR,fs))';
                        NewAR{i} =  (resample(double(AR{i}),NSR,fs))';

                        if left
                            Attended{TRL,gentel5} = NewAL{i};
                            Unattended{TRL, gentel5} = NewAR{i};
                        else
                            Attended{TRL, gentel5} = NewAR{i};
                            Unattended{TRL, gentel5} = NewAL{i};
                        end
                    end

                    E = New; AL = NewAL; AR = NewAR;          
                    New = {}; NewAL = {}; NewAR = {};
                else
                    New2 = [];
                    for ch = 1:24 % EEG
                        [EEGaad(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                        New2(ch,:) = resample(double(EEGaad(ch,:)),NSR,EEGsampleRate);
                    end
                    EEGaad = New2;
                    EEGsampleRate = NSR;
                    New2 = [];

                    if subband
                        nofsubbands = length(gammatonefreqs);
                        for sb = 1:nofsubbands
                            if audioLP
                                [envelopeL(:,sb)] = Rbp(LP, HP, 6, fs, envelopeL(:,sb)); % Band pass
                                [envelopeR(:,sb)] = Rbp(LP, HP, 6, fs, envelopeR(:,sb)); % Band pass
                            else
                                [envelopeL(:,sb)] = Rbp([], HP, 6, fs, envelopeL(:,sb)); % Band pass
                                [envelopeR(:,sb)] = Rbp([], HP, 6, fs, envelopeR(:,sb)); % Band pass                
                            end
                        end
                        % Sum of subbands so that there only is one envelope again
                        envelopeL = (sum(envelopeL,2))';
                        envelopeR = (sum(envelopeR,2))';
                    else
                        if audioLP
                            [envelopeL(:)] = Rbp(LP, HP, 6, fs, envelopeL(:,1)); % Band pass
                            [envelopeR(:)] = Rbp(LP, HP, 6, fs, envelopeR(:,1)); % Band pass
                        else
                            [envelopeL(:)] = Rbp([], HP, 6, fs, envelopeL(:,1)); % Band pass
                            [envelopeR(:)] = Rbp([], HP, 6, fs, envelopeR(:,1)); % Band pass                
                        end
                    end
                    check = 1;

                    envelopeL =  resample(double(envelopeL),NSR,fs);
                    envelopeR =  resample(double(envelopeR),NSR,fs);

                    fs = NSR;

                    % CREATE TRIALS
                    window = triallength*NSR; % seconds * Sampling-rate
                    [trs] = mod(size(EEGaad,2),window);
                    trn = (size(EEGaad,2)-trs)/window; % number of trials possible.

                    % EEG
                    s1 = 0; E = []; AL = []; AR = [];
                    for i = 1:trn
                        for ch = 1:24
                            E{i}(ch,:) = EEGaad(ch,s1+1:s1+window);
                        end
                        % AUDIO
                        AL{i} = envelopeL(s1+1:s1+window);
                        AR{i} = envelopeR(s1+1:s1+window);
                        gentel5 = gentel5+1;
                        if left
                            Attended{TRL,gentel5} = AL{i};
                            Unattended{TRL,gentel5} = AR{i};
                        else
                            Attended{TRL,gentel5} = AR{i};
                            Unattended{TRL,gentel5} = AL{i};
                        end
                        s1 = s1+window;
                        gentel4 = gentel4+1;
                    end
                    check = 2
                end
                
                check = 3

                % Lags
                test = {}; Rxx = []; Rxy_at = []; x = [];
                for trials = 1:trn

                    gentel6 = gentel6 + 1
                    x=squeeze(E{trials})';
                    start = floor(lag(1)/1e3*NSR); %convert milliseconds to samples 100
                    fin = ceil(lag(2)/1e3*NSR); %convert milliseconds to samples 400
                    nofsamples = size(x,1);
                    [start,fin] = deal(-fin,-start);
                    [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples

                    test{trials} = X;
                    EE{TRL,gentel6} = X;

                    % Covariance Matrices
                    Rxx(trials,:,:) = (X'*X)/nofsamples;
                    Rxx_all{TRL}(gentel6,:,:) = Rxx(trials,:,:);
                    if left
                        Rxy_at(trials,:) = (X'*AL{trials})/nofsamples;
                        At = AL;
                        Unat = AR;
                    else
                        Rxy_at(trials,:) = (X'*AR{trials})/nofsamples;
                        At = AR; 
                        Unat = AL; 
                    end
                    Rxy_at_all{TRL}(gentel6,:) = Rxy_at(trials,:)
                end
                
                check = 4

                % Session Accuracy
                train = (1:trn); CA = []; CUA = [];
                for HH = 1:trn % repeat for all trials - leave one out
                    clear RXX RXY_ATT
                    trials2 = train(1:end-1);

                    RXX = squeeze(mean(Rxx(trials2,:,:),1)); % plain averaging of cov matrices
                    RXY_ATT = squeeze(mean(Rxy_at(trials2,:),1));
                    
                    if mimmick_online == 1 && TF == 2
                        load(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec')
                        avg_dec_ses{TRL, session, train(end)} = avg_dec;
                    else
                        avg_dec_ses{TRL, session, train(end)} = RXX \ RXY_ATT'; % LS solution. 
                    end
                    
                    recenv = test{train(end)} * avg_dec_ses{TRL, session, train(end)};

                    CA(train(end)) = corr2(recenv,At{train(end)}); % Attended
                    CUA(train(end)) = corr2(recenv,Unat{train(end)}); % Unattended
                    train = [train(end) train(1:end-1)];
                end
                session_CA{TRL, TF, session} = CA;
                session_CUA{TRL, TF, session} = CUA;
                session_difference{TRL, TF, session} = CA - CUA;
                session_Accuracy(TRL, session, TF) = (sum(CA-CUA> 0)/trn)*100            
            end
            check = 6
        end
        
        check = 7

        %Subject Accuracy
        clear CA CUA train
        noftrials = gentel6;
        train = (1:noftrials);
        for HH = 1:noftrials % repeat for all trials - leave one out
            clear RXX RXY_ATT recenv
            trials2 = train(1:end-1);

            RXX = squeeze(mean(Rxx_all{TRL}(trials2,:,:),1)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(Rxy_at_all{TRL}(trials2,:),1));

            avg_dec_subj{TRL, train(end)} = RXX \ RXY_ATT'; % LS solution.
            recenv = EE{TRL, train(end)} * avg_dec_subj{TRL, train(end)};

            CA(train(end)) = corr2(recenv,Attended{TRL,train(end)}); % Attended
            CUA(train(end)) = corr2(recenv,Unattended{TRL, train(end)}); % Unattended
            train = [train(end) train(1:end-1)];
        end
        subject_CA{TRL} = CA;
        subject_CUA{TRL} = CUA;
        subject_difference{TRL} = CA - CUA;
        subject_Accuracy(TRL) = (sum(CA-CUA> 0)/noftrials)*100            

    end


    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')

    save(['Accuracies_S ' num2str(S) '_normal.mat'], 'EE', 'Rxx_all', 'Rxy_at_all', 'Attended', 'Unattended', 'subject_CA', 'subject_CUA', 'subject_difference', 'subject_Accuracy', 'session_CA', 'session_CUA', 'session_difference', 'session_Accuracy', 'avg_dec_subj', 'avg_dec_ses' ) % , '-v7.3') Vanaf NSR > 50 nodig

    clear EE Attended Unattended

end
% end
% end
% end
% end

%% Subject Accuracies samenvoegen
subjects = 1:7

for S = subjects
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
    load(['Accuracies_S ' num2str(S) '_cut_first.mat'], 'subject_Accuracy')
    Subject_Accuracy(S, :) = subject_Accuracy
end

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
save(['Subject_Accuracy_' mat2str(subjects) '_cut_first.mat'], 'Subject_Accuracy')
%% Feedback accuracy chosen decoder compared to leave one out
subjects = 1:10
trll = [10 30 60];

for S = subjects
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
    load(['Accuracies_S ' num2str(S) '.mat'], 'EE', 'Attended', 'Unattended', 'subject_difference', 'subject_CA', 'subject_CUA')
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec')
    for TRL = 1:length(trll)
        nofFB = size(Attended{TRL},1)/2;
        train = (nofFB +1 :nofFB*2);
        CA = []; CUA =[];
        for HH = 1: nofFB
            
            recenv = EE{TRL}(:,:,train(end)) * avg_dec;

            CA(train(end)-nofFB) = corr2(recenv,Attended{TRL}(train(end),:)'); % Attended
            CUA(train(end)-nofFB) = corr2(recenv,Unattended{TRL}(train(end),:)'); % Unattended
            
            train = [train(end) train(1:end-1)];
        end
        % 'best decoder'
        FB_CA{TRL}(1,:) = CA;
        FB_CUA{TRL}(1,:) = CUA;
        FB_difference{TRL}(1,:) = CA - CUA;
        FB_Accuracy(TRL,1) = (sum(CA-CUA> 0)/nofFB)*100
        
        % leave one out
        FB_CA{TRL}(2,:) = subject_CA{TRL}(nofFB +1 :nofFB*2);
        FB_CUA{TRL}(2,:) = subject_CUA{TRL}(nofFB +1 :nofFB*2);
        FB_difference{TRL}(2,:) = subject_difference{TRL}(nofFB +1 :nofFB*2);;
        FB_Accuracy(TRL,2) = (sum(FB_difference{TRL}(2,:)> 0)/nofFB)*100
    end
    
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
    save(['FB_Accuracies_S' num2str(S) '.mat'], 'FB_CUA', 'FB_CA', 'FB_difference', 'FB_Accuracy')

end

%% Feedback Accuracies samenvoegen
subjects = 1:10
for S = subjects
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
    load(['FB_Accuracies_S' num2str(S) '.mat'], 'FB_Accuracy')
    FB_Accuracy_All(S, :,:) = FB_Accuracy
end

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\gdf accuracies')
save(['FB_Accuracy_All_' mat2str(subjects) '.mat'], 'FB_Accuracy_All')

%% Subject accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Test_general')
Subject_Accuracy = [];
nofsubj = 1;
%for LP = [0 1 2 3 4]
    %for HP = [6 7 8 9 10]
%for nsr = 3:3
%for isr = 8:8
for TRL = 1:3
    
    %load(['LP' num2str(LP) '_HP' num2str(HP) '_TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    %load(['TRL' num2str(TRL) '_ISR' num2str(isr) '_NSR' num2str(nsr) '.mat'])
    load(['General_S17' num2str(TRL) '_Subjects12.mat'])

    for S = 1:nofsubj
        
        noftrials = size(Attended{S},1);
        train = (1:noftrials);
    
        for HH = 1:noftrials % repeat for all trials - leave one out
            clear RXX RXY_ATT avg_dec
            trials2 = train(1:end-1);
        
            RXX = squeeze(mean(covar.Rxx_all(S,trials2,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,trials2,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            recenv = EE{S}(:,:,train(end)) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{S}(train(end),:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{S}(train(end),:)'); % Unattended
            train = [train(end) train(1:end-1)];
        end
         difference{TRL}(S,:) = CA - CUA;
        
        Subject_Accuracy(TRL,S) = (sum(CA-CUA> 0)/(noftrials))*100; % Final accuracy on all trials
        CA = []; CUA = [];
        
%         noftrlses = noftrials/8;
%         for i = 1:8
%         Session_Accuracy(S,i) = sum(difference(S,noftrlses*(i-1)+1:noftrlses*i)>0)/noftrlses*100;
%         end
        

    end

end
% end
% end
Subject_Accuracy
%Session_Accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Test_general')
        save(['Result_sessions_locs2anders.mat'], 'Subject_Accuracy', 'difference') %, 'Session_Accuracy')

        %save(['Result_LP' num2str(LP) '_HP' num2str(HP) '.mat'], 'Subject_Accuracy')
%    end
%end

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

