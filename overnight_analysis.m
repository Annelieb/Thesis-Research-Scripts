%% Initial parameters
close all 
clear all

subjects = 1:12;
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];


%%
% Get the covariances
envelopemethod = 'powerlaw'; %or hilbert
power = 0.6;

for onr = 10 % veranderde naar trll alleen 10 30 en 60
    if onr == 1 % normal
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 20 30 40 50 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 2 % mimmick online
        mimmick_online = true;
        trll = 10
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
         seperate_lags = false;
    elseif onr == 3 % audioLP
        lag1 = 100;
        lag2 = 400;                                                   
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 20 30 40 50 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = true;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 4 % LP & HP
        lag1 = 100;
        lag2 = 400;                                               
        lp = [0.5 1 2]; % in pilot 2 tot 9 Hz bekeken
        hp = [7 8 9];
        nsr = 20
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 5 % NSR at lags 100 to 500 
        lag1 = [0];
        lag2 = [500];
        nsr = [20]; 
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 6 % subband
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 30 60];
        subband = true;
         if subband
            spacing = 1.5; % OPM: How is this defined? Ask Neetha...
            gammatonefreqs = erbspacebw(150,4000,spacing); % gammatone filter centerfrequencies % OPM: What is all off this? How are the inputs defined? Ask Neetha...
            gammatonebetamul = spacing*ones(size(gammatonefreqs)); % multiplier for gammatone filter bandwidths
            cut_off = 4000;
         end 
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 7 % cut first
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = true;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 8 % start online
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = true;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 9 % cut first & start online
        lag1 = 100;
        lag2 = 400;                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = true;
        audioLP = false;
        start_online = true;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 10 % lags
        lag1 = [150];
        lag2 = [350 400 450 500];                                                    
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        nsr = 20
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
         seperate_lags = false;
    elseif onr == 11 % seperate lags
        lag1 = [0];
        lag2 = [500];
        nsr = [10 20 30 40]; % 20 hoeft niet gedaan te worden, want die heb ik. 10 is gedaan voor gestopt, dus dat er even uitgehaald...
        lp = 1; % in pilot 2 tot 9 Hz bekeken
        hp = 8;
        trll = [10 30 60];
        subband = false;
        intSamp = 500;
        cutfirst = false;
        audioLP = false;
        start_online = false;
        mimmick_online = false;
        seperate_lags = true;
    end

for l1 = lag1
for l2 = lag2
%     if l1 == 150 && l2 == 350
%         subjects = 7:12;
%     else 
%         subjects = 1:12;
%     end
for LP = lp    
for HP = hp
for NSR = nsr
intermediateSampleRate = intSamp

for S = subjects % Subject loop
        EE = {}; Attended = {}; Unattended = {}; Rxx_all ={}; Rxy_at_all = {};
    
%     if onr == 4 && LP == 0.5 && S <8
%     elseif onr == 4 && LP == 0.5 && S > 7 && HP < 9
%     else        
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
                        envelopeL = (sum(envelopeL,2));
                        envelopeR = (sum(envelopeR,2));
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
                    start = floor(l1/1e3*NSR); %convert milliseconds to samples 100
                    fin = ceil(l2/1e3*NSR); %convert milliseconds to samples 400
                    nofsamples = size(x,1);
                    [start,fin] = deal(-fin,-start);
                    [X] = aad_LagGenerator(x,start:fin); % Time lags of 2-500ms in samples

                    test{trials} = X;
                    EE{TRL,gentel6} = X;

                    % Covariance Matrices
                    if seperate_lags
                        noflag = length(start:fin);
                        for l = 1:noflag
                            Rxx(trials,l,:,:) = (X(:, (l-1)*24 +1 : l*24)'*X(:, (l-1)*24 +1 : l*24))/nofsamples;
                            Rxx_all{TRL}(gentel6,l,:,:) = Rxx(trials,l,:,:);
                            
                            if left
                                Rxy_at(trials,l,:) = (X(:, (l-1)*24 +1 : l*24)'*AL{trials})/nofsamples;
                                At = AL;
                                Unat = AR;
                            else
                                Rxy_at(trials,l,:) = (X(:, (l-1)*24 +1 : l*24)'*AR{trials})/nofsamples;
                                At = AR; 
                                Unat = AL; 
                            end
                            Rxy_at_all{TRL}(gentel6,l,:) = Rxy_at(trials,l,:);
                        end
                            
                    else
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
                        Rxy_at_all{TRL}(gentel6,:) = Rxy_at(trials,:);
                    end
                end
                
                check = 4
                
                if seperate_lags
                    train = (1:trn); CA = []; CUA = [];
                    for l = 1:noflag
                        for HH = 1:trn % repeat for all trials - leave one out
                            clear RXX RXY_ATT
                            trials2 = train(1:end-1);

                            RXX = squeeze(mean(Rxx(trials2,l,:,:),1)); % plain averaging of cov matrices
                            RXY_ATT = squeeze(mean(Rxy_at(trials2,l,:),1));

                            avg_dec_ses{TRL, session,l, train(end)} = RXX \ RXY_ATT; % LS solution. 

                            recenv = test{train(end)}(:,(l-1)*24 +1 : l*24) * avg_dec_ses{TRL, session,l, train(end)};

                            CA(train(end)) = corr2(recenv,At{train(end)}); % Attended
                            CUA(train(end)) = corr2(recenv,Unat{train(end)}); % Unattended
                            train = [train(end) train(1:end-1)];
                        end
                        session_CA{TRL, TF, session,l} = CA;
                        session_CUA{TRL, TF, session,l} = CUA;
                        session_difference{TRL, TF, session,l} = CA - CUA;
                        session_Accuracy(TRL, session, TF,l) = (sum(CA-CUA> 0)/trn)*100 
                    end
                else
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
            end
            check = 6
        end
        
        check = 7
        
        if seperate_lags
            %Subject Accuracy
            clear CA CUA train
            noftrials = gentel6;
            train = (1:noftrials);
            for l = 1:noflag
                for HH = 1:noftrials % repeat for all trials - leave one out
                    clear RXX RXY_ATT recenv
                    trials2 = train(1:end-1);

                    RXX = squeeze(mean(Rxx_all{TRL}(trials2,l,:,:),1)); % plain averaging of cov matrices
                    RXY_ATT = squeeze(mean(Rxy_at_all{TRL}(trials2,l,:),1));

                    avg_dec_subj{TRL,l, train(end)} = RXX \ RXY_ATT; % LS solution.
                    recenv = EE{TRL, train(end)}(:,(l-1)*24 +1 : l*24) * avg_dec_subj{TRL,l, train(end)};

                    CA(train(end)) = corr2(recenv,Attended{TRL,train(end)}); % Attended
                    CUA(train(end)) = corr2(recenv,Unattended{TRL, train(end)}); % Unattended
                    train = [train(end) train(1:end-1)];
                end
                subject_CA{TRL,l} = CA;
                subject_CUA{TRL,l} = CUA;
                subject_difference{TRL,l} = CA - CUA;
                subject_Accuracy(TRL,l) = (sum(CA-CUA> 0)/noftrials)*100
            end
        else

            %Subject Accuracy % Note not assessed for mimick online! 
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
    end
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
    if onr == 10
        save(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_lag1_' num2str(l1) '_lag2_' num2str(l2) '_onr' num2str(onr) '.mat'], 'EE', 'Rxx_all', 'Rxy_at_all', 'Attended', 'Unattended', 'subject_CA', 'subject_CUA', 'subject_difference', 'subject_Accuracy', 'session_CA', 'session_CUA', 'session_difference', 'session_Accuracy', 'avg_dec_subj', 'avg_dec_ses' )
    else
        save(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'EE', 'Rxx_all', 'Rxy_at_all', 'Attended', 'Unattended', 'subject_CA', 'subject_CUA', 'subject_difference', 'subject_Accuracy', 'session_CA', 'session_CUA', 'session_difference', 'session_Accuracy', 'avg_dec_subj', 'avg_dec_ses' )
    end
    clear EE Attended Unattended Rxx_all Rxy_at_all subject_CA subject_CUA subject_difference subject_Accuracy session_CA session_CUA session_difference session_Accuracy avg_dec_subj avg_dec_ses
    %end
end
end
end
end
end
end
end