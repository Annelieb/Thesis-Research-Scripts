% S9 compare on-line and mimick

%% Get the EEG data of online
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];
S = 9

EEG_online = []; EEG_mimick = []; At_online = []; At_mimick = []; Ut_online = []; Ut_mimick = [];

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2')
load('Accuracies_S 9_NSR20_HP8_LP1_onr2.mat', 'Attended', 'Unattended')
for k = 1:312
    At_mimick = [At_mimick, Attended{k}'];
    Ut_mimick = [Ut_mimick, Unattended{k}'];
end

for TF = 1:2
    if startleft(S) == 1
        ear = [1 1 0 0; 0 0 1 1];
    else
        ear = [0 0 1 1; 1 1 0 0];
    end
    for session = 1:4
        cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject9'])
        if ear(TF, session) == 1
            EAR = 'left';
        else
            EAR = 'right';
        end
        if TF == 1
            EEG = pop_biosig(['Training_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
            load(['Decoder_Params_Subject9_Session' num2str(session) '_' EAR '.mat'])
        else
            EEG = pop_biosig(['Feedback_Subject_' num2str(S) '_Session_' num2str(session) '_' EAR '.gdf']);
            load(['Feedback_Subject9_Session' num2str(session) '_' EAR '.mat'])
        end
        
        % GDF file processing
        locs(1) = EEG.event(1).latency;
        if session == 4
            locs(2) = locs(1) + ((395)*500);
        else
            locs(2) = EEG.event(2).latency; %Rob doet dit iets anders... => maakt geen verschil
        end
        EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
        
        WIN(cellfun('isempty',WIN)) = []; %Remove empty cells
        trn = 39;
        E = {};
        for i = 1:trn
            for ch = 1:24
                E{i}(ch,:) = EEGaad(ch,WIN{i}+1:WIN{i+1});
            end
        end
        
        % Put everything back together
        % EEG data & Attended & Unattended
        data2(1:2) = [];
        data2(40:end) = [];
        for l = 1:39
            EEG_online = [EEG_online, data2{l}(1:24,:)];
            EEG_mimick = [EEG_mimick, E{l}];
            At_online = [At_online, Attended{l}'];
            Ut_online = [Ut_online, Unattended{l}'];
        end
    end
end

%%

% Attended
figure
subplot(3,1,1)
plot(At_online, 'b')
xlim([0 62480])
xlabel('Data points')
ylabel('Envelope amplitude')
title('Online attended audio envelope')

subplot(3,1,2)
plot(At_mimick, 'r')
xlim([0 62480])
xlabel('Data points')
ylabel('Envelope amplitude')
title('Mimicked attended audio envelope')

subplot(3,1,3)
diffA = (At_mimick-At_online)./((At_mimick+At_online)/2);
plot(diffA, 'k')
xlim([0 62480])
xlabel('Data points')
ylabel('Relative difference [%]')
title('Relative difference')

suptitle(['Comparing the online and mimicked \newline attended audio envelope'])

% Unattended
figure
subplot(3,1,1)
plot(Ut_online)
xlim([0 62480])
xlabel('Data points')
ylabel('Envelope amplitude')
title('On-line Unattended audio envelope')

subplot(3,1,2)
plot(Ut_mimick)
xlim([0 62480])
xlabel('Data points')
ylabel('Envelope amplitude')
title('Mimiked Unattended audio envelope')

subplot(3,1,3)
diffU = (Ut_mimick-Ut_online)./((Ut_mimick+Ut_online)/2);
plot(diffU)
xlim([0 62480])
xlabel('Data points')
ylabel('Relative difference [%]')
title('Relative difference')

suptitle(['Comparing the on-line and mimicked \newline unattended audio envelope'])

% EEG
figure
subplot(3,1,1)
plot(EEG_online(17,:), 'b')
xlim([0 1560064])
xlabel('Data points')
ylabel('EEG amplitude')
title('Online EEG data')

subplot(3,1,2)
plot(EEG_mimick(17,:), 'r')
xlim([0 1560064])
xlabel('Data points')
ylabel('EEG amplitude')
title('Mimiked EEG data')

subplot(3,1,3)
diffE = (EEG_mimick(17,:)-EEG_online(17,:))./((EEG_mimick(17,:)+EEG_online(17,:))/2);
plot(diffE, 'k')
xlim([0 1560064])
xlabel('Data points')
ylabel('Relative difference [%]')
title('Relative difference')

suptitle('Comparing the online and mimicked EEG data')

