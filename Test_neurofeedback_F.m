%% 1a)
close all
clear all
%% 1b) Settings % These should be the same for all subjects and experiments!

NSR = 20;
LP = 1; % According to my findings setting LP = 1 yields better results; The decoder must be adjusted accordingly!
HP = 8; 
lags = [100 400]; % Set the lag values. Perhaps to 0-500 or 50-500 or 50-400.

bkrcolor = [0.2 0.2 0.2]; % Backround color


%% 2) Fill out subject, session and ear

S = 1; % Which subject? Change for names!
sess = 4; % Which session is the subject listening to? (1, 2, 3 or 4)
ear = 'right'; % Which ear is the subject attending to? ('left' or 'right')
CD = ['C:\Users\Rob XPS\Desktop\THESIS\Experiments\Subject' num2str(S)]; %Fill in saving directory
save(['Settings'], 'NSR', 'LP', 'HP', 'lags', 'bkrcolor', 'CD')

%% 3) Sessions to build decoder; no neurofeedback
load(['Settings'], 'NSR', 'LP', 'HP', 'lags', 'bkrcolor', 'CD')

if strcmp(ear, 'left')
    left = true; 
elseif strcmp(ear, 'right')
    left = false;
end

% instantiate the library

START = false;
GO = false;
stopp = false;

trial = 1;

% Loading and making the audio envelopes
data22 = [];
data = [];
[data22,fs] = audioread(['C:\Users\Rob XPS\Desktop\THESIS\Thesis Student 2015\Data\part' num2str(sess) '_track1_dry.wav']);
[data33,fs] = audioread(['C:\Users\Rob XPS\Desktop\THESIS\Thesis Student 2015\Data\part' num2str(sess) '_track2_dry.wav']);
data(:,1) = data22(:,1); % Left audio dry
data(:,2) = data33(:,1); % Right audio dry
if sess < 3
    AudioL = data(:,1); 
    AudioR = data(:,2);
elseif sess > 2
    AudioR = data(:,1); 
    AudioL = data(:,2);
end
envelopeL = abs(AudioL).^0.6; % Powerlaw Envelope
envelopeR = abs(AudioR).^0.6; % Powerlaw 0.6
envelopeL =  resample(double(envelopeL),500,fs);
envelopeR =  resample(double(envelopeR),500,fs);

%figure
x = [0 0 1 1 0];
y = [0 1 1 0 0];
hfig = figure('name','Performance','numbertitle','off','unit','normalized','outerposition',[0 0 1 1],'Color',bkrcolor,'menubar','none');
rectangle('Position', [0 0 1 1], 'FaceColor', bkrcolor, 'EdgeColor', bkrcolor)
ax = gca;
set(ax,'XTick',[]),set(ax,'YTick',[])
set(ax,'XColor',bkrcolor) , set(ax,'YColor',bkrcolor)
axis square
kleur = [0.4 0.4 0.4];
rectangle('Position', [0.4 0.4 0.2 0.2], 'Curvature', [1 1], 'FaceColor', kleur);


% loading the lsl library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');

result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','Stimulations'); 
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

tmr = timer('ExecutionMode', 'FixedRate', ...
    'Period', 10, ...
    'TimerFcn',...
'[chunk,stamps] = inlet.pull_chunk();data2{trial} = chunk;if trial > 1;if max(data2{trial}(25,:)) > 20;START = true;win1 = 1;[pks,locs] = findpeaks(double(data2{trial}(25,:)));if stopp;stop(tmr);clc;close;GO = false;end;end;if GO;clear New New2 AL AR envelopeL2 envelopeR2;for ch = 1:24;New(ch,:) = double(data2{trial}(ch,:));[New(ch,:)] = Rbp(LP, HP, 10, 500, New(ch,:));New2(:,ch) = resample(double(New(ch,:)),NSR,500);end;envelopeLq = envelopeL(win1:win1+size(data2{trial},2)-1);envelopeRq = envelopeR(win1:win1+size(data2{trial},2)-1);envelopeL2(1,:) =  double(envelopeLq);envelopeR2(1,:) =  double(envelopeRq);[AL(:,1)] = Rbp([], HP,10, 500, envelopeL2(1,:));[AR(:,1)] = Rbp([], HP,10, 500, envelopeR2(1,:));AL =  resample(double(AL),NSR,500);AR =  resample(double(AR),NSR,500);win1 = win1+size(data2{trial},2);WIN{trial} = win1;x=New2(:,:);start = floor(lags(1)/1e3*NSR);fin = ceil(lags(2)/1e3*NSR);nofsamples = size(x,1);[start,fin] = deal(-fin,-start);[X] = aad_LagGenerator(x,start:fin);EE{trial} = X;Rxx(trial,:,:)= (X''*X)/nofsamples;if left;Rxy_at(trial,:) = (X''*AL)/nofsamples;Attended{trial} = AL;Unattended{trial} = AR;else;Rxy_at(trial,:) = (X''*AR)/nofsamples;Attended{trial} = AR;Unattended{trial} = AL;end;end;if START;GO = true;if GO == 1&& stopp == 0;begin = trial;end;win1 = round(win1+size(data2{trial},2)-locs);START = false;stopp = true;WIN{trial} = win1;end;end;trial = trial +1;')

start(tmr);

%% stop tmr and save raw
stop(tmr)

cd(CD)
save(['Training_Raw_data_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'data2', 'Rxx', 'Rxy_at', 'Attended', 'Unattended', 'EE', 'WIN')
%% 4) Save the parameters to build decoder
% Some processing before saving
% Delete the empty cells.
EE(cellfun('isempty',EE)) = [];
Attended(cellfun('isempty', Attended)) = [];
Unattended(cellfun('isempty', Unattended)) = [];

Rxx(1:begin,:,:) = [];
Rxy_at(1:begin, :) = [];

if size(EE,2) > 39
    EE(40:end) = [];
    Attended(40:end) = [];
    Unattended(40:end) = [];
    Rxx(40:end,:,:) = [];
    Rxy_at(40:end,:) = [];
end

% Rxx(cellfun('isempty',Rxx)) = [];
% Rxy_at(cellfun('isempty',Rxy_at)) = [];
% 
% % Make Rxx and Rxy_at matrices
% Rxx_new = []; Rxy_new = [];
% noftrials = size(Rxy_at, 2);
% for i = 1:noftrials
%     Rxx_new = cat(3, Rxx_new, Rxx{i});
%     Rxy_new = cat(2, Rxy_new, Rxy_at{i});
% end
% clear Rxx Rxy_at
% 
% Rxx = Rxx_new; Rxy_at = Rxy_new; 
% 
% clear Rxx_new Rxy_new
% Save
cd(CD)
save(['Decoder_Params_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'data2', 'Rxx', 'Rxy_at', 'Attended', 'Unattended', 'EE', 'WIN')

%% 4b) [Optional] Check the session accuracy
cd(CD)
load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'])
noftrials = size(Attended,2);

train = (1:noftrials);


for HH = 1:noftrials % repeat for all trials - leave one out
    trials2 = train(1:end-1);
        
    RXX = squeeze(mean(Rxx(trials2,:,:),1)); % plain averaging of cov matrices
    RXY_ATT = squeeze(mean(Rxy_at(trials2,:),1))';
        
    avg_dec = RXX \ RXY_ATT; % LS solution.
    recenv = EE{train(end)} * avg_dec;
        
    CA(train(end)) = corr2(recenv,Attended{train(end)}); % Attended
    CUA(train(end)) = corr2(recenv,Unattended{train(end)}); % Unattended
    train = [train(end) train(1:end-1)];
end

difference = CA - CUA;
Session_Accuracy = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials

save(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'Session_Accuracy', 'difference','NSR', 'LP', 'HP', 'lags', 'CUA', 'CA')

%% Clear and close all and Restart from point 2 for more training sessions. Continue to point 5 to build the decoder. 

clear all 
close all

%% 5a) Build the decoder: Which sessions have been used?
S = 12; % Subject
sessions = [1 2 7 8]; % 1 2 3 4 = left sessions; 5 6 7 8 = right sessions
startleft = true; % First session was left -> true; or right -> false
%% 5b) Build all possible combinations of decoders and test them through leave one out
for nofsessions = 1:length(sessions) % loop over number of sessions that can be used to build the decoder
    combinations = nchoosek(sessions, nofsessions);
    nofcomb = size(combinations,1);
    for comb = 1:nofcomb
        difference_self =[]; difference_other =[];
        Rxx_dec = []; Rxy_dec = []; Attended_dec = {}; Unattended_dec = {}; EE_dec = {};
        Attended_other = {}; Unattended_other = {}; EE_other = {};
        
        others = sessions; 
        
        for decsession = combinations(comb, :)          
            others(find(others == decsession)) = [];
            if startleft
                if decsession < 5
                    sess = decsession;
                    ear = 'left';
                elseif decsession >4
                    ear = 'right';
                    sess = decsession -4;
                end
            else
                if decsession < 5
                    sess = decsession;
                    ear = 'right';
                elseif decsession >4
                    ear = 'left';
                    sess = decsession -4;
                end
            end
            load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'Rxx', 'Rxy_at', 'Attended', 'Unattended', 'EE')
            Rxx_dec = cat(1, Rxx_dec, Rxx(1:39,:,:)); 
            Rxy_dec = cat(1, Rxy_dec, Rxy_at(1:39,:));
            Attended_dec = [Attended_dec, Attended{1:39}];
            Unattended_dec = [Unattended_dec, Unattended{1:39}];
            EE_dec = [EE_dec, EE{1:39}];
            clear Rxx Rxy_at Attended Unattended EE
            decsession
        end
        
        % Leave one out accuracy on self
        nofself = size(EE_dec,2);
        train = (1:nofself);
CA=[];CUA=[];
        for HH = 1:nofself % repeat for all trials - leave one out
            trials2 = train(1:end-1);

            RXX = squeeze(mean(Rxx_dec(trials2,:,:),1)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(Rxy_dec(trials2,:),1));

            avg_dec = RXX \ RXY_ATT'; % LS solution.
            recenv = EE_dec{train(end)} * avg_dec;

            CA(train(end)) = corr2(recenv,Attended_dec{train(end)}); % Attended
            CUA(train(end)) = corr2(recenv,Unattended_dec{train(end)}); % Unattended
            train = [train(end) train(1:end-1)];
        end

        difference_self = CA - CUA;
        if isempty(others)
            nofother = 0;
            difference_other = [];
            RXX = squeeze(mean(Rxx_dec,1)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(Rxy_dec,1));
            avg_dec = RXX \ RXY_ATT'; % LS solution.
            decoders{nofsessions,comb} = avg_dec;
        else
            for othersession = others;
                if startleft
                    if othersession < 5
                        sess = othersession;
                        ear = 'left';
                    elseif othersession >4
                        ear = 'right';
                        sess = othersession -4;
                    end
                else
                    if othersession < 5
                        sess = othersession;
                        ear = 'right';
                    elseif othersession >4
                        ear = 'left';
                        sess = othersession -4;
                    end
                end
                load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'],'Attended', 'Unattended', 'EE')
                Attended_other = [Attended_other, Attended{1:39}];
                Unattended_other = [Unattended_other, Unattended{1:39}];
                EE_other = [EE_other, EE{1:39}];
            end

            % Decoder & accuracy on others
            nofother = size(EE_other,2);
            RXX = squeeze(mean(Rxx_dec,1)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(Rxy_dec,1));
            avg_dec = RXX \ RXY_ATT'; % LS solution.
            decoders{nofsessions,comb} = avg_dec;
                    CA=[];CUA=[];

            for i = 1: nofother
                recenv = EE_other{i} * avg_dec;
                CA(i) = corr2(recenv,Attended_other{i}); % Attended
                CUA(i) = corr2(recenv,Unattended_other{i}); % Unattended
            end
            difference_other = CA - CUA;
        end
        noftrials = nofself + nofother;
        % Total accuracy
        difference_total(nofsessions, comb, :) = [difference_self, difference_other]; %(sum(CA-CUA> 0)/(noftrials))*100
        Total_Accuracy(nofsessions, comb) = (sum(squeeze(difference_total(nofsessions,comb,:)) > 0)/noftrials)*100;
    end
end

% Find the best decoder
[Max_Accuracy, I] = max(Total_Accuracy(:))
[nofsessions, comb] = ind2sub(size(Total_Accuracy),I);
combinations = nchoosek(sessions, nofsessions);
best_comb = combinations(comb, :)

% save
save(['All_decoders_Subject' num2str(S) '.mat'], 'difference_total', 'Total_Accuracy', 'Max_Accuracy', 'best_comb', 'decoders', 'CUA', 'CA')
avg_dec = decoders{nofsessions, comb};
difference = squeeze(difference_total(nofsessions, comb, :));
save(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec', 'difference')

%% 5b) Build the decoder
% trl = 0;
%  RXX = []; RXY_ATT = []; Rxx_all = []; Rxy_all = []; Attended_all = {}; Unattended_all = {}; EE_all = {};
% for sess = sessions
%     for e = ears
%         if e ==1
%             ear = 'left';
%         else
%             ear = 'right';
%         end
%         cd(CD)
%         load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'])
%         Rxx_all = cat(1, Rxx_all, Rxx);
%         Rxy_all = cat(1, Rxy_all, Rxy_at);
%         Attended_all = [Attended_all, Attended];
%         Unattended_all = [Unattended_all, Unattended];
%         EE_all = [EE_all, EE]
%     end
% end
% 
% RXX = squeeze(mean(Rxx_all,1)); % plain averaging of cov matrices
% RXY_ATT = squeeze(mean(Rxy_all,1));
% 
% avg_dec = RXX \ RXY_ATT';
% 
% cd(CD)
% save(['Decoder_Subject' num2str(S) '.mat'], 'avg_dec')
% save(['All4SubjAccTrainingSessions__Subject' num2str(S) '.mat'], 'Rxx_all', 'Rxy_all', 'Attended_all', 'Unattended_all', 'EE_all')

%% 6) Determine the subject accuracies for training sessions
% cd(CD)
% load(['All4SubjAccTrainingSessions__Subject' num2str(S) '.mat'])
% 
% noftrials = size(Attended_all,2);
% train = (1:noftrials);
% 
% for HH = 1:noftrials % repeat for all trials - leave one out
%     trials2 = train(1:end-1);
%         
%     RXX = squeeze(mean(Rxx_all(trials2,:,:),1)); % plain averaging of cov matrices
%     RXY_ATT = squeeze(mean(Rxy_all(trials2,:),1));
%         
%     avg_dec = RXX \ RXY_ATT'; % LS solution.
%     recenv = EE_all{train(end)} * avg_dec;
%         
%     CA(HH) = corr2(recenv,Attended_all{train(end)}); % Attended
%     CUA(HH) = corr2(recenv,Unattended_all{train(end)}); % Unattended
%     train = [train(end) train(1:end-1)];
% end
% 
% difference = CA - CUA;
% Accuracy = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials
% 
% cd(CD)
% save(['TrainingSessions_Accuracy_Subject' num2str(S) '.mat' ], 'Accuracy', 'difference')

%% 6) Determine the thresholds based on the difference statistics
S = 12

%%
load(['Best_decoder_Subject' num2str(S) '.mat'], 'difference')

% Thresholds based on percentiles
good = difference(difference > 0);
bad = difference(difference < 0);

th1 = prctile(good, 50) 
th2 = 0; % threshold for correct or incorrect decoding
th3 = prctile(bad,50)

save(['Thresholds_Subject' num2str(S) '.mat'], 'th1', 'th2', 'th3')

%% PART 2 of experiment
clear all 
close all 

%% 7) Explain the feedback to the subject & show the different colors
kleuren = [0.2 0.7 0.2; 0.6 1 0.6; 1 0.6 0.6; 0.7 0.2 0.2];
bkrcolor = [0.2 0.2 0.2];


hfig = figure('name','Performance','numbertitle','off','unit','normalized','outerposition',[0 0 1 1],'Color',bkrcolor,'menubar','none');
rectangle('Position', [0 0 1 1], 'FaceColor', bkrcolor, 'EdgeColor', bkrcolor)
ax = gca;
set(ax,'XTick',[]),set(ax,'YTick',[])
set(ax,'XColor',bkrcolor) , set(ax,'YColor',bkrcolor)
axis square
kleur = [0.4 0.4 0.4];
rectangle('Position', [0.4 0.4 0.2 0.2], 'Curvature', [1 1], 'FaceColor', kleur);

for i = 1:4
    pause
    kleur = kleuren(i,:);
    rectangle('Position', [0.4 0.4 0.2 0.2], 'Curvature', [1 1], 'FaceColor', kleur);
end
%% 8) Settings
clear all 
close all 
%%

S = 12;
sess = 4; % Which session is the subject listening to? (1, 2, 3 or 4)
ear = 'left'; % Which ear is the subject attending to? ('left' or 'right')
FB = false
SH = false

CD = ['C:\Users\Rob XPS\Desktop\THESIS\Experiments\Subject' num2str(S)];; %Fill in saving directory
cd(CD)
load(['Settings'], 'NSR', 'LP', 'HP', 'lags', 'bkrcolor', 'CD')

%% 9) Run the experiment
% instantiate the library
load(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec')

load(['Thresholds_Subject' num2str(S) '.mat'], 'th1', 'th2', 'th3')

if strcmp(ear, 'left')
    left = true; 
elseif strcmp(ear, 'right')
    left = false;
end


% instantiate the library
START = false;
GO = false;
stopp = false;

trial = 1;

% Loading and making the audio envelopes
data22 = [];
data = [];
[data22,fs] = audioread(['C:\Users\Rob XPS\Desktop\THESIS\Thesis Student 2015\Data\part' num2str(sess) '_track1_dry.wav']);
[data33,fs] = audioread(['C:\Users\Rob XPS\Desktop\THESIS\Thesis Student 2015\Data\part' num2str(sess) '_track2_dry.wav']);
data(:,1) = data22(:,1); % Left audio dry
data(:,2) = data33(:,1); % Right audio dry
if sess < 3
    AudioL = data(:,1); 
    AudioR = data(:,2);
elseif sess > 2
    AudioR = data(:,1); 
    AudioL = data(:,2);
end
envelopeL = abs(AudioL).^0.6; % Powerlaw Envelope
envelopeR = abs(AudioR).^0.6; % Powerlaw 0.6
envelopeL =  resample(double(envelopeL),500,fs);
envelopeR =  resample(double(envelopeR),500,fs);

%figure
x = [0 0 1 1 0];
y = [0 1 1 0 0];
hfig = figure('name','Performance','numbertitle','off','unit','normalized','outerposition',[0 0 1 1],'Color',bkrcolor,'menubar','none');
rectangle('Position', [0 0 1 1], 'FaceColor', bkrcolor, 'EdgeColor', bkrcolor)
ax = gca;
set(ax,'XTick',[]),set(ax,'YTick',[])
set(ax,'XColor',bkrcolor) , set(ax,'YColor',bkrcolor)
axis square
kleur = [0.4 0.4 0.4];
rectangle('Position', [0.4 0.4 0.2 0.2], 'Curvature', [1 1], 'FaceColor', kleur);

% loading the lsl library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');

result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','Stimulations'); 
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

tmr = timer('ExecutionMode', 'FixedRate', ...
    'Period', 10, ...
    'TimerFcn',...
'[chunk,stamps] = inlet.pull_chunk();data2{trial} = chunk;if trial > 1;if max(data2{trial}(25,:)) > 20;START = true;win1 = 1;[pks,locs] = findpeaks(double(data2{trial}(25,:)));if stopp;stop(tmr);clc;close;GO = false;end;end;if GO;clear New New2 AL AR envelopeL2 envelopeR2;for ch = 1:24;New(ch,:) = double(data2{trial}(ch,:));[New(ch,:)] = Rbp(LP, HP, 10, 500, New(ch,:));New2(:,ch) = resample(double(New(ch,:)),NSR,500);end;envelopeLq = envelopeL(win1:win1+size(data2{trial},2)-1);envelopeRq = envelopeR(win1:win1+size(data2{trial},2)-1);envelopeL2(1,:) =  double(envelopeLq);envelopeR2(1,:) =  double(envelopeRq);[AL(:,1)] = Rbp([], HP,10, 500, envelopeL2(1,:));[AR(:,1)] = Rbp([], HP,10, 500, envelopeR2(1,:));AL =  resample(double(AL),NSR,500);AR =  resample(double(AR),NSR,500);win1 = win1+size(data2{trial},2);WIN{trial} = win1;x=New2(:,:);start = floor(lags(1)/1e3*NSR);fin = ceil(lags(2)/1e3*NSR);nofsamples = size(x,1);[start,fin] = deal(-fin,-start);[X] = aad_LagGenerator(x,start:fin);recenv = X(:,:) * avg_dec;EE{trial} = X;Rxx(trial,:,:)= (X''*X)/nofsamples;if left;Rxy_at(trial,:) = (X''*AL)/nofsamples;Attended{trial} = AL;Unattended{trial} = AR;CUA = corr2(recenv,AR(:,:));CA = corr2(recenv,AL(:,:));else;Rxy_at(trial,:) = (X''*AR)/nofsamples;Attended{trial} = AR;Unattended{trial} = AL;CUA = corr2(recenv,AL(:,:));CA = corr2(recenv,AR(:,:));end;difference(trial) = CA-CUA;if FB; if SH; else; if difference(trial) > th1;kleur = [0.2 0.7 0.2];elseif difference(trial) > th2;kleur = [0.6 1 0.6];elseif difference(trial) > th3;kleur = [1 0.6 0.6];else;kleur = [0.7 0.2 0.2];end;rectangle(''Position'', [0.4 0.4 0.2 0.2], ''Curvature'', [1 1], ''FaceColor'', kleur);end;end;end;if START;GO = true;if GO == 1 && stopp == 0;begin = trial;end;win1 = round(win1+size(data2{trial},2)-locs);START = false;stopp = true;WIN{trial} = win1;end;end;trial = trial +1;')

start(tmr);

%% 
stop(tmr)
%%
save(['Feedback_fout_Raw_data_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'data2', 'Rxx', 'Rxy_at', 'Attended', 'Unattended', 'EE', 'WIN')

%% 10) Accuracy & saving
difference(1:begin) = [];

% Some processing before saving
% Delete the empty cells. 
EE(cellfun('isempty',EE)) = [];
Attended(cellfun('isempty', Attended)) = [];
Unattended(cellfun('isempty', Unattended)) = [];

Rxx(1:begin,:,:) = [];
Rxy_at(1:begin, :) = [];

if size(EE,2) > 39
    EE(40:end) = [];
    Attended(40:end) = [];
    Unattended(40:end) = [];
    Rxx(40:end,:,:) = [];
    Rxy_at(40:end,:) = [];
    difference(40:end) = [];
end

Session_Accuracy = (sum(difference>0)/length(difference))*100
% Rxx(cellfun('isempty',Rxx)) = [];
% Rxy_at(cellfun('isempty',Rxy_at)) = [];
% 
% % Make Rxx and Rxy_at matrices
% noftrials = size(Rxy_at, 2);
% for i = 1:noftrials
%     Rxx_new = cat(3, Rxx_new, Rxx{i});
%     Rxy_new = cat(2, Rxy_new, Rxy_at{i});
% end
% delete Rxx Rxy_at
% 
% Rxx = Rxx_new; Rxy_at = Rxy_new; 
% 
% delete Rxx_new Rxy_new

% Save
cd(CD)
save(['Feedback_Subject' num2str(S) '_Session' num2str(sess) '_' ear '.mat'], 'data2', 'Rxx', 'Rxy_at', 'Attended', 'Unattended', 'EE','WIN', 'difference', 'Session_Accuracy', 'CUA', 'CA')
