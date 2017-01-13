%% Analyse feedback

%% Matrix plot per subject
subjects = 1:12
adjmat_all = zeros(4,4)
adjmat_all_FB = zeros(4,4)
adjmat_all_NFB = zeros(4,4)
total_all = zeros(1,4)
total_all_FB = zeros(1,4)
total_all_NFB = zeros(1,4)
for S = 1:12
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Best_decoder_Subject' num2str(S) '.mat'], 'difference')
    %load(['Thresholds_Subject' num2str(S) '.mat'])

    % Thresholds based on percentiles
    good = difference(difference > 0);
    bad = difference(difference < 0);

    th1 = prctile(good, 50) 
    th2 = 0; % threshold for correct or incorrect decoding
    th3 = prctile(bad,50)
%     
    clear difference
    for session = 1:4
        if S == 3 && session == 1
            load(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference')
        elseif S == 8 && session == 2
            load(['FB_S8_ses2_R_corrected_diff&Acc.mat'], 'difference')
        else
            try
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'difference');
            catch
                load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'difference')
            end
        end
        if length(difference) > 39
            difference(40:end) = []
        end
        total{S, session} = zeros(1,4)
        for trial = 1:length(difference)
            if difference(trial) > th1;
                kleur{S,session}(trial) = 4; %groen
                total{S, session}(4) = total{S, session}(4) + 1
            elseif difference(trial) > th2;
                kleur{S,session}(trial) = 3; % licht groen
                total{S, session}(3) = total{S, session}(3) + 1
            elseif difference(trial) > th3; % licht rood
                kleur{S,session}(trial) = 2;
                total{S, session}(2) = total{S, session}(2) + 1
            else % donker rood
                kleur{S,session}(trial) = 1;
                total{S, session}(1) = total{S, session}(1) + 1
            end;
        end
        adjmat_ses{S,session} = zeros(4,4);
        for i = 1:length(difference)-1
            k = kleur{S,session}(i);
            l = kleur{S,session}(i+1);
            adjmat_ses{S,session}(k,l) = adjmat_ses{S,session}(k,l) +1
        end 

    end
    adjmat_subj{S} = adjmat_ses{S,1} + adjmat_ses{S,2} + adjmat_ses{S,3} + adjmat_ses{S,4} % sum sessions
    total_subj{S} = total{S,1} + total{S,2}+ total{S,3}+ total{S,4}
    adjmat_all = adjmat_all + adjmat_subj{S} % sum subjects
    total_all = total_all + total_subj{S}
    if S == 1 || S == 2 || S == 5  || S == 6 || S == 9 || S == 10
        adjmat_all_FB = adjmat_all_FB + adjmat_subj{S}
        total_all_FB = total_all_FB + total_subj{S}
    else
        adjmat_all_NFB = adjmat_all_NFB + adjmat_subj{S}
        total_all_NFB = total_all_NFB + total_subj{S}
    end
    % weighting
    for k = 1:4
        totrij = sum(adjmat_subj{S}(k,:))
        if totrij ~= 0
            adjmat_subj_perc{S}(k,:) = adjmat_subj{S}(k, :)/totrij*100;
        end
    end
    perc_subj = total_subj{S}/sum(total_subj{S})*100;
    for l = 1:4
        if perc_subj(l) ~= 0
            adjmat_subj_w{S}(:,l) = adjmat_subj_perc{S}(:,l)/perc_subj(l)
            adjmat_subj_diff{S}(:,l) = (adjmat_subj_perc{S}(:,l)-perc_subj(l))/perc_subj(l)
        end
    end
end

% weighting
for k = 1:4
    totrij = sum(adjmat_all(k,:))
    if totrij ~= 0
        adjmat_all_perc(k,:) = adjmat_all(k, :)/totrij*100;
    end
    totrij = sum(adjmat_all_FB(k,:))
    if totrij ~= 0
        adjmat_all_FB_perc(k,:) = adjmat_all_FB(k, :)/totrij*100;
    end
    totrij = sum(adjmat_all_NFB(k,:))
    if totrij ~= 0
        adjmat_all_NFB_perc(k,:) = adjmat_all_NFB(k, :)/totrij*100;
    end
end

perc_all = total_all/sum(total_all)*100;
perc_all_FB = total_all_FB/sum(total_all_NFB)*100;
perc_all_NFB = total_all_NFB/sum(total_all_NFB)*100;

for l = 1:4
    if perc_all(l) ~= 0
        adjmat_all_w(:,l) = adjmat_all_perc(:,l)/perc_all(l)
        adjmat_all_diff(:,l) = (adjmat_all_perc(:,l)-perc_all(l))/perc_all(l)
    end
    if perc_all_FB(l) ~= 0
        adjmat_all_FB_w(:,l) = adjmat_all_FB_perc(:,l)/perc_all_FB(l)
        adjmat_all_FB_diff(:,l) = (adjmat_all_FB_perc(:,l)-perc_all_FB(l))/perc_all_FB(l)
    end
    if perc_all_NFB(l) ~= 0
        adjmat_all_NFB_w(:,l) = adjmat_all_NFB_perc(:,l)/perc_all_NFB(l)
        adjmat_all_NFB_diff(:,l) = (adjmat_all_NFB_perc(:,l)-perc_all_NFB(l))/perc_all_NFB(l)
    end
end

%% If all all_FB and all_NFB were made by averaging the the results of of the 12 subjects... 
SFB = 0; SNFB = 0;
for S = 1:12
    adjmat_all_subj(S,:,:) = adjmat_subj_w{S}
    if S == 1 || S == 2 || S == 5  || S == 6 || S == 9 || S == 10
        SFB = SFB + 1;
        adjmat_all_subj_FB(SFB,:,:) = adjmat_subj_w{S}
    else
        SNFB = SNFB + 1;
        adjmat_all_subj_NFB(SNFB,:,:) = adjmat_subj_w{S}
    end
end

% Not subjects 1234
for S = 5:12
    adjmat_all_subj(S,:,:) = adjmat_subj_w{S}
    if S == 1 || S == 2 || S == 5  || S == 6 || S == 9 || S == 10
        SFB = SFB + 1;
        adjmat_all_subj_FB(SFB,:,:) = adjmat_subj_w{S}
    else
        SNFB = SNFB + 1;
        adjmat_all_subj_NFB(SNFB,:,:) = adjmat_subj_w{S}
    end
end

figure
subplot(1,3,1)
image(squeeze(mean(adjmat_all_subj,1)), 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('all')
caxis([0.5 2])

subplot(1,3,2)
image(squeeze(mean(adjmat_all_subj_FB,1)), 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('FB')
caxis([0.5 2])

subplot(1,3,3)
image(squeeze(mean(adjmat_all_subj_NFB,1)), 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('NFB')
caxis([0.5 2])

suptitle('Weighted Adjacency Matrix over all subjects')

% Significance
for l = 1:4
    for k = 1:4
        [h p] = ttest2(squeeze(adjmat_all_subj_FB(:,k,l)), squeeze(adjmat_all_subj_NFB(:,k,l)))
        P(k,l) = p
    end
end

%% For difference...If all all_FB and all_NFB were made by averaging the the results of of the 12 subjects... 
SFB = 0; SNFB = 0;
for S = 1:4
    adjmat_subj_diff{S}(2,:) = 0
end
for S = 1:12
    adjmat_all_subj(S,:,:) = adjmat_subj_diff{S}
    if S == 1 || S == 2 || S == 5  || S == 6 || S == 9 || S == 10
        SFB = SFB + 1;
        adjmat_all_subj_FB(SFB,:,:) = adjmat_subj_diff{S}
    else
        SNFB = SNFB + 1;
        adjmat_all_subj_NFB(SNFB,:,:) = adjmat_subj_diff{S}
    end
end

% % Per subj
% for S = 1:12
% figure
% subplot(1,2,1)
% image(squeeze(adjmat_all_subj_FB(S,:,:)), 'CDataMapping', 'scaled')
% colormap(parula)
% colorbar
% title('FB')
% 
% subplot(1,2,2)
% image(squeeze(adjmat_all_subj_NFB(S,:,:)), 'CDataMapping', 'scaled')
% colormap(parula)
% colorbar
% title('NFB')
% 
% suptitle(['Weighted Adjacency Matrix S' num2str(S)])
% end

% Mean
figure
subplot(1,2,1)
image(squeeze(mean(adjmat_all_subj_FB(3:6,:,:),1)), 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('Feedback Group')
caxis([-15 15])


subplot(1,2,2)
image(squeeze(mean(adjmat_all_subj_NFB(3:6,:,:),1)), 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('No Feecback Group')
caxis([-15 15])

suptitle('Adjacency Matrices')

%% %%%%%%%%%%%

figure
subplot(1,2,1)

matFB = squeeze(mean(adjmat_all_subj_FB(3:6,:,:),1))*100;           %# A 5-by-5 matrix of random values from 0 to 1
imagesc(matFB);            %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
caxis([-50 100])
title('Feedback Group')

textStrings = num2str(matFB(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:4);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(matFB(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'XTick',1:4,...                         %# Change the axes tick marks
        'XTickLabel',{'Very Bad','Bad','Good','Very Good'},...  %#   and tick labels
        'YTick',1:4,...
        'YTickLabel',{'Very Bad','Bad','Good','Very Good'},...
        'TickLength',[0 0]);
    
subplot(1,2,2)
    matNFB = squeeze(mean(adjmat_all_subj_NFB,1))*100;
imagesc(matNFB);
colormap(flipud(gray));
caxis([-50 100])
title('No Feecback Group')


textStrings = num2str(matNFB(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:4);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(matNFB(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'XTick',1:4,...                         %# Change the axes tick marks
        'XTickLabel',{'Very Bad','Bad','Good','Very Good'},...  %#   and tick labels
        'YTick',1:4,...
        'YTickLabel',{'Very Bad','Bad','Good','Very Good'},...
        'TickLength',[0 0]);
    
suptitle('Adjacency Matrices')

%%
for S = 5:12
        TOT(S,:) = total{S,1} + total{S,2} + total{S,3} + total{S,4}
        Distri(S,:) = TOT(S,:)/sum(TOT(S,:))
end

Didi = Distri * 100
%% Significance
for l = 1:4
    for k = 1:4
        [h p] = ttest2(squeeze(adjmat_all_subj_FB(3:6,k,l)), squeeze(adjmat_all_subj_NFB(:,k,l)))
        P(k,l) = p
    end
end
%% figures
% figure
% subplot(1,3,1)
% image(adjmat_all_w, 'CDataMapping', 'scaled')
% colormap(parula)
% colorbar
% title('all')
% caxis([0.5 2])

figure

subplot(1,2,1)
image(adjmat_all_FB_diff, 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('FB')
caxis([-15 15])

subplot(1,2,2)
image(adjmat_all_NFB_diff, 'CDataMapping', 'scaled')
colormap(parula)
colorbar
title('NFB')
caxis([-15 15])

suptitle('Weighted Adjacency Matrix over all subjects')

%% Matrix plot per subject

for S = 1:12
    figure
    image(adjmat_subj_w{S}, 'CDataMapping', 'scaled')
    colormap(parula)
    colorbar
    title(['Subject ' num2str(S)])
    caxis([0.5 2])
end


%% Super graphs => were not very clear
SM = {zeros(4,3), zeros(4,3), zeros(4,3), zeros(4,3);...
      zeros(4,3), zeros(4,3), zeros(4,3), zeros(4,3);...
      zeros(4,3), zeros(4,3), zeros(4,3), zeros(4,3);...
      zeros(4,3), zeros(4,3), zeros(4,3), zeros(4,3)}

SM_FB = {zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2)}
SM_NFB = {zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2);...
        zeros(3,2), zeros(3,2), zeros(3,2), zeros(3,2)}

for k = 1:4
    for l = 1:4
        SFB = 0; SNFB = 0;
        for S = 1:12
            SM{k,l}(S) = adjmat_subj_w{S}(k,l)
            if S == 1 || S == 2 || S == 5  || S == 6 || S == 9 || S == 10
                SFB = SFB +1;
                SM_FB{k,l}(SFB) = adjmat_subj_w{S}(k,l)
            else
                SNFB = SNFB + 1;
                SM_NFB{k,l}(SNFB) = adjmat_subj_w{S}(k,l)
            end
                
        end
    end
end



M = [SM{1,1}, SM{1,2}, SM{1,3}, SM{1,4};...
    SM{2,1}, SM{2,2}, SM{2,3}, SM{2,4}; ...
    SM{3,1}, SM{3,2}, SM{3,3}, SM{3,4}; ...
    SM{4,1}, SM{4,2}, SM{4,3}, SM{4,4}]

M_FB = [SM_FB{1,1}, SM_FB{1,2}, SM_FB{1,3}, SM_FB{1,4};...
        SM_FB{2,1}, SM_FB{2,2}, SM_FB{2,3}, SM_FB{2,4}; ...
        SM_FB{3,1}, SM_FB{3,2}, SM_FB{3,3}, SM_FB{3,4}; ...
        SM_FB{4,1}, SM_FB{4,2}, SM_FB{4,3}, SM_FB{4,4}]

M_NFB = [SM_NFB{1,1}, SM_NFB{1,2}, SM_NFB{1,3}, SM_NFB{1,4};...
        SM_NFB{2,1}, SM_NFB{2,2}, SM_NFB{2,3}, SM_NFB{2,4}; ...
        SM_NFB{3,1}, SM_NFB{3,2}, SM_NFB{3,3}, SM_NFB{3,4}; ...
        SM_NFB{4,1}, SM_NFB{4,2}, SM_NFB{4,3}, SM_NFB{4,4}]


% figures
figure
subplot(1,3,1)
image(M, 'CDataMapping', 'scaled')
colormap(hot)
colorbar
title('all')
caxis([0 5])

subplot(1,3,2)
image(M_FB, 'CDataMapping', 'scaled')
colormap(hot)
colorbar
title('FB')
caxis([0 5])

subplot(1,3,3)
image(M_NFB, 'CDataMapping', 'scaled')
colormap(hot)
colorbar
title('NFB')
caxis([0 5])

suptitle('Weighted Adjacency Matrix over all subjects')

% => too confusing to use! 


%% Color ifv trial % OPM even aangepast naar training ipv feedback!
subjects = 1:12
groen = zeros(1,312)
lgroen = zeros(1,312)
lrood = zeros(1,312)
rood = zeros(1,312)

FBgroen = zeros(1,312)
FBlgroen = zeros(1,312)
FBlrood = zeros(1,312)
FBrood = zeros(1,312)

NFBgroen = zeros(1,312)
NFBlgroen = zeros(1,312)
NFBlrood = zeros(1,312)
NFBrood = zeros(1,312)

for S = subjects
    trl = 0
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Best_decoder_Subject' num2str(S) '.mat'], 'difference')

    % Thresholds based on percentiles
    good = difference(difference > 0);
    bad = difference(difference < 0);

    th1 = prctile(good, 50) 
    th2 = 0; % threshold for correct or incorrect decoding
    th3 = prctile(bad,50)
    
    clear difference
    for session = 1:4
%         if S == 3 && session == 1
%             load(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference')
%         end
        try
            load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'difference');
        catch
            load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'difference')
        end
        if length(difference) > 39
            difference(40:end) = []
        end
        total{S, session} = zeros(1,4)
        for trial = 1:length(difference)
            trl = trl + 1
            if S == 1 || S == 2 || S == 5 || S == 6 || S == 9 || S == 10 %fb
                if difference(trial) > th1;
                    groen(trl) = groen(trl)+1
                    FBgroen(trl) = FBgroen(trl)+1
                elseif difference(trial) > th2;
                    lgroen(trl) = lgroen(trl)+1
                    FBlgroen(trl) = FBlgroen(trl)+1
                elseif difference(trial) > th3;
                    lrood(trl) = lrood(trl)+1
                    FBlrood(trl) = FBlrood(trl)+1
                else;
                    rood(trl) = rood(trl)+1
                    FBrood(trl) = FBrood(trl)+1
                end;
            else
                if difference(trial) > th1;
                    groen(trl) = groen(trl)+1
                    NFBgroen(trl) = NFBgroen(trl)+1
                elseif difference(trial) > th2;
                    lgroen(trl) = lgroen(trl)+1
                    NFBlgroen(trl) = NFBlgroen(trl)+1
                elseif difference(trial) > th3;
                    lrood(trl) = lrood(trl)+1
                    NFBlrood(trl) = NFBlrood(trl)+1
                else;
                    rood(trl) = rood(trl)+1
                    NFBrood(trl) = NFBrood(trl)+1
                end;
            end
        end
    end
end


%Cummulative
cgroen(1) = groen(1) 
clgroen(1) = lgroen(1)
clrood(1) = lrood(1)
crood(1) = rood(1)

FBcgroen(1) = FBgroen(1) 
FBclgroen(1) = FBlgroen(1)
FBclrood(1) = FBlrood(1)
FBcrood(1) = FBrood(1)

NFBcgroen(1) = NFBgroen(1) 
NFBclgroen(1) = NFBlgroen(1)
NFBclrood(1) = NFBlrood(1)
NFBcrood(1) = NFBrood(1)
for trl = 2:312
    cgroen(trl) = groen(trl) + cgroen(trl-1)
    clgroen(trl) = lgroen(trl) + clgroen(trl-1)
    clrood(trl) = lrood(trl) + clrood(trl-1)
    crood(trl) = rood(trl) + crood(trl-1)
    
    FBcgroen(trl) = FBgroen(trl) + FBcgroen(trl-1)
    FBclgroen(trl) = FBlgroen(trl) + FBclgroen(trl-1)
    FBclrood(trl) = FBlrood(trl) + FBclrood(trl-1)
    FBcrood(trl) = FBrood(trl) + FBcrood(trl-1)
    
    NFBcgroen(trl) = NFBgroen(trl) + NFBcgroen(trl-1)
    NFBclgroen(trl) = NFBlgroen(trl) + NFBclgroen(trl-1)
    NFBclrood(trl) = NFBlrood(trl) + NFBclrood(trl-1)
    NFBcrood(trl) = NFBrood(trl) + NFBcrood(trl-1)
    
end
figure
plot(1:156, cgroen(1:156), 'g')
hold on
plot(1:156, clgroen(1:156), 'y')
plot(1:156, clrood(1:156), 'm')
plot(1:156, crood(1:156), 'r')
hold off

figure
plot(1:156, FBcgroen(1:156), '-- g')
hold on
plot(1:156, FBclgroen(1:156), '-- y')
plot(1:156, FBclrood(1:156), '-- m')
plot(1:156, FBcrood(1:156), '-- r')
hold off

figure
plot(1:156, NFBcgroen(1:156), 'g')
hold on
plot(1:156, NFBclgroen(1:156), 'y')
plot(1:156, NFBclrood(1:156), 'm')
plot(1:156, NFBcrood(1:156), 'r')
hold off

%% Diference and CUA and CA analysis: put these variables of all subjects together
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

for S = 1:12
    ses = 0;
    for TF = 1:2
        if startleft(S) == 1
           ear = [1 1 0 0; 0 0 1 1];  % left = 1; right = 0
        else
           ear = [0 0 1 1; 1 1 0 0];
        end
        for session = 1:4
            ses = ses + 1;
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
            elseif TF == 2
                if S == 3 && session == 1
                    load(['FB_S3_ses1_R_corrected_diff&Acc.mat'], 'difference')
                elseif S == 8 && session == 2
                    load(['FB_S8_ses2_R_corrected_diff&Acc.mat'], 'difference')
                else
                    load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                end
            end
            if length(difference) > 39
                difference(40:end) = [];
            end
            D(S,ses,:) = difference;
        end
    end
end                          
 
%% First second of 'changes' !!!! OPM: OPNIEUW DOEN WANT S8 was ook een foutje!!! 
subjects = 1:12;
sort_trials = cell(12, 4, 4, 4);

for S = subjects
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Thresholds_Subject' num2str(S) '.mat'])
    for session = 1:4
        if S == 3 && session == 1
            load(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference')
        end
        % if S = 8 ! ...
        try
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'difference');
        catch
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'difference')
        end
        if length(difference) > 39
            difference(40:end) = []
        end
        for trial = 1:length(difference)
            if difference(trial) > th1;
                kleur{S,session}(trial) = 4;
            elseif difference(trial) > th2;
                kleur{S,session}(trial) = 3;
            elseif difference(trial) > th3;
                kleur{S,session}(trial) = 2;
            else;
                kleur{S,session}(trial) = 1;
            end;
        end
        for trial = 2:length(difference)
            k = kleur{S,session}(trial-1);
            l = kleur{S,session}(trial);
            sort_trials{S, session, k,l} = [sort_trials{S, session , k,l}, trial]
        end 
    end
end
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save('Sorted_trials', 'sort_trials')

%% Get the EEG data 100 previous samples and 200 after sections... 
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load('Sorted_trials')

for S = 1:12
    num = 0;
    for session = 1:4
        cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
        try
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'data2');
        catch
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'data2')
        end
        data2(cellfun('isempty',data2)) = []; %Remove empty cells
        if size(data2,2)>39
            data2(40:end) = [];
        end
        for k = 1:4
            for l = 1:4
                noftrials = length(sort_trials{S,session,k,l})
                for trl = 1: noftrials
                    num = num+1;
                    trial = sort_trials{S,session,k,l}(trl)+1; % +1 Because we want the 10s AFTER the change was shown to the subject!
                    if trial < 40
                        EEG_sorted_12{S, k,l, num} = [data2{trial-1} data2{trial}];
                        z = size(data2{trial-1},2)
                        for ch = 1:24
                            EEG_sorted_filtered{S, k,l,num}(ch,:) =  Rbp(1, 30, 6, 500, squeeze(EEG_sorted_12{S, k,l,num}(ch,:)));
                        end
                        EEG_sorted_cut{S,k,l}(:,:,num) = EEG_sorted_filtered{S,k,l,num}(:,z-99 : z + 200)
                        earref = (EEG_sorted_cut{S,k,l}(18,:,num)+ EEG_sorted_cut{S,k,l}(19,:,num))/2
                        for ch = 1:24
                            EEG_sorted_reref {S,k,l}(ch,:,num)= EEG_sorted_cut{S,k,l}(ch,:,num)-earref
                            first100 = mean(squeeze(EEG_sorted_reref {S,k,l}(ch,1:100,num)))
                            EEG_sorted_final{S,k,l}(ch,:,num) = squeeze(EEG_sorted_reref {S,k,l}(ch,:,num)) - first100
                        end
                    end
                end
            end
        end
    end
end

EEG_sorted_FB = EEG_sorted_final([1 2 5 6 9 10],:,:);
EEG_sorted_NFB = EEG_sorted_final([3 4 7 8 11 12],:,:);
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['EEG_sorted'], 'EEG_sorted_12', 'EEG_sorted_reref', 'EEG_sorted_final', 'EEG_sorted_filtered', 'EEG_sorted_cut', 'EEG_sorted_FB', 'EEG_sorted_NFB', '-v7.3')

%% Put all together for comparing FB and NFB
%cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
%load(['EEG_sorted'])

EEG_sorted_FB_all = cell(4,4);
EEG_sorted_NFB_all = cell(4,4);
EEG_sorted_FB_rtg = []; EEG_sorted_FB_gtr = []; % red to green and green to red
EEG_sorted_NFB_rtg = []; EEG_sorted_NFB_gtr = [];
EEG_sorted_FB_nochange = []; EEG_sorted_NFB_nochange = [];

for k = 1:4
    for l = 1:4
        for S = 1:6
            if isempty(EEG_sorted_FB{S,k,l})
            else
                EEG_sorted_FB_all{k,l} = cat(3,EEG_sorted_FB_all{k,l}, squeeze(EEG_sorted_FB{S,k,l}));
                if k < 3 && l > 2
                    EEG_sorted_FB_rtg = cat(3,EEG_sorted_FB_rtg, squeeze(EEG_sorted_FB{S,k,l}));
                elseif k > 2 && l < 3
                    EEG_sorted_FB_gtr = cat(3,EEG_sorted_FB_gtr, squeeze(EEG_sorted_FB{S,k,l}));
                elseif k == l
                    EEG_sorted_FB_nochange = cat(3,EEG_sorted_FB_nochange, squeeze(EEG_sorted_FB{S,k,l}));
                end
            end
            if isempty(EEG_sorted_NFB{S,k,l})
            else
                EEG_sorted_NFB_all{k,l} = cat(3,EEG_sorted_NFB_all{k,l}, squeeze(EEG_sorted_NFB{S,k,l}));
                if k < 3 && l > 2
                    EEG_sorted_NFB_rtg = cat(3,EEG_sorted_NFB_rtg, squeeze(EEG_sorted_NFB{S,k,l}));
                elseif k > 2 && l < 3
                    EEG_sorted_NFB_gtr = cat(3,EEG_sorted_NFB_gtr, EEG_sorted_NFB{S,k,l});
                elseif k == l
                    EEG_sorted_NFB_nochange = cat(3,EEG_sorted_NFB_nochange, squeeze(EEG_sorted_NFB{S,k,l}));
                end
            end
        end
    end
end

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['More_sorted_EEG.mat'], 'EEG_sorted_FB_all', 'EEG_sorted_FB_rtg', 'EEG_sorted_FB_gtr','EEG_sorted_FB_nochange','EEG_sorted_NFB_all', 'EEG_sorted_NFB_rtg', 'EEG_sorted_NFB_gtr', 'EEG_sorted_NFB_nochange', 'EEG_sorted_FB', 'EEG_sorted_NFB')

%% Averaging & plots

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load(['More_sorted_EEG.mat'], 'EEG_sorted_FB_all', 'EEG_sorted_FB_rtg', 'EEG_sorted_FB_gtr','EEG_sorted_FB_nochange','EEG_sorted_NFB_all', 'EEG_sorted_NFB_rtg', 'EEG_sorted_NFB_gtr', 'EEG_sorted_NFB_nochange', 'EEG_sorted', 'EEG_sorted_FB', 'EEG_sorted_NFB')

FB_rtg = squeeze(mean(EEG_sorted_FB_rtg,3));
st_FB_rtg = std(EEG_sorted_FB_rtg,[],3);
FB_gtr = squeeze(mean(EEG_sorted_FB_gtr,3));
st_FB_gtr = std(EEG_sorted_FB_gtr,[],3);
FB_nochange = squeeze(mean(EEG_sorted_FB_nochange,3));
st_FB_nochange = std(EEG_sorted_FB_nochange,[],3);

% NFB_rtg = squeeze(mean(EEG_sorted_NFB_rtg,3));
% st_NFB_rtg = std(EEG_sorted_NFB_rtg,[],3);
% NFB_gtr = squeeze(mean(EEG_sorted_NFB_gtr,3));
% st_NFB_gtr = std(EEG_sorted_NFB_gtr,[],3);
% NFB_nochange = squeeze(mean(EEG_sorted_NFB_nochange,3));
% st_NFB_nochange = std(EEG_sorted_NFB_nochange,[],3);
% 
for ch = 1:24
    figure
    plot(-99:200, FB_rtg(ch,:), 'g')
    hold on
    plot(-99:200, FB_gtr(ch,:), 'r')
    plot(-99:200, FB_nochange(ch,:), 'k')
    title(['Channel ' num2str(ch)])
    hold off
end

%%
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load(['More_sorted_EEG.mat'], 'EEG_sorted_FB_all', 'EEG_sorted_FB_rtg', 'EEG_sorted_FB_gtr','EEG_sorted_FB_nochange','EEG_sorted_NFB_all', 'EEG_sorted_NFB_rtg', 'EEG_sorted_NFB_gtr', 'EEG_sorted_NFB_nochange', 'EEG_sorted', 'EEG_sorted_FB', 'EEG_sorted_NFB')

for ch = 1:24
    for i = 1:size(EEG_sorted_FB_gtr,3)
        fEEG_sorted_FB_gtr(ch,:,i) = Rbp(1, 30, 6, 500, squeeze(EEG_sorted_FB_gtr(ch,:,i)));
    end
    for i = 1:size(EEG_sorted_FB_rtg,3)
        fEEG_sorted_FB_rtg(ch,:,i) = Rbp(1, 30, 6, 500, squeeze(EEG_sorted_FB_rtg(ch,:,i)));
    end
    for i = 1:size(EEG_sorted_FB_nochange,3)
        fEEG_sorted_FB_nochange(ch,:,i) = Rbp(1, 30, 6, 500, squeeze(EEG_sorted_FB_nochange(ch,:,i)));
    end
end

fFB_rtg = squeeze(mean(fEEG_sorted_FB_rtg,3));
fst_FB_rtg = std(fEEG_sorted_FB_rtg,[],3);
fFB_gtr = squeeze(mean(fEEG_sorted_FB_gtr,3));
fst_FB_gtr = std(fEEG_sorted_FB_gtr,[],3);
fFB_nochange = squeeze(mean(fEEG_sorted_FB_nochange,3));
fst_FB_nochange = std(fEEG_sorted_FB_nochange,[],3);

% fNFB_rtg = squeeze(mean(EEG_sorted_NFB_rtg,3));
% fst_NFB_rtg = std(EEG_sorted_NFB_rtg,[],3);
% fNFB_gtr = squeeze(mean(EEG_sorted_NFB_gtr,3));
% fst_NFB_gtr = std(EEG_sorted_NFB_gtr,[],3);
% fNFB_nochange = squeeze(mean(EEG_sorted_NFB_nochange,3));
% fst_NFB_nochange = std(EEG_sorted_NFB_nochange,[],3);

for ch = 1:24
    figure
    subplot(1,3,1)
    plot(1:500, fFB_rtg(ch,:), 'g')
    title('red to green')
    subplot(1,3,2)
    plot(1:500, fFB_gtr(ch,:), 'r')
    title('green to red')
    subplot(1,3,3)
    plot(1:500, fFB_nochange(ch,:), 'k')
    title('no change in color')
    suptitle(['Channel ' num2str(ch)])
end

