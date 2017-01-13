%% Build the sham feedback
clear all 
Subjects = 1:4 % Which subject should we take into account?
startleft = [1 0 1 0] % Fill in for each subject
%%
for S = Subjects % subjectloop
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Thresholds_Subject' num2str(S) '.mat'])
    if startleft(S) == 1
        ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
    else
        ear = [0 0 1 1; 1 1 0 0];
    end
    for session = 1:4
        % load feedback results
        if ear(2, session) == 1
            EAR = 'left'; 
        else
            EAR = 'right';
        end
        load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Feedback_Subject' num2str(S) '_Session' num2str(session) '_' num2str(EAR) '.mat'], 'difference')
        for i = 1:39
            if difference(i) > th1;
                color(S,session, ear(2,session) +1, i) = 4;
            elseif difference(i) > th2;
                color(S,session, ear(2,session) +1, i) = 3;
            elseif difference(i) > th3;
                color(S,session, ear(2,session) +1, i) = 2;
            else;
                color(S,session, ear(2,session) +1, i) = 1;
            end;
        end
    end
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Best_decoder_Subject' num2str(S) '.mat'], 'difference')
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\All_decoders_Subject' num2str(S) '.mat'], 'best_comb')
    %self
    others = [1 2 7 8];
    A = 0;
    for k = 1:length(best_comb)
        others(others == best_comb(k)) = [];
    end
    order = [best_comb, others]
    for o = order
        if o <5
            sess = o
        elseif o > 4
            sess = o - 4
        end
        for i = 1:39
            if difference(i+A) > th1;
                color(S,sess, ear(1,sess) +1, i) = 4;
            elseif difference(i+A) > th2;
                color(S,sess, ear(1,sess) +1, i) = 3;
            elseif difference(i+A) > th3;
                color(S,sess, ear(1,sess) +1, i) = 2;
            else;
                color(S,sess, ear(1,sess) +1, i) = 1;
            end;
        end
        A = A + 39
    end
end

%%
colour = squeeze(cat(2, color(:,:,1,:), color(:,:,2,:)));
kleur = [];
for i = 1:8
kleur = cat(2, kleur, squeeze(colour(:,i,:)));
end
kleur = squeeze(kleur)'
