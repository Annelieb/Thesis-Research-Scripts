%% Determine significance levels at 0.05 chance level.
noftrials = [312, 152, 104, 72, 56, 48];

for TRL = 1:6
    N = noftrials(TRL); 
    siglevel(TRL) = binoinv(0.95, N, 0.5)/N;
end

siglevel
%cd('C:\Users\Annelies\Documents\Studie\Thesis\script\Significance levels')
%save('siglevel.mat', 'siglevel')