EEG_test = EEG.data;
A = data2{3}(1:24,1);
for i = 1:size(EEG_test,2)
    EEG_min(:,i) = EEG_test(:,i) - A;
    if EEG_min(:,i) == zeros(24,1)
        EEG_min(:,i) = ones(24,1);
        location = i;
        break
    end        
end

