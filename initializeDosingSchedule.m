function drug = initializeDosingSchedule(drug,start_day,censor_date)

drug.dose_index = 1; % indexes which dose is given next

drug.times = drug.start_day + drug.days_between*(0:drug.n_doses-1)';

weekend_dose_ind = find(mod(start_day + drug.times,7) >= 5,1); % earliest therapy scheduled for a weekend

while ~isempty(weekend_dose_ind)
    change = mod(-start_day-drug.times(weekend_dose_ind,1),7); % 1 or 2 
    drug.times(weekend_dose_ind:end) = drug.times(weekend_dose_ind:end) + change;
    
    weekend_dose_ind = find(mod(start_day + drug.times,7) >= 5,1); % earliest therapy scheduled for a weekend
end

drug.dose_vals = drug.dose_val * ones(drug.n_doses,1);

%% remove doses after censor date
delete_ind = drug.times >= censor_date;
drug.times(delete_ind,:) = [];
drug.dose_vals(delete_ind,:) = [];