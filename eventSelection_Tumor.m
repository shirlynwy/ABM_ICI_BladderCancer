function [M,out] = eventSelection_Tumor(M)

rate_matrix = computeRateMatrix(M);

prob_matrix = 1-exp(-rate_matrix.*[max(0,M.dt-M.tumors(:,M.I.proliferation_timer)),M.dt*ones(M.NT,1)]); % use M.dt as time step except for proliferation, use M.dt - remaining time to wait before next proliferation

event_log_array = rand(size(prob_matrix))<prob_matrix; % choose whether events happen

inactive_log = all(~event_log_array,2); % index of all tumor cells that do nothing
[out.active_ind,out.events] = find(event_log_array); % index of all tumor cells that do something and the out.events they will do

final_event = zeros(M.NT,1,'uint8');
final_event(inactive_log) = size(rate_matrix,2)+1; % all inactive tumor cells do nothing (which is marked by final event index + 1)

u = rand(length(out.events),1); % random number to decide when event occurs
min_wait = (out.events==1).*M.tumors(out.active_ind,M.I.proliferation_timer); % minimum time waited until proliferation event can possibly occur
rates = rate_matrix(sub2ind(size(rate_matrix),out.active_ind,out.events)); % rates at which all the selected out.events happen
dts = M.dt - min_wait; % adjust total dts for proliferation

out.time_to_event = min_wait + (- log(1-u.*(1-exp(-rates.*dts)))./rates); % given that the event happens in this M.dt update, this determines when it occurs in that update step; this is why it's not just min_wait-log(rand())./rates

[out.time_to_event,event_order] = sort(out.time_to_event);
out.active_ind = out.active_ind(event_order);
out.events = out.events(event_order);

%% remove proliferations that happen after apoptosis

apop_event_ind = find(out.events==2);

remove_ind = false(length(out.active_ind),1);
for i = 1:length(apop_event_ind)
    ai = apop_event_ind(i); % find event index corresponding to this apoptosis   
    remove_ind(ai+find(out.active_ind(ai+1:end)==out.active_ind(ai))) = true; % get ready to remove all out.events by this apoptotic cell that occur after apoptosis
end

out.active_ind(remove_ind) = []; % remove these post-apoptotic events
out.events(remove_ind) = []; % remove these post-apoptotic events
out.time_to_event(remove_ind) = []; % remove these post-apoptotic events

final_event(out.active_ind) = out.events; % this should set the event status to the latest-occuring event (really, this should mean that if both prolif and apop happen in that order, then the event will be labeled apop)
    
M.tumors(:,M.I.event) = final_event; % label the tumor based on the final event it performs

%% if plotting with prob_matrix, put that here
if M.plot_pars.plot_fig && mod(M.i,M.plot_pars.plot_every)==M.plot_pars.plot_offset
    plotFunction_EndEventSelectionTumor(M, prob_matrix);
end
