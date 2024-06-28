function M = updateImmune(M)

if M.NI==0
    return;
end

M.immune_stimulatory_factor.gradient = zeros([M.grid.size,3]);
[M.immune_stimulatory_factor.gradient(:,:,:,2),M.immune_stimulatory_factor.gradient(:,:,:,1),M.immune_stimulatory_factor.gradient(:,:,:,3)] = gradient(M.immune_stimulatory_factor.concentration); % MATLAB's gradient assumes the "vertical" direction in the matrix is the y direction; I use the second dimension (columns, i.e. "horizontal" in the matrix) for y
M.immune_stimulatory_factor.gradient = reshape(M.immune_stimulatory_factor.gradient,[],3)'; % so we use the index of cells to grab the gradient

%% disengage immune cells that had been targeting a tumor cell that just apoptosed naturally
dead_targets = M.tumors(:,M.I.event)==2 & M.tumors(:,M.I.tumor_clearance(2))>0; % M.tumors that were being targeted by CTLs but died naturally
if any(dead_targets)
    M.tumors(dead_targets,M.I.tumor_clearance) = 0; % these cells no longer have progress towards immune clearance, nor are making progress towards that
    attacking_immunes = find(M.immunes(:,M.I.type)>0); % all immune cells that were attacking
    immunes_done_lysing = attacking_immunes(ismember(M.immunes(attacking_immunes,M.I.immune_target),M.tumors(dead_targets,M.I.ind))); % search across attacking immune cells for those that were targeting tumor cells that are now dead
    M.immunes(immunes_done_lysing,[M.I.type,M.I.immune_target]) = repmat([0,NaN],length(immunes_done_lysing),1); % reset these immune cells as to type (0=not engaged and not deactivated) and location inds of the cell they are targeting
end

for ImIter=1:M.Nsteps_Imm
    alive_log = M.immunes(:,M.I.event)~=2; % logical of which immune cells are alive
    alive_ind = find(alive_log); % indices of alive immune cells
    unengaged_and_activated_subalive_log = M.immunes(alive_ind,M.I.type)==0; % of the alive immune cells, those that are both unengaged and activated (not deactivated yet)
    activated_subalive_log = M.immunes(alive_ind,M.I.type)~=-1; % of the alive immune cells, those that are activated
    engaged_subalive_log = activated_subalive_log & ~unengaged_and_activated_subalive_log; % immune cells currently conjugated with a tumor cell and so susceptible to pd1-induced exhaustion/deactivation
    n_alive = length(alive_ind);

    % I could at some point say deactivated/exhausted ctls undergo AICD;
    % this could be interesting under the assumption that ICI therapy can
    % allow ctls to recover from exhaustion

    M = updateCheckpoint(M);

    %%%%%% old by Daniel, for deactivating T cells  that are conjugated with a tumor cell
    temp_pd1 = zeros(n_alive,1);
    temp_pd1(engaged_subalive_log) = M.pd1_pdl1_equilibrium(M.immunes(engaged_subalive_log,M.I.region_checkpoint));

    temp_pd1_deactivation = (temp_pd1 / M.immune_pars.deactivation_ec50_pd1).^M.immune_pars.deactivation_hill_coefficient_pd1;
    temp_deactivation_rate = M.immune_pars.deactivation_function(temp_pd1_deactivation,M);

    %%%%%% new by Shirlyn, for suppressing antigen-stimulated CTL proliferation
    temp_alive_pd1 = M.pd1_pdl1_equilibrium(M.immunes(alive_ind,M.I.region_checkpoint));

    temp_pd1_deactivation2 = (temp_alive_pd1 / M.immune_pars.deactivation_ec50_pd1).^M.immune_pars.deactivation_hill_coefficient_pd1;
    immune_suppress_factor = 1 ./ (1 + temp_pd1_deactivation2);
    %%%%%%

    temp_isf_conc = M.immune_stimulatory_factor.concentration(M.immunes(alive_ind,M.I.ind));
    temp_isf_prolif = (temp_isf_conc / M.immune_pars.isf_prolif_ec50).^M.immune_pars.isf_prolif_hill_coefficient;
    
    %%%%%%% new by Shirlyn: for CTL death from conjugating with tumor cells
    engaged_and_fast_log = (M.immunes(alive_ind,M.I.type)==2); % of the alive immune cells, those that are both unengaged and activated (not deactivated yet)
    engaged_and_slow_log = (M.immunes(alive_ind,M.I.type)==1);
    ctl_death_from_fast = zeros(n_alive,1);
    ctl_death_from_fast(engaged_and_fast_log) = M.immune_pars.death_by_fast;
    ctl_death_from_slow = zeros(n_alive,1);
    ctl_death_from_slow(engaged_and_slow_log) = M.immune_pars.death_by_slow;
    ctl_death_from_tumor = ctl_death_from_fast + ctl_death_from_slow;

    if length(ctl_death_from_tumor) ~= n_alive
        keyboard;
    end

%     if ~isempty(find(engaged_and_fast_log,1)) ||  ~isempty(find(engaged_and_slow_log,1))
%         keyboard;
%     end

   %line below: immune cell proliferation changed from Daniel's code
   rate_matrix_imm = [unengaged_and_activated_subalive_log .* (M.immune_pars.prolif_rate  + immune_suppress_factor .* M.immune_pars.isf_prolif_max .* temp_isf_prolif./(1+temp_isf_prolif) ),... % only proliferate if no timer within [0,M.dt_Imm)
        M.immune_pars.apop_rate * ones(n_alive,1) + ctl_death_from_tumor,... % immune apoptosis
        unengaged_and_activated_subalive_log .* M.immune_pars.move_rate/M.immune_pars.steps_per_move,... % random movement
        unengaged_and_activated_subalive_log * M.immune_pars.conjugation_rate,... % conjugation attempt probability
        temp_deactivation_rate,... % deactivation of immune cells
        (unengaged_and_activated_subalive_log & M.immunes(alive_ind,M.I.immune_seek_time) > M.immune_pars.time_to_seek) * M.immune_pars.aicd_rate]; % rate of AICD due to not finding tumor cells

    M.immunes(unengaged_and_activated_subalive_log,M.I.immune_seek_time) = M.immunes(unengaged_and_activated_subalive_log,M.I.immune_seek_time) + M.dt_Imm; % all activated and unengaged ctls increase their seek time (those that find a tumor right now will have it set to 0 in performImmuneEvents)

    in.p_matrix_imm = 1-exp(-rate_matrix_imm.*[max(0,M.dt_Imm-M.immunes(alive_ind,M.I.proliferation_timer)),M.dt_Imm*ones(n_alive,size(rate_matrix_imm,2)-1)]);

    assert(all(sum(in.p_matrix_imm,2)<=1) && (isempty(in.p_matrix_imm) || min(in.p_matrix_imm,[],'all')>=0)) % make sure the probabilities do not add up to more than 1

    [M.immunes(alive_ind,M.I.event),~] = find(diff(rand(n_alive,1)<cumsum([zeros(n_alive,1),in.p_matrix_imm,ones(n_alive,1)],2),1,2)');

    active_ind = find(alive_log & M.immunes(:,M.I.event)~=(size(in.p_matrix_imm,2)+1));
    perm_order = randperm(length(active_ind));
    in.order = active_ind(perm_order);

    M = performImmuneEvents(M,in);

    %% update proliferation timers for all non-proliferating immune cells
    non_prolif_ind = M.immunes(:,M.I.event)>1;
    M.immunes(non_prolif_ind,M.I.proliferation_timer) = max(0,M.immunes(non_prolif_ind,M.I.proliferation_timer)-M.dt_Imm); % update immune cell proliferation probabilities for all but those that just proliferated (1) and those that just arrived or were born (0)

    %% increment engaged timers
    M.tumors(:,M.I.tumor_clearance(1)) = M.tumors(:,M.I.tumor_clearance(1))+M.dt_Imm*M.tumors(:,M.I.tumor_clearance(2)); % column 6 is the rate (per day) and so multiply by the time step
    lysing_tumors_log = M.tumors(:,M.I.tumor_clearance(1))>=1;
    if any(lysing_tumors_log)
        M.tumors(lysing_tumors_log,M.I.event) = 2; % these tumor cells are dying, so set their prob clearance, clearance rate to 0; set to apoptosis
        M.tumors(lysing_tumors_log,M.I.tumor_clearance) = 0; % these tumor cells are dying, so set their prob clearance, clearance rate to 0; set to apoptosis
        
        M.L(M.tumors(lysing_tumors_log,M.I.ind)) = M.val.tum_apop;

        % (time step, mut, antigen)
        M.tracked.imm_cleared(M.i,:,:) = M.tracked.imm_cleared(M.i,:,:) + sum(M.tumors(lysing_tumors_log,M.I.tumor_mut)==[0,1] & M.tumors(lysing_tumors_log,M.I.type)==reshape([0,1],1,1,2),1);

        %% change engage status for immune cells that are now done lysing
        attacking_immunes = find(M.immunes(:,M.I.type)>0);
        immunes_done_lysing = attacking_immunes(ismember(M.immunes(attacking_immunes,M.I.immune_target),M.tumors(lysing_tumors_log,M.I.ind)));
        M.immunes(immunes_done_lysing,M.I.type) = 0;
        M.immunes(immunes_done_lysing,M.I.immune_target) = NaN;
    end
    
end %end of ImIter

apoptosis_log = M.immunes(:,M.I.event)==2;

if any(apoptosis_log)
    M.L(M.immunes(apoptosis_log,M.I.ind)) = 0;

    M = checkpointChangeRemoveImmune(M,find(apoptosis_log));

end

M.immunes(apoptosis_log,:)=[];

%% if plotting with prob_matrix, put that here
if M.plot_pars.plot_fig && mod(M.i,M.plot_pars.plot_every)==M.plot_pars.plot_offset
    M = plotFunction_EndUpdateImmune(M,in,mean(temp_pd1_deactivation(temp_pd1_deactivation~=0))); % storing order in this function is the ONLY reason I need to output M
end

