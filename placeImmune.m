function M = placeImmune(M,n_newI)

potential_immune_ind = find(M.blood_vessels);
potential_immune_ind = setdiff(potential_immune_ind,M.tumors(:,M.I.ind));
potential_immune_ind = setdiff(potential_immune_ind,M.immunes(:,M.I.ind));

n_inds = length(potential_immune_ind);
if n_inds == 0
    return;
end
n_newI = min(n_newI,n_inds);

if n_inds>1
    new_immune_inds = randsample(potential_immune_ind,n_newI,false);
else
    new_immune_inds = potential_immune_ind;
end
M.tracked.imm_recruit(M.i) = length(new_immune_inds);

new_immune_subs = zeros(n_newI,3);
[new_immune_subs(:,1),new_immune_subs(:,2),new_immune_subs(:,3)] = ind2sub(M.grid.size,new_immune_inds);

%% immune cells arrive
new_immunes = zeros(n_newI,size(M.immunes,2));
new_immunes(:,M.I.subs) = new_immune_subs;
new_immunes(:,M.I.type) = 0; % current immune killing type: -1=deactivated; 0=unengaged; 1=only FasL (slow); 2=perforin + FasL (fast)
new_immunes(:,M.I.event) = 0; % immune event
new_immunes(:,M.I.immune_target) = NaN; % locations of tumor cells being targeted
new_immunes(:,M.I.proliferation_timer) = 0; % time until the cell can proliferate again

new_immunes(:,M.I.ind) = new_immune_inds;

new_immunes(:,M.I.region_fgfr3) = M.fgfr3.regions(new_immunes(:,M.I.ind));
new_immunes(:,M.I.region_checkpoint) = M.checkpoint.regions(new_immunes(:,M.I.ind));

new_immunes(:,M.I.immune_seek_time) = -M.immune_pars.extra_seek_time_on_arrival; % they have not yet found a tumor cell, so don't let them become exhausted so fast

new_inds = size(M.immunes,1) + (1:n_newI);

M.L(new_immunes(:,M.I.ind)) = M.val.imm;
M.immunes = cat(1,M.immunes,new_immunes);

M = checkpointChangePlaceImmune(M,new_inds);