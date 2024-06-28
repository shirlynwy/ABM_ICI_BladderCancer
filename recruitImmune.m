function M = recruitImmune(M)

if M.flags.fgfr3_affects_immune_recruit
    rate = M.dt*M.immune_pars.immune_recruit_rate*M.NT*(M.immune_pars.min_imm_recruit_prop+(1-M.immune_pars.min_imm_recruit_prop)/(1+mean(M.phiD)/M.immune_pars.min_imm_recruit_prop_ec50));
else
    %     rate = M.dt*M.immune_pars.immune_recruit_rate*M.NT;
    if M.NI == 0
        rate = M.dt*M.immune_pars.immune_recruit_rate2;
    else
        alive_log = M.immunes(:,M.I.event)~=2; % logical of which immune cells are alive
        alive_ind = find(alive_log); % indices of alive immune cells
        temp_alive_pd1 = M.pd1_pdl1_equilibrium(M.immunes(alive_ind,M.I.region_checkpoint));
        temp_pd1_deactivation2 = (temp_alive_pd1 / M.immune_pars.deactivation_ec50_pd1).^M.immune_pars.deactivation_hill_coefficient_pd1;
        immune_suppress_factor = 1 ./ (1 + temp_pd1_deactivation2);
        mean_factor = mean(immune_suppress_factor);

        rate = M.dt*M.immune_pars.immune_recruit_rate2 * mean(immune_suppress_factor);
    end
end
n_newI = poissrnd(rate);
if n_newI>0

    M = placeImmune(M,n_newI);

    %     new_immune_neighbors_inds = immunes(1:n_newI,ind_ind)'+rel_pos_ind_VN;
    %     tumor_neighbors_log = reshape(any(reshape(L(new_immune_neighbors_inds),[],1)==tum_vals,2),6,[]); % I probably need to do some reshaping and stuff to get this check right

end
