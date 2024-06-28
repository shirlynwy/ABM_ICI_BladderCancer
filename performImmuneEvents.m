function M = performImmuneEvents(M,in)

for ord_ind=1:length(in.order)
    j = in.order(ord_ind);
    switch M.immunes(j,M.I.event)
        case 1 % immune proliferates
            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.immunes(j,M.I.ind),M.immunes(j,M.I.subs),M); % neighbor indices
            if (nnz(M.L(n_ind)) + (26-length(n_ind)))<=M.immune_pars.occmax % check how many M.pars.neighbors are occupied
                % ind = randsample(length(n_ind),1,true,(M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors)); % weight by whether sites are empty and by the reciprocal of the distance
                weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                ind = find(sum(weights)*rand() < cumsum(weights),1);
                rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                M.immunes(end+1,:) = 0;
                M.immunes(end,M.I.immune_target) = NaN;
                M.immunes(end,M.I.type) = M.immunes(j,M.I.type); % inherits parent engagement status...
                M.immunes(end,M.I.subs) = M.immunes(j,M.I.subs)+rel_loc;  % store new array-specific locations
                M.immunes(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                M.L(n_ind(ind)) = M.val.imm; % set value at lattice site based on original cell (this should work instead of above commented out line)
                dt_cell = M.dt_Imm-M.immunes(j,M.I.proliferation_timer); % how long the cell actually had to proliferate
%                 M.immunes([j,end],M.I.proliferation_timer) = M.pars.min_prolif_wait - (dt_cell + log(1-rand()*(1-exp(-in.p_matrix_imm(j,1)*dt_cell)))/in.p_matrix_imm(j,1)); % this requires tracking more stuff and it's more complicated...will just assume the event happens uniformly across this small interval
                M.immunes([j,end],M.I.proliferation_timer) = M.pars.min_prolif_wait - rand()*dt_cell; % must wait M.pars.min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [M.immunes(j,M.I.proliferation_timer),M.dt] so some progress towards next prolif has happened
                M.immunes([j,end],M.I.immune_seek_time) = 0; % reset their seek clocks because they just proliferated, they're feeling great!

                M.immunes(end,M.I.region_fgfr3) = M.fgfr3.regions(M.immunes(end,M.I.ind));
                M.immunes(end,M.I.region_checkpoint) = M.checkpoint.regions(M.immunes(end,M.I.ind));
                
                M = checkpointChangePlaceImmune(M,size(M.immunes,1));

                M.tracked.imm_prolif(M.i) = M.tracked.imm_prolif(M.i) + 1;

            else
                M.immunes(j,M.I.proliferation_timer) = 0; % if not enough empty space, then allow this cell to try proliferating again
                M.tracked.imm_contact_inhibition(M.i) = M.tracked.imm_contact_inhibition(M.i)+1;
            end

        case 2 %immune cells death (taken out of ABM)
            
            if M.immunes(j,M.I.type)<=0 % then the cell is not targeting a tumor cell so there is nothing else to do here, deactivated or unengaged
                continue;
            end
            % this deals with slowing the rate of the tumor cell this
            % immune cell was attacking
            tum_ind = find(M.immunes(j,M.I.immune_target)==M.tumors(:,M.I.ind),1);
            
            if M.immunes(j,M.I.type) == 1 % was slow-killing
                M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.slow_kill_rate;
            elseif M.immunes(j,M.I.type) == 2 % was fast-killing
                M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.fast_kill_rate;
            else
                disp("unengaged or deactivated cell case(bug)")
                keyboard;
            end
            %check if no more CTL is attacking this target tumor cell
            if M.tumors(tum_ind,M.I.tumor_clearance(2)) < .5*M.immune_pars.slow_kill_rate % in case rounding errors make this hard to detect
                if rand()<M.tumors(tum_ind,M.I.tumor_clearance(1))
                    M.tumors(tum_ind,M.I.event) = 2;
                    M.L(M.tumors(tum_ind,M.I.ind)) = M.val.tum_apop;
                    M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1) = M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1)+1;
                else
                    M.tumors(tum_ind,M.I.tumor_clearance) = 0; %reset damage to 0
                end
            end

        case 3 %Im cells migrates

            % pick direction as combination of random direction and
            % gradient direction (towards ISF) with weights chosen based on
            % (norm of) ISF gradient fed into Hill function

            grad = M.immune_stimulatory_factor.gradient(:,M.immunes(j,M.I.ind));
            if all(grad==0)
                gradient_bias = 0;
            else
                temp = sqrt(sum(grad.^2));
                grad = grad/temp; % normalize the gradient for combining below
                gradient_bias = 1/(1+(M.immune_pars.isf_gradient_ec50/temp)^M.immune_pars.isf_gradient_hill_coefficient);
            end

            random_direction = uniformOnUnitSphere();
            direction = random_direction * (1-gradient_bias) + grad * gradient_bias; % no need to normalize this vector since that just scales all the weights below

            for step_ind = 1:M.immune_pars.steps_per_move
                [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.immunes(j,M.I.ind),M.immunes(j,M.I.subs),M); % neighbor indices
                if (nnz(M.L(n_ind)) + (26-length(n_ind)))<=M.immune_pars.occmax_move % check how many M.pars.neighbors are occupied
                    weights = max(0,(M.L(n_ind)==0).*...
                        (M.pars.neighbors(non_border_neighbors,:)*direction).*M.pars.neighbor_weights(non_border_neighbors).^2);
                    if any(weights>0)
                        % ind = randsample(length(n_ind),1,true,weights);
                        ind = find(sum(weights)*rand() < cumsum(weights),1);
                        rel_loc = M.pars.neighbors(non_border_neighbors(ind),:);
                        if step_ind==1
                            M.L(M.immunes(j,M.I.ind)) = 0;
                            M = checkpointChangeRemoveImmune(M,j);
                        end
                        M.immunes(j,M.I.subs) = M.immunes(j,M.I.subs) + rel_loc;
                        M.immunes(j,M.I.ind) = n_ind(ind);

                    else % no direction along gradient to move
                        M.tracked.imm_movement_halts(M.i) = M.tracked.imm_movement_halts(M.i) + 1;
                        break;
                    end
                else % no open spots to move to
                    M.tracked.imm_movement_halts(M.i) = M.tracked.imm_movement_halts(M.i) + 1;
                    break;
                end
            end
            if step_ind > 1 % otherwise it did not move
                
                M.L(M.immunes(j,M.I.ind)) = M.val.imm; % this should be very expensive
                M.immunes(j,M.I.region_fgfr3) = M.fgfr3.regions(M.immunes(j,M.I.ind)); % this should be very expensive
                M.immunes(j,M.I.region_checkpoint) = M.checkpoint.regions(M.immunes(j,M.I.ind)); % this should be very expensive

                M = checkpointChangePlaceImmune(M,j);
            end

        case 4 % immune cell conjugates with tumor cell

            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.immunes(j,M.I.ind),M.immunes(j,M.I.subs),M); % neighbor indices
            n_ind_with_tumor_log = any(M.L(n_ind)==M.val.tum,2);
            if any(n_ind_with_tumor_log)
                n_ind_with_tumor_inds = find(n_ind_with_tumor_log);
%                 temp_ind = randsample(length(n_ind_with_tumor_inds),1,true,M.pars.neighbor_weights(non_border_neighbors(n_ind_with_tumor_inds)));
                weights = M.pars.neighbor_weights(non_border_neighbors(n_ind_with_tumor_inds));
                temp_ind = find(sum(weights)*rand() < cumsum(weights),1);
                target_ind = n_ind(n_ind_with_tumor_inds(temp_ind));
                tum_ind = find(M.tumors(:,M.I.ind)==target_ind);
                
                if M.flags.fgfr3_affects_cytotoxicity && rand()>1/(1+M.phiD(tum_ind)/M.fgfr3.imm_evasion_ec50)
                    continue; % if phiD is affecting immune clearance, then include chance of evading immune system; bigger phiD==>greater evasiveness
                end

                M.immunes(j,M.I.immune_seek_time) = 0; % reset the seek time clock
                    
                M.immunes(j,M.I.immune_target) = M.tumors(tum_ind,M.I.ind); % the immune cell tracks which tumor cell it is attacking
                                
                temp_prob = rand;
                if M.tumors(tum_ind,M.I.type)==0 %low antigen cell
                    if temp_prob < M.immune_pars.p2
                        M.immunes(j,M.I.type) = 2; %1for slow kill, 2 for fast kill
                    else
                        M.immunes(j,M.I.type) = 1;
                    end
                else % high antigen cell
                    if temp_prob < M.immune_pars.p1
                        M.immunes(j,M.I.type) = 2; %1for slow kill, 2 for fast kill
                    else
                        M.immunes(j,M.I.type) = 1;
                    end
                end
              
%                 if M.tumors(tum_ind,M.I.type)==0 % immune cell is gonna kill it slowly
                if M.immunes(j,M.I.type) == 1 % immune cell is gonna kill it slowly

                    M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) + M.immune_pars.slow_kill_rate;
                    
%                     DO NOT DELETE THESE LINES! THEY EXPLAIN THE UNCOMMENTED LINE BELOW.
%                     time_to_start_attack = - log(1-rand().*(1-exp(-M.immune_pars.conjugation_rate.*M.dt_Imm)))./M.immune_pars.conjugation_rate;
%                     time_spent_attacking = M.dt_Imm - time_to_start_attack;
%                     attack_from_this_cell_during_this_step = time_spent_attacking * M.immune_pars.slow_kill_rate;
%                     attack_added_outside_this_function = M.dt_Imm * M.immune_pars.slow_kill_rate;
%                     amount_to_add_now = attack_from_this_cell_during_this_step - attack_added_outside_this_function;
%                     M.tumors(tum_ind,M.I.tumor_clearance(1)) = M.tumors(tum_ind,M.I.tumor_clearance(1)) + amount_to_add_now;  % assume that the immune cell engages at random time during this step, this pre-emptively removes the time NOT spent engaged
                    M.tumors(tum_ind,M.I.tumor_clearance(1)) = M.tumors(tum_ind,M.I.tumor_clearance(1)) + log(1-rand().*(1-exp(-M.immune_pars.conjugation_rate.*M.dt_Imm))) * M.immune_pars.slow_kill_rate/M.immune_pars.conjugation_rate; % assume that the immune cell engages at random time during this step, this pre-emptively removes the time NOT spent engaged                    
%                     M.immunes(j,M.I.type) = 1; % record that this cell is not using perforin
                elseif M.immunes(j,M.I.type) == 2 % immune cell is gonna kill it fast
                    M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) + M.immune_pars.fast_kill_rate;

%                     DO NOT DELETE THESE LINES! THEY EXPLAIN THE UNCOMMENTED LINE BELOW.
%                     time_to_start_attack = - log(1-rand().*(1-exp(-M.immune_pars.conjugation_rate.*M.dt_Imm)))./M.immune_pars.conjugation_rate;
%                     time_spent_attacking = M.dt_Imm - time_to_start_attack;
%                     attack_from_this_cell_during_this_step = time_spent_attacking * M.immune_pars.fast_kill_rate;
%                     attack_added_outside_this_function = M.dt_Imm * M.immune_pars.fast_kill_rate;
%                     amount_to_add_now = attack_from_this_cell_during_this_step - attack_added_outside_this_function;
%                     M.tumors(tum_ind,M.I.tumor_clearance(1)) = M.tumors(tum_ind,M.I.tumor_clearance(1)) + amount_to_add_now;  % assume that the immune cell engages at random time during this step, this pre-emptively removes the time NOT spent engaged
                    M.tumors(tum_ind,M.I.tumor_clearance(1)) = M.tumors(tum_ind,M.I.tumor_clearance(1)) + log(1-rand().*(1-exp(-M.immune_pars.conjugation_rate.*M.dt_Imm))) * M.immune_pars.fast_kill_rate/M.immune_pars.conjugation_rate; % assume that the immune cell engages at random time during this step, this pre-emptively removes the time NOT spent engaged                    %add assert command for debugging
%                     M.immunes(j,M.I.type) = 2; % record that this cell is using perforin
                else
                    error("Immune cell kill bug")
%                     keyboard;
                end
               
                M.tracked.imm_attack_hits(M.i) = M.tracked.imm_attack_hits(M.i) + 1;

            else % attack did not find a tumor cell
                M.tracked.imm_attack_misses(M.i) = M.tracked.imm_attack_misses(M.i) + 1;
            end % end of if

        case 5 % deactivation due to pd1-pdl1 checkpoint, exhaustion (still stick around)
            
            if M.immunes(j,M.I.type)>0
                % this deals with slowing the rate of the tumor cell this
                % immune cell was attacking
                tum_ind = find(M.immunes(j,M.I.immune_target)==M.tumors(:,M.I.ind),1);
                
                if M.immunes(j,M.I.type) == 1 % was slow-killing
                    M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.slow_kill_rate;
                else % was fast-killing
                    M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.fast_kill_rate;
                end
                
                if M.tumors(tum_ind,M.I.tumor_clearance(2)) < .5*M.immune_pars.slow_kill_rate % in case rounding errors make this hard to detect
                    if rand()<M.tumors(tum_ind,M.I.tumor_clearance(1))
                        M.tumors(tum_ind,M.I.event) = 2;
                        M.L(M.tumors(tum_ind,M.I.ind)) = M.val.tum_apop;
                        M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1) = M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1)+1;
                    else
                        M.tumors(tum_ind,M.I.tumor_clearance) = 0;
                    end
                end
            end
            %this part different from case 2 immune cell de
            temp_log = M.immunes(j,M.I.type)==[0,1,2];
            M.tracked.deactivations(M.i,temp_log) = M.tracked.deactivations(M.i,temp_log) + 1;
            M.immunes(j,M.I.type) = -1;

        case 6 % AICD

            M.immunes(j,M.I.event) = 2; % mark the immune cell to be removed by apoptosis
            M.tracked.imm_aicd(M.i) = M.tracked.imm_aicd(M.i)+1;

            if M.immunes(j,M.I.type)<=0 % then the cell is not targeting a tumor cell so there is nothing else to do here
                continue;
            else
                error("AICD should really only be for immune cells that are not engaged with a tumor cell")
            end
            % this deals with slowing the rate of the tumor cell this
            % immune cell was attacking
            tum_ind = find(M.immunes(j,M.I.immune_target)==M.tumors(:,M.I.ind),1);
            
            if M.immunes(j,M.I.type) == 1 % was slow-killing
                M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.slow_kill_rate;
            else % was fast-killing
                M.tumors(tum_ind,M.I.tumor_clearance(2)) = M.tumors(tum_ind,M.I.tumor_clearance(2)) - M.immune_pars.fast_kill_rate;
            end
            
            if M.tumors(tum_ind,M.I.tumor_clearance(2)) < .5*M.immune_pars.slow_kill_rate % in case rounding errors make this hard to detect
                if rand()<M.tumors(tum_ind,M.I.tumor_clearance(1))
                    M.tumors(tum_ind,M.I.event) = 2;
                    M.L(M.tumors(tum_ind,M.I.ind)) = M.val.tum_apop;
                    M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1) = M.tracked.imm_cleared(M.i,M.tumors(tum_ind,M.I.tumor_mut)+1,M.tumors(tum_ind,M.I.type)+1)+1;
                else
                    M.tumors(tum_ind,M.I.tumor_clearance) = 0;
                end
            end

        otherwise
            error("WHAT IS THIS?!")
            
    end % end of switch
    
end % end of j for