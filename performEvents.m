function M = performEvents(M,in)

for ord_ind=1:length(in.active_ind)

    j = in.active_ind(ord_ind);
    mut_ind = M.tumors(j,M.I.tumor_mut)+1;
    ant_ind = M.tumors(j,M.I.type)+1;

    switch in.events(ord_ind)
        case 1 % proliferation
            %%
%             [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumors(j,:),M); % neighbor indices
            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumors(j,M.I.ind),M.tumors(j,M.I.subs),M); % neighbor indices
            if (nnz(M.L(n_ind)) + (26-length(n_ind)))<=M.pars.occmax % check how many M.pars.neighbors are occupied
                weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                %                 ind = randsample(length(n_ind),1,true,weights); % weight by whether sites are empty and by the reciprocal of the distance
                ind = find(sum(weights)*rand() < cumsum(weights),1);
                rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                M.tumors(end+1,:) = 0;
                M.tumors(end,M.I.type) = M.tumors(j,M.I.type); % inherits parent antigen expression...
                M.tumors(end,M.I.tumor_mut) = M.tumors(j,M.I.tumor_mut); % ...and mutation burden;
                M.tumors(end,M.I.subs) = M.tumors(j,M.I.subs)+rel_loc;  % store new array-specific locations
                M.tumors(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                M.L(n_ind(ind)) = M.L(M.tumors(j,M.I.ind)); % set value at lattice site based on original cell (this should work instead of above commented out line)

                M.tumors([j,end],M.I.proliferation_timer) = M.pars.min_prolif_wait + in.time_to_event(ord_ind); % must wait M.pars.min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [M.tumors(j,M.I.proliferation_timer),M.dt] so some progress towards next prolif has happened (will subtract off M.dt with all cells in simForward)
                M.tumors([j,end],M.I.tumor_receptors(1)) = M.tumors(j,M.I.tumor_receptors)*[1;1;0;.5;1;1]; % monomer concentrations stay the same; all complexes get split evenly between the two, so to compensate half the number of monomers in them is added back to the monomer species
                M.tumors([j,end],M.I.tumor_receptors([2,4:end])) = .5*M.tumors([j,j],M.I.tumor_receptors([2,4:end])); % complexes end up split evenly between the two

                M.phiD([j,end+1]) = M.tumors([j,end],M.I.tumor_receptors(2))/M.fgfr3.RT; % set  phiD value for all tumor cells

                M.tumors(end,M.I.region_fgfr3) = M.fgfr3.regions(M.tumors(end,M.I.ind));
                M.tumors(end,M.I.tumor_receptors(3)) = M.fgfr3.concentration(M.tumors(end,M.I.region_fgfr3));
                M.tumors(end,M.I.region_checkpoint) = M.checkpoint.regions(M.tumors(end,M.I.ind));

                M.tracked.tum_prolif(M.i,mut_ind,ant_ind) = M.tracked.tum_prolif(M.i,mut_ind,ant_ind)+1;

                [xx,yy,zz,xind,yind,zind] = computeISFIndices(M,M.tumors(end,M.I.subs));
                if M.tumors(end,M.I.type)==0
                    M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) + M.immune_pars.low_antigen_isf_factor * M.immune_stimulatory_factor.stamp(xind,yind,zind);
                else
                    M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) + M.immune_stimulatory_factor.stamp(xind,yind,zind);
                end
            else
                M.tumors(j,M.I.proliferation_timer) = 0; % if not enough empty space, then allow this cell to try proliferating again
                M.tracked.tum_contact_inhibition(M.i,mut_ind,ant_ind) = M.tracked.tum_contact_inhibition(M.i,mut_ind,ant_ind)+1;
            end
            
        case 2 % spontaneous apoptosis
            %%
            M.L(M.tumors(j,M.I.ind)) = M.val.tum_apop;
            M.tracked.tum_apop(M.i,mut_ind,ant_ind) = M.tracked.tum_apop(M.i,mut_ind,ant_ind)+1;
        otherwise
            %%
            error('should not do nothing')
            
    end % end of switch
end % end of j for