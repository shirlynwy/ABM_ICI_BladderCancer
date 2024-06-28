function M = updateCheckpoint(M)

ic = [M.checkpoint.PD1;... % col 1
      M.checkpoint.aPD1;... % col 2
      M.checkpoint.PD1aPD1]; % col 3
ic = [ic(:);M.checkpoint.aPD1_circulation];

assert(all(ic>=0))

M.checkpoint.pd1_proportion = M.checkpoint.volumes.pd1./M.checkpoint.volumes.regions;

sol = ode15s(@(t,x) globalCheckpointODE(x,M),[0,M.dt_Imm],ic(:));
Y = sol.y(:,end);

assert(all(Y>=0,'all'))

M.checkpoint.aPD1_circulation = Y(end);

Y = reshape(Y(1:end-1),M.checkpoint.sz);

M.checkpoint.PD1 = Y(:,1);
M.checkpoint.aPD1 = Y(:,2);
M.checkpoint.PD1aPD1 = Y(:,3);

M.pd1_pdl1_equilibrium = 0.5*((M.checkpoint.PD1+M.checkpoint.pdl1_on_tumor) + ...
    M.checkpoint.Kd_pd1_pdl1 - ...
    sqrt((M.checkpoint.PD1+M.checkpoint.pdl1_on_tumor+M.checkpoint.Kd_pd1_pdl1).^2-4*M.checkpoint.PD1*M.checkpoint.pdl1_on_tumor));


% this block explains the below uncommented lines. DO NOT DELETE
% WITHOUT ADDING EXPLANATION
%     neighbor_types = L(immunes(onGrid,ind_ind)'+rel_pos_ind); % look at the types of cells that neighbor all immune cells
%     neighbor_is_tumor = any(neighbor_types(:)==tum_vals_and_tum_apop_val,2); % find all the neighbors that are tumor cells that will contribute PD1PDL1; include apoptotic tumor cells because they can still have PDL1  that could interact with PD1 on CTLs (and because in the code here apoptotic tumor cells do count towards the reacting volume of PD1)
%     nvals(onGrid) = sum(reshape(neighbor_is_tumor,6,[]),1);

%% get lattice indices for all neighbor sites; for those on the border, use the center site
% neighbor_inds = M.immunes(alive_ind,M.I.ind)'+M.rel_pos_ind_VN;
% cells_on_left = M.immunes(alive_ind,M.I.subs(1))==1;
% neighbor_inds(M.pars.left_neighbors_VN_ind,cells_on_left) = M.immunes(cells_on_left,M.I.ind);
% cells_on_right = M.immunes(alive_ind,M.I.subs(1))==M.grid.size(1);
% neighbor_inds(M.pars.right_neighbors_VN_ind,cells_on_right) = M.immunes(cells_on_right,M.I.ind);
% cells_on_front = M.immunes(alive_ind,M.I.subs(2))==1;
% neighbor_inds(M.pars.front_neighbors_VN_ind,cells_on_front) = M.immunes(cells_on_front,M.I.ind);
% cells_on_back = M.immunes(alive_ind,M.I.subs(2))==M.grid.size(2);
% neighbor_inds(M.pars.back_neighbors_VN_ind,cells_on_back) = M.immunes(cells_on_back,M.I.ind);
% cells_on_bottom = M.immunes(alive_ind,M.I.subs(3))==1;
% neighbor_inds(M.pars.bottom_neighbors_VN_ind,cells_on_bottom) = M.immunes(cells_on_bottom,M.I.ind);
% cells_on_top = M.immunes(alive_ind,M.I.subs(3))==M.grid.size(3);
% neighbor_inds(M.pars.top_neighbors_VN_ind,cells_on_top) = M.immunes(cells_on_top,M.I.ind);


% neighbor_inds = getNeighbors_VN_OffToCenter(M,M.immunes(alive_ind,M.I.ind),M.immunes(alive_ind,M.I.subs));
%     
% number_tumor_neighbors = sum(any(M.L(neighbor_inds)==reshape(M.val.tum,[1,1,2]),3),1)'; % see commented out block above for a step-by-step breakdown of all this
% 
% if any(number_tumor_neighbors>0)
%     disp('')
% end
% 
% M.pd1_signaling_on_alive_immunes = number_tumor_neighbors.*M.pd1_pdl1_equilibrium(M.immunes(alive_ind,M.I.region_checkpoint));
