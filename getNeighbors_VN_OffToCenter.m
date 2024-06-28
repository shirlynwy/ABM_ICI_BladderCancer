function neighbor_inds = getNeighbors_VN_OffToCenter(M,ind,subs)

% neighbors that are off the grid get replaced with the center ind

neighbor_inds = ind'+M.rel_pos_ind_VN;

cells_on_left = subs(:,1)==1;
neighbor_inds(M.pars.left_neighbors_VN_ind,cells_on_left) = ind(cells_on_left);
cells_on_right = subs(:,1)==M.grid.size(1);
neighbor_inds(M.pars.right_neighbors_VN_ind,cells_on_right) = ind(cells_on_right);
cells_on_front = subs(:,2)==1;
neighbor_inds(M.pars.front_neighbors_VN_ind,cells_on_front) = ind(cells_on_front);
cells_on_back = subs(:,2)==M.grid.size(2);
neighbor_inds(M.pars.back_neighbors_VN_ind,cells_on_back) = ind(cells_on_back);
cells_on_bottom = subs(:,3)==1;
neighbor_inds(M.pars.bottom_neighbors_VN_ind,cells_on_bottom) = ind(cells_on_bottom);
cells_on_top = subs(:,3)==M.grid.size(3);
neighbor_inds(M.pars.top_neighbors_VN_ind,cells_on_top) = ind(cells_on_top);




% out = ind+M.rel_pos_ind_VN;
% 
% border_neighbors = false(26,1);
% 
% switch subs(1)
%     case 1 % then on left boundary
%         border_neighbors(M.pars.left_neighbors_VN_ind) = true;
%     case M.grid.size(1) % then on right boundary
%         border_neighbors(M.pars.right_neighbors_VN_ind) = true;
% end
% 
% switch subs(2)
%     case 1 % then on front boundary
%         border_neighbors(M.pars.front_neighbors_VN_ind) = true;
%     case M.grid.size(2) % then on back boundary
%         border_neighbors(M.pars.back_neighbors_VN_ind) = true;
% end
% 
% switch subs(3)
%     case 1 % then on bottom boundary
%         border_neighbors(M.pars.bottom_neighbors_VN_ind) = true;
%     case M.grid.size(3) % then on top boundary
%         border_neighbors(M.pars.top_neighbors_VN_ind) = true;
% end
% 
% out(border_neighbors) = ind; % point to the center spot for any index off the grid
