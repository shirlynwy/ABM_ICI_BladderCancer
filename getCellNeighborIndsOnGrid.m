function [out,non_border_neighbors] = getCellNeighborIndsOnGrid(cell_ind,cell_subs,M)

out = cell_ind+M.rel_pos_ind;

% if M.setup.ndims == 3
    border_neighbors = false(26,1);

    switch cell_subs(1)
        case 1 % then on left boundary
            border_neighbors(M.pars.left_neighbors) = true;
        case M.grid.size(1) % then on right boundary
            border_neighbors(M.pars.right_neighbors) = true;
    end

    switch cell_subs(2)
        case 1 % then on front boundary
            border_neighbors(M.pars.front_neighbors) = true;
        case M.grid.size(2) % then on back boundary
            border_neighbors(M.pars.back_neighbors) = true;
    end

    switch cell_subs(3)
        case 1 % then on bottom boundary
            border_neighbors(M.pars.bottom_neighbors) = true;
        case M.grid.size(3) % then on top boundary
            border_neighbors(M.pars.top_neighbors) = true;
    end

% else
%     border_neighbors = false(8,1);
% 
%     switch cell_subs(1)
%         case 1 % then on left boundary
%             border_neighbors(M.pars.left_neighbors) = true;
%         case M.grid.size(1) % then on right boundary
%             border_neighbors(M.pars.right_neighbors) = true;
%     end
% 
%     switch cell_subs(2)
%         case 1 % then on front boundary
%             border_neighbors(M.pars.front_neighbors) = true;
%         case M.grid.size(2) % then on back boundary
%             border_neighbors(M.pars.back_neighbors) = true;
%     end
% 
% end

out(border_neighbors) = [];
non_border_neighbors = find(~border_neighbors);
%% old version
% function [out,non_border_neighbors] = getCellNeighborIndsOnGrid(cell_row,M)
%
% out = cell_row(M.I.ind)+M.rel_pos_ind;
%
% border_neighbors = false(26,1);
%
% switch cell_row(M.I.subs(1))
%     case 1 % then on left boundary
%         border_neighbors(M.pars.left_neighbors) = true;
%     case M.grid.size(1) % then on right boundary
%         border_neighbors(M.pars.right_neighbors) = true;
% end
%
% switch cell_row(M.I.subs(2))
%     case 1 % then on front boundary
%         border_neighbors(M.pars.front_neighbors) = true;
%     case M.grid.size(2) % then on back boundary
%         border_neighbors(M.pars.back_neighbors) = true;
% end
%
% switch cell_row(M.I.subs(3))
%     case 1 % then on bottom boundary
%         border_neighbors(M.pars.bottom_neighbors) = true;
%     case M.grid.size(3) % then on top boundary
%         border_neighbors(M.pars.top_neighbors) = true;
% end
%
% out(border_neighbors) = [];
% non_border_neighbors = find(~border_neighbors);