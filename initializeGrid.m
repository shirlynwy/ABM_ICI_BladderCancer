function M = initializeGrid(M)

%% setting up grid variables
M.grid.size = 1+ceil(M.setup.grid_size_microns/M.pars.cell_width); % increase the grid to size to the next multiple of cell width
M.grid.xx = 1:M.grid.size(1);
M.grid.yy = 1:M.grid.size(2);
M.grid.zz = 1:M.grid.size(3);
M.grid.center = 0.5*(1+M.grid.size); % center of the grid

M.rel_pos_ind = M.pars.neighbors*cumprod([1,M.grid.size(1:2)])'; % for Moore neighborhood
M.rel_pos_ind_VN = M.rel_pos_ind(sum(abs(M.pars.neighbors),2)==1); % for von Neumann neighborhood
M.V_tot = prod(M.grid.size); % volume of TME as volume of lattice (in cells)