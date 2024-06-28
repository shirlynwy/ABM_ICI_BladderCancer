function M = finishParameterSetup_Patient(M)

M.setup.grid_size_microns = [M.setup.grid_size_microns_x,M.setup.grid_size_microns_y,M.setup.grid_size_microns_z];

%% neighbor stuff
M.pars.neighbors = allCombos(-1:1,-1:1,-1:1,'matlab');
M.pars.neighbors(all(M.pars.neighbors==0,2),:) = []; % don't count self as neighbor
M.pars.neighbors_VN = M.pars.neighbors;
M.pars.neighbors_VN(sum(abs(M.pars.neighbors_VN),2)>1,:) = [];

M.pars.neighbor_weights = 1./sqrt(sum(M.pars.neighbors.^2,2));
M.pars.left_neighbors = M.pars.neighbors(:,1) == -1;
M.pars.right_neighbors = M.pars.neighbors(:,1) == 1;
M.pars.front_neighbors = M.pars.neighbors(:,2) == -1;
M.pars.back_neighbors = M.pars.neighbors(:,2) == 1;
M.pars.bottom_neighbors = M.pars.neighbors(:,3) == -1;
M.pars.top_neighbors = M.pars.neighbors(:,3) == 1;

M.pars.left_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == -1);
M.pars.right_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == 1);
M.pars.front_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == -1);
M.pars.back_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == 1);
M.pars.bottom_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == -1);
M.pars.top_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == 1);

%% index and value stuff
M.I = buildIndices();
M.val = buildLatticeVals();

%% cell/fgfr3 stuff
M.pars.V = M.pars.cell_width^3;

M.fgfr3.rectangle_compartment_size = [M.fgfr3.rectangle_compartment_size_x,M.fgfr3.rectangle_compartment_size_y,M.fgfr3.rectangle_compartment_size_z];

M.fgfr3.RT = M.fgfr3.RT_num / M.pars.V;

M.fgfr3.receptors_without_inhibitor = (M.fgfr3.kr+M.fgfr3.kp - sqrt((M.fgfr3.kr+M.fgfr3.kp)^2+8*M.fgfr3.kf*(M.fgfr3.kr+M.fgfr3.kp)*M.fgfr3.RT))/(-4*M.fgfr3.kf); % if no anti-fgfr3 therapy, then all tumor cells quickly reach this concentration of dimer
M.fgfr3.all_receptors_without_inhibitor = [M.fgfr3.receptors_without_inhibitor,(M.fgfr3.RT-M.fgfr3.receptors_without_inhibitor)/2,zeros(1,4)];   % if no anti-fgfr3 therapy, then all tumor cells quickly reach this concentration of fgfr3 stuff


%% events
M = initializeEvents(M);

%% immune stuff
M.immune_pars.move_rate = M.immune_pars.move_rate_microns / M.pars.cell_width;

%% immune stimulatory factor
r = M.immune_stimulatory_factor_pars.reach;
l = 2*r+1; % length of each dimension in stamp
M.immune_stimulatory_factor.interval = -r:r;

d = allCombos(1:l,1:l,1:l,'matlab') - r - 1;
d = sqrt(sum(d.^2,2)); % distance to center point
d = reshape(d,l,l,l);

k = M.pars.cell_width * log(2)/M.immune_stimulatory_factor_pars.length_scale; % set decay so that at length scale away from center (tumor site), ISF is at half center value
M.immune_stimulatory_factor.stamp = exp(-k*d);
M.immune_stimulatory_factor.stamp = M.immune_stimulatory_factor.stamp / M.immune_stimulatory_factor.stamp(r+2,r+1,r+1); % normalize so one space away from tumor cell has value of 1

M.immune_pars.isf_gradient_ec50 = 0.5*diff(M.immune_stimulatory_factor.stamp(r+[4,2],r+1,r+1)); % set the ec50 for gradient sensitivity to the gradient at two spaces away from tumor cell

%% checkpoint stuff
M.checkpoint.rectangle_compartment_size = [M.checkpoint.rectangle_compartment_size_x,M.checkpoint.rectangle_compartment_size_y,M.checkpoint.rectangle_compartment_size_z];

M.checkpoint.diffusion_factor_apd1 = 6*M.checkpoint.diffusivity_apd1 / (M.pars.cell_width^2);
M.checkpoint.diffusion_factor_apdl1 = 6*M.checkpoint.diffusivity_apdl1 / (M.pars.cell_width^2);

M.checkpoint.pd1_on_immune = (M.checkpoint.pd1_on_immune_nmols ...
    /M.pars.V)... % divide by total volume of reaction space to get concentration of PDL1 in nmol/um^3
    *1e15; % multiply by 1e15um^3/L to get concentration of PDL1 in nM

M.checkpoint.pdl1_on_tumor = (M.checkpoint.pdl1_on_tumor_nmols ...
    /M.pars.V)... % divide by total volume of reaction space to get concentration of PDL1 in nmol/um^3
    *1e15; % multiply by 1e15um^3/L to get concentration of PDL1 in nM

v1 = M.checkpoint.pd1_on_immune;
v2 = M.checkpoint.pdl1_on_tumor;
M.checkpoint.Kd_pd1_pdl1 = M.checkpoint.kr_pd1_pdl1/M.checkpoint.kf_pd1_pdl1;
pd1pdl1_equilibrium = 0.5*((v1+v2)+ M.checkpoint.Kd_pd1_pdl1 - sqrt((v1+v2+M.checkpoint.Kd_pd1_pdl1)^2-4*v1*v2));

if M.checkpoint.kr_pd1_pdl1~= 0
    M.immune_pars.deactivation_ec50_pd1 = M.immune_pars.deactivation_ec50_pd1_proportion * pd1pdl1_equilibrium;
else
    M.immune_pars.deactivation_ec50_pd1 = 0.01;
end


