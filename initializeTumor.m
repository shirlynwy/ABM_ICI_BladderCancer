function M = initializeTumor(M)


radius = nthroot((3/4)*(1/pi())*M.setup.N0,3); % radius of sphere with volume equal to the number of cells
all_subs = allCombos(M.grid.xx,M.grid.yy,M.grid.zz);

M.tumors = zeros(M.setup.N0,M.I.tumor_receptors(end));
d = sqrt(sum((all_subs-M.grid.center).^2,2));

M.tumors(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

ha_ind = randsample(M.setup.N0,floor(M.setup.prop_ha0*M.setup.N0),false);
M.ha_log = false(M.setup.N0,1);
M.ha_log(ha_ind) = true;
M.tumors(ha_ind,M.I.type) = 1; % this will contain cell-type information (0=low antigen tumor cell; 1=high antigen tumor cell)
mut_ind = randsample(M.setup.N0,floor(M.setup.prop_mut0*M.setup.N0),false);
M.mut_log = false(M.setup.N0,1);
M.mut_log(mut_ind) = true;
M.tumors(mut_ind,M.I.tumor_mut) = 1; % 1=mutated in fgfr3, 0=no mutation
M.tumors(:,M.I.tumor_clearance(1)) = 0; % proportion cleared by immune cell (if it reaches 1, then the tumor is marked for apoptosis; if its clearance rate is set to 0, then it is cleared with this probability)
M.tumors(:,M.I.tumor_clearance(2)) = 0; % rate at which the tumor is being cleared by the immune system.
M.tumors(:,M.I.event) = 0; % event index

% M.tumors(:,M.I.tumor_receptors(1)) = random(RT_dist,M.setup.N0,1).*(M.tumors(:,M.I.tumor_mut)==[0,1])*[non_mut_RT_factor;1]; % inactive monomers of FGFR3 (nmol/m^3)
M.tumors(~M.mut_log,M.I.tumor_receptors(1)) = M.fgfr3.non_mut_RT_factor * M.fgfr3.RT; % inactive monomers of FGFR3 (nmol/m^3)
M.tumors(M.mut_log,M.I.tumor_receptors(1)) = M.fgfr3.RT; % inactive monomers of FGFR3 (nmol/m^3)
M.tumors(:,M.I.tumor_receptors(2)) = 0; % active dimers of signaling FGFR3 (nmol/m^3)
M.tumors(:,M.I.tumor_receptors(3)) = 0; % inhibitor (nmol/m^3)
M.tumors(:,M.I.tumor_receptors(4)) = 0; % inhibitor-receptor complex (nmol/m^3)
M.tumors(:,M.I.tumor_receptors(5)) = 0; % inhibitor-receptor dimer (nmol/m^3)
M.tumors(:,M.I.tumor_receptors(6)) = 0; % inhibitor-dimer complex (nmol/m^3)

%% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
prolif_rate = M.pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
u = rand(M.setup.N0,1);
i1 = u>=1/(1+M.pars.min_prolif_wait*prolif_rate);
i2 = ~i1;
x = zeros(M.setup.N0,1);
x(i1) = u(i1)*(1+M.pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.pars.min_prolif_wait;
x(i2) = log((1+M.pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.pars.min_prolif_wait;
M.tumors(:,M.I.proliferation_timer) = max(0,x+M.pars.min_prolif_wait); % time until next possible proliferation (days)

%% figure out where tumors belong on grid
% (tumor location relative to start of lattice)/(length of each step) =
% num of lattice "right" of start; add 1 because start has index 1
M.tumors(:,M.I.ind) = sub2ind(M.grid.size,M.tumors(:,M.I.subs(1)),M.tumors(:,M.I.subs(2)),M.tumors(:,M.I.subs(3))); % linear indices
M.tumors(:,M.I.region_fgfr3) = M.fgfr3.regions(M.tumors(:,M.I.ind));

M.NT = size(M.tumors,1);