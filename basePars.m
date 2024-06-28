function pars = basePars()

pars.max_dt = 9 / 24; % number of days per step
pars.max_dt_Imm = (1/6) / 24; % number of days per immune step


pars.cell_width = 20; % in micrometers; cell is about 20micrometers in diameter
pars.min_prolif_wait = 9/24; % number of days all cells must wait at minimum between proliferations
pars.max_tumor_size = Inf; % if tumor exceeds this, stop simulation
pars.min_tumor_size = 2; %if tumor is <=2 cells, stop simulation
pars.max_tumor_boundary = 54; %if number of tumor cells on the boundary exceeds this, stop simulation

%% neighbor parameters
pars.occmax = 20; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate

%% tumor parameters

pars.prolif_rate = 1.5/(1-1.5*9/24); % proliferation rate of tumor cells (per day); from Global/Step18/test_determine_apoptosis.m saw that proliferation rate (r4 there) should be about 1.5, but want the minimum time before proliferation to be 9 hours, so need to increase proliferation rate here to make the average still be 1.5 proliferations / day
pars.apop_rate = .05; % max death rate of tumor cells (per day)


