function pars = basePars()

pars.max_dt = 1/4 / 24; % number of days per tumor step
pars.max_dt_Imm = 1/8/24; %number of days per immun step


pars.cell_width = 20; % in micrometers; cell is about 20micrometers in diameter
pars.min_prolif_wait = 9/24; % number of days all cells must wait at minimum between proliferations
pars.max_tumor_size = 60000; % if tumor exceeds this, stop simulation
pars.min_tumor_size = 2; %if tumor is <=2 cells, stop simulation
pars.max_tumor_boundary = 36; %if number of tumor cells on the boundary exceeds this, stop simulation

%% neighbor parameters
pars.occmax = 20; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate

%% tumor parameters

pars.prolif_rate = 2.6 ; %proliferation rate of tumor cells 1 per day
pars.apop_rate = .05; % max death rate of tumor cells (per day)


