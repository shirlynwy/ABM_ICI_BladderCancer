clearvars;

% addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

M.setup.prop_ha0 = 1; % initial proportion HA
M.setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation
%%
M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.flags.fgfr3_affects_cytotoxicity = false;
M.flags.fgfr3_affects_immune_recruit = false;

M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

M.setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation 

M.pars.max_dt = 2 / 24; % number of days per step

M.pars.prolif_rate = 2; %proliferation rate of tumor cells 1 per day

M.setup.censor_date = 20;
M.setup.N0 = 20;

M.fgfr3.start_day = Inf;
M.fgfr3.n_doses = 0;
M.fgfr3.dose_val = 0; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)

M.plot_pars.plot_fig = true;
M.plot_pars.plot_location = true;
M.plot_pars.use_carveout = true;

% M.immune_stimulatory_factor_pars.length_scale = 60;
M.immune_pars.steps_per_move = 4;
M.immune_pars.aicd_rate = 0;

%tumor escapes
M.immune_pars.apop_rate = 0;
% immune_pars.immune_recruit_rate  = 0.05;

%recruitment like ODE
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day
M.immune_pars.immune_recruit_rate2 = 0; % immune cells recrufited per day

%checkpoint therapy parameters
M.checkpoint.start_day = Inf;
% M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
% M.checkpoint.dose_val = 100;

%set to 0 for checkpoint blockade
% M.immune_pars.deactivation_rate_pd1 = 0;

%fast killing probability for high and low antigen cells 
M.immune_pars.p1 = 1; % high antigen cells
M.immune_pars.p2 = 0; % low antigen cells



M = simPatient(M);
