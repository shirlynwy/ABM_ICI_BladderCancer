clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 3;
cohort_pars.min_parfor_num = 4e5;
cohort_pars.last_dose_is_no_dose = true;

%%
M = allBaseParameters();
%%

M.flags.fgfr3_affects_cytotoxicity = [false;true];
M.flags.fgfr3_affects_immune_recruit = [false;true];

M.fgfr3.start_day = [1;8];
M.checkpoint.start_day = [1;8];

M.save_pars.dt = 0.125;

M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

M.pars.max_dt = 0.25 / 24; % number of days per step

M.pars.prolif_rate = 2;
M.pars.tum_prolif_up = 0.5;

M.setup.censor_date = 50;
M.setup.N0 = 100;

M.fgfr3.n_doses = 500;
M.fgfr3.dose_val = 1e3; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;

M.immune_stimulatory_factor_pars.length_scale = 20;
M.immune_pars.steps_per_move = 4;

M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
M.checkpoint.dose_val = 100;

%%
simCohort(M,cohort_pars);
