function runBatchSample(I)
rng shuffle  % Daniel's: use a different sequence of random numbers for each sim
% runs sample i for the cohort defined here

%% cohort structure
cohort_pars.nsamps_per_condition = 10;
cohort_pars.last_dose_is_no_dose = false;


%%
M = allBaseParameters();
%%
M.save_pars.dt = 0.25;

M.flags.fgfr3_affects_cytotoxicity = [false;true];
M.flags.fgfr3_affects_immune_recruit = [false;true];

M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

M.setup.prop_ha0 = [0;0.5;1];
M.setup.prop_mut0 = [0;0.5;1];

M.pars.max_dt = 0.25 / 24; % number of days per step

M.pars.prolif_rate = 2;
M.pars.tum_prolif_up = 0.5;

M.setup.censor_date = 0.5;
M.setup.N0 = 100;

M.fgfr3.n_doses = 500;
M.fgfr3.dose_val = 1e3; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;

M.immune_stimulatory_factor_pars.reach = 5;
M.immune_stimulatory_factor_pars.length_scale = 20;

M.immune_pars.immune_recruit_rate = 0.025;
M.immune_pars.apop_rate = 0.2;
M.immune_pars.aicd_rate = 0.7;
M.immune_pars.move_rate_microns = 2 * (24*60);
M.immune_pars.steps_per_move = 4;

M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
M.checkpoint.dose_val = 100;

%%
simBatchSample(M,cohort_pars,I);