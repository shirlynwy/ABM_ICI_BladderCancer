clearvars;

% addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 3; %patients per parameter vector
cohort_pars.min_parfor_num = 4e5;
cohort_pars.last_dose_is_no_dose = true;

%%
M = allBaseParameters();
%%
M.save_pars.dt = Inf; % how frequent I save output

M.flags.fgfr3_affects_cytotoxicity = [false];
M.flags.fgfr3_affects_immune_recruit = [false];

M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

M.setup.prop_ha0 = [0.5;1];
M.setup.prop_mut0 = [0];

M.pars.max_dt = 0.25 / 24; % number of days per step

M.pars.prolif_rate = 2;
M.pars.tum_prolif_up = 0.5;

M.setup.censor_date = 10;
M.setup.N0 = 100;

M.fgfr3.n_doses = 500;
M.fgfr3.dose_val = 1e3; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)
M.fgfr3.start_day = Inf;

M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
M.checkpoint.dose_val = 100;
M.checkpoint.start_day = [Inf];

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;

M.immune_stimulatory_factor_pars.reach = 5;
M.immune_stimulatory_factor_pars.length_scale = 20;

M.immune_pars.immune_recruit_rate = 0.025;
M.immune_pars.apop_rate = 0.2;
M.immune_pars.aicd_rate = 0.7;
M.immune_pars.move_rate_microns = 2 * (24*60);
M.immune_pars.steps_per_move = 4;


%%
simCohort(M,cohort_pars);

load gong.mat
sound(y,1*8192)

% if cohort.nsamps_per_condition>=cohort.min_parfor_num
%     F(1:cohort.nsamps_per_condition) = parallel.FevalFuture;
%     ppool = gcp;
%     cohort.num_workers = ppool.NumWorkers;
%     for i = 1:cohort.nsamps_per_condition
%         F(i) = parfeval(ppool,@simPatient,1,M);
%     end
% else
%     cohort.num_workers = 1;
% end
% cohort.mu_n = 0;
% cohort.start = tic;
% cohort.batch_start = tic;
% 
% for si = cohort.nsamps_per_condition:-1:1
%     if cohort.nsamps_per_condition>=cohort.min_parfor_num
%         [idx,out_temp] = fetchNext(F);
%     else
%         idx = si;
%         out_temp = simPatient(M);
%     end
%     cohort = updateCohortTimer(cohort,cohort.nsamps_per_condition-si+1);
%     cohort.tracked(idx) = out_temp.tracked;
%     cohort.ids(idx) = out_temp.save_pars.sim_identifier;
% end
% 
% 
% save(sprintf("data/cohort_%d",tic),'-struct',"cohort")
