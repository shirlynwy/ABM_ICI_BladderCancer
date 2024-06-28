clearvars;

% addpath("~/Documents/MATLAB/myfunctions/")
% Use semicolon ";" for parameter batches, not "," 
%% cohort structure
cohort_pars.nsamps_per_condition = 6 ; %patients per parameter vector
cohort_pars.min_parfor_num = 1e4;
cohort_pars.last_dose_is_no_dose = true;

M = allBaseParameters();

%% slow and fast killing parameters
M.immune_pars.p1 = [0.92]; % fast killing probability for high antigen cells
M.immune_pars.p2 = 0.33; % fast killing probability for low antigen cells

% immune_pars.slow_kill_rate = 1 / (2/24); % rate of slow kill in days^-1 (should be once per 2 hours)
M.immune_pars.fast_kill_rate = [1 / (1/24); 1 / (0.5/24); 1 / (0.3/24)]; % rate of fast kill in days^-1 (should be once per half hour)

%% recruitment like ODE
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day
M.immune_pars.immune_recruit_rate2 = [  50; 100; 200]; % immune cells recrufited per day

%%
M.pars.prolif_rate = 3.34; %proliferation rate of tumor cells 1 per day
M.pars.occmax = 19;

%%
% M.plot_pars.plot_fig = true;
% M.plot_pars.plot_fig = false;
% M.plot_pars.plot_location = true;
% M.plot_pars.use_carveout = true;

M.setup.N0 = 20;
M.setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation 
M.setup.prop_ha0 = 1; % initial proportion HA
M.setup.censor_date = 20;
M.pars.max_dt = 2 / 24; % number of days per step
M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.setup.grid_size_microns_x = 1600;
M.setup.grid_size_microns_y = 1600;
M.setup.grid_size_microns_z = 1600;
%%
M.flags.fgfr3_affects_cytotoxicity = false;
M.flags.fgfr3_affects_immune_recruit = false;

M.fgfr3.start_day = Inf;
M.fgfr3.n_doses = 0;
M.fgfr3.dose_val = 0; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)

%%
M.immune_stimulatory_factor_pars.length_scale = 60;
M.immune_pars.steps_per_move = 4;
M.immune_pars.aicd_rate = 0;

%% checkpoint therapy parameters
M.checkpoint.start_day = Inf;
% M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
% M.checkpoint.dose_val = 100;

%set to 0 for checkpoint blockade
% M.immune_pars.deactivation_rate_pd1 = 0;

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
