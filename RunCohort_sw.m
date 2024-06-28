clearvars;

% addpath("~/Documents/MATLAB/myfunctions/")
is_cpactive = 0;

%% cohort structure
cohort_pars.nsamps_per_condition = 6; %patients per parameter vector
cohort_pars.min_parfor_num = 4e5;
cohort_pars.last_dose_is_no_dose = true;

%%
M = allBaseParameters();

M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;
M.plot_pars.use_carveout = false;

%% set up

% M.flags.fgfr3_affects_cytotoxicity = [false;true];
% M.flags.fgfr3_affects_immune_recruit = [false;true];

M.setup.prop_ha0 = [0.1;0.4;0.7;0.95]; % initial proportion HA
M.setup.N0 = 1000;
M.setup.NI0 = 0;
M.setup.censor_date = 50;
% M.setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation
M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

%% basic parameters
% M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.prolif_rate = 2.36; %proliferation rate of tumor cells 1 per day
M.pars.occmax =      20;
M.pars.max_dt = 2 / 24; % number of days per step

% M.pars.tum_prolif_up = 0.5;

% 
% M.fgfr3.n_doses = 500;
% M.fgfr3.dose_val = 1e3; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)
% M.fgfr3.start_day = [1;Inf];

%% immune parameters
M.immune_stimulatory_factor_pars.length_scale = 60;


M.immune_pars.steps_per_move = 4;
M.immune_pars.aicd_rate = 0;
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day

%immune kill rates
M.immune_pars.slow_kill_rate = 1 / (2/24); % rate of slow kill in days^-1 (should be once per 2 hours)
M.immune_pars.fast_kill_rate = 1 / (.5/24); % rate of fast kill in days^-1 (should be once per half hour)
%%% fast killing probability for high and low antigen cells
% M.immune_pars.p1 = 0.95; % high antigen cells
% M.immune_pars.p2 = 0.05; % low antigen cells
%%% immune cell death due to interaction with tumor cells
immune_pars.death_by_fast = 0.01; % number of immune cells per day
immune_pars.death_by_slow = 0.01; % number of immune cells per day

%recruitment like ODE
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day
M.immune_pars.immune_recruit_rate2 = [60] ; % immune cells recruited per day
M.immune_pars.isf_prolif_hill_coefficient = 2;
M.immune_pars.deactivation_hill_coefficient_pd1 = 1;
%set to so immune cells are not deactivated in conjugation (b/c not in ODE)
M.immune_pars.deactivation_rate_pd1 = 0;

M.immune_pars.isf_prolif_ec50 = 10;
M.immune_pars.isf_prolif_max = [0.15; 0.3]; % maximum antigen stimulated proliferation rate

%% checkpoint therapy parameters
M.checkpoint.start_day = Inf;
M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
M.checkpoint.dose_val = 100;

%%%% checkpoint blocked
if ~is_cpactive
    % M.checkpoint.kf_pd1_pdl1 = 0; %can't set this to 0 bc can't /0
    M.checkpoint.kr_pd1_pdl1 = 0;
    M.checkpoint.pd1_on_immune_nmols = 0;
    M.checkpoint.pdl1_on_tumor_nmols = 0;
end

%%
simCohort(M,cohort_pars);

% load gong.mat
% sound(y,1*8192)

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
