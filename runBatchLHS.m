function runBatchLHS(I)
% rng(I) % worked on hpc

%testing below
% tempstr = string(datetime("now","Format","MMddHHmmss"));
tempstr = string(datetime("now","Format","yyMMddHHmm"));
tempnum = str2double(tempstr);
s1 = RandStream.create('mt19937ar','seed', tempnum+I);
RandStream.setGlobalStream(s1)  % use a different sequence of random numbers for each sim

is_cpactive = 0;
%% cohort structure
cohort_pars.nsamps_per_condition = 10;
cohort_pars.last_dose_is_no_dose = false;


%%
M = allBaseParameters();
%%
M.save_pars.dt = 0.5; % set to Inf to not save anything; otherwise dt is in days
M.save_pars.batch_name = strcat('lhs_cpblk_10par_', string(datetime("now","Format","yyMMddHH")));

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;
M.plot_pars.use_carveout = false;

%% set up

M.setup.prop_ha0 = 0.5; %linspace(0.1, 1, 10); % initial proportion HA
M.setup.N0 = 20;
M.setup.NI0 = 0;
M.setup.censor_date = 150;
% M.setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation
M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

%% basic parameters
M.pars.max_dt = 1/4 / 24; % number of days per tumor step
M.pars.max_dt_Imm = 1/8/24;%number of days per immune step
M.pars.prolif_rate = 2.6; %proliferation rate of tumor cells 1 per day
M.pars.occmax =  20;
M.pars.max_tumor_size = 60000; 
M.pars.max_tumor_boundary = 36; %if number of tumor cells on the boundary exceeds this, stop simulation
M.pars.min_tumor_size = 2;

% M.pars.tum_prolif_up = 0.5;


%% immune parameters
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
% M.immune_pars.death_by_fast = 0.01; % number of immune cells per day
% M.immune_pars.death_by_slow = 0.01; % number of immune cells per day

%recruitment like ODE
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day
M.immune_pars.immune_recruit_rate2 = 10; % immune cells recruited per day
% M.immune_pars.isf_prolif_hill_coefficient = 2;
% M.immune_pars.deactivation_hill_coefficient_pd1 = 1;
%set to so immune cells are not deactivated in conjugation (b/c not in ODE)
M.immune_pars.deactivation_rate_pd1 = 0;

M.immune_pars.isf_prolif_ec50 = 10;
M.immune_pars.isf_prolif_max = 0.15; % maximum antigen stimulated proliferation rate

%% checkpoint therapy parameters
M.checkpoint.start_day = Inf;
% M.checkpoint.diffusivity_apd1 = 0.1 * 60*60*24;
M.checkpoint.dose_val = 0;

%%%% checkpoint blocked
if ~is_cpactive
    % M.checkpoint.kf_pd1_pdl1 = 0; %can't set this to 0 bc can't /0
    M.checkpoint.kr_pd1_pdl1 = 0;
    M.checkpoint.pd1_on_immune_nmols = 0;
    M.checkpoint.pdl1_on_tumor_nmols = 0;
end

%%
lhs_filename = 'LHSsamples_abm_03-27-24_001.mat';
simBatchLHS(M,cohort_pars,I, lhs_filename);


