clearvars;
tic;
% addpath("~/Documents/MATLAB/myfunctions/")
is_cpactive = 0;
isplot = 1;
issave = 0;
%%
M = allBaseParameters();

if issave
    M.save_pars.dt = 0.5; % set to Inf to not save anything; otherwise dt is in days
else
    M.save_pars.dt = Inf;
end

M.save_pars.batch_name = strcat('test_', string(datetime("now","Format","yyMMddHH")));
%% set up
M.setup.prop_ha0 = 0.5; % initial proportion HA
M.setup.N0 = 20;
M.setup.NI0 = 0;
M.setup.censor_date = 20;
M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;

%% basic parameters
M.pars.max_dt = 1/4 / 24; % number of days per tumor step
M.pars.max_dt_Imm = 1/8/24;%number of days per immun step
M.pars.prolif_rate = 2.6 ; %proliferation rate of tumor cells 1 per day
M.pars.occmax =  20; %pairs of prolif rate and occmax from SMorePars: (2.36, 20.3), (3.34, 19)
M.pars.max_tumor_size = 60000; 
M.pars.min_tumor_size = 2;
M.pars.max_tumor_boundary = 36; %if number of tumor cells on the boundary exceeds this, stop simulation

%% immune parameters
M.immune_pars.steps_per_move = 4;
M.immune_pars.aicd_rate = 0;
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day

%immune kill rates
M.immune_pars.slow_kill_rate = 1 / (2/24); % rate of slow kill in days^-1 (should be once per 2 hours)
M.immune_pars.fast_kill_rate = 1 / (.5/24); % rate of fast kill in days^-1 (should be once per half hour)
%%% fast killing probability for high and low antigen cells
% M.immune_pars.p1 = 0.93; % high antigen cells
% M.immune_pars.p2 = 0.33; % low antigen cells
%%% immune cell death due to interaction with tumor cells
% M.immune_pars.death_by_fast = 0.01; % number of immune cells per day
% M.immune_pars.death_by_slow = 0.01; % number of immune cells per day

%recruitment like ODE
M.immune_pars.prolif_rate = 0; % immune proliferation rate per day
M.immune_pars.immune_recruit_rate2 = 8; % immune cells recruited per day

%set to so immune cells are not deactivated in conjugation (b/c not in ODE)
M.immune_pars.deactivation_rate_pd1 = 0;

%ISF stimulated CTL proliferation
M.immune_pars.isf_prolif_max = 0.15; % maximum antigen stimulated proliferation rate

%% checkpoint blockade therapy parameters, no dosing in this ABM
M.checkpoint.start_day = Inf;
M.checkpoint.dose_val = 0;

%% checkpoint blocked
if ~is_cpactive
    % M.checkpoint.kf_pd1_pdl1 = 0; %can't set this to 0 bc can't /0
    M.checkpoint.kr_pd1_pdl1 = 0;
    M.checkpoint.pd1_on_immune_nmols = 0;
    M.checkpoint.pdl1_on_tumor_nmols = 0;
end

%% plotting
if isplot
    M.plot_pars.plot_fig = true;
    M.plot_pars.plot_location = true;
    M.plot_pars.use_carveout = true;
    plot_pars.plot_every = 1;
end

M = simPatient(M);

toc;

