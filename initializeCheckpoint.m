function M = initializeCheckpoint(M)

M = initializeSubstrateRegions(M,"checkpoint");
M.tumors(:,M.I.region_checkpoint) = M.checkpoint.regions(M.tumors(:,M.I.ind));
M.immunes(:,M.I.region_checkpoint) = M.checkpoint.regions(M.immunes(:,M.I.ind));

%% concentrations
M.checkpoint.PD1 = M.checkpoint.pd1_on_immune * ones(M.checkpoint.n_regions,1); % PD1 concentration in each region
M.checkpoint.aPD1 = zeros(M.checkpoint.n_regions,1); % aPD1 concentration in each region
M.checkpoint.PD1aPD1 = zeros(M.checkpoint.n_regions,1); % PD1aPD1 concentration in each region
M.checkpoint.aPD1_circulation = 0;

M.checkpoint.sz = [M.checkpoint.n_regions,3];

%% intra region volumes and proportions

M.checkpoint.neighbor_region_props = M.checkpoint.neighbors_when_empty ./ sum(M.checkpoint.neighbors_when_empty,1);
M.checkpoint.volumes.regions = accumarray(M.checkpoint.regions(:),1,[M.checkpoint.n_regions,1]);
M.checkpoint.volumes.pd1 = accumarray(M.immunes(:,M.I.region_checkpoint),1,[M.checkpoint.n_regions,1]); % volumes of each region that has pd1 (i.e. has a ctl)





