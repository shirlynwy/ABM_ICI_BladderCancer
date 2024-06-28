function dx = globalCheckpointODE(x,M)

% this will try to track drugs near reacting, non-reacting, and empty
% locations

dx = zeros(numel(x),1);

aPD1_circulation = x(end);
x = reshape(x(1:end-1),M.checkpoint.sz);
dx_tme = zeros(M.checkpoint.sz);

pd1_apd1_reaction   = -M.checkpoint.kf_pd1_apd1*x(:,1).*x(:,2) + M.checkpoint.kr_pd1_apd1*x(:,3);

blood_exchange_apd1 = zeros(M.checkpoint.n_regions,1);
blood_exchange_apd1(M.checkpoint.regions_with_blood_vessels)  = M.checkpoint.influx_apd1 * aPD1_circulation  - M.checkpoint.efflux_apd1 *x(M.checkpoint.regions_with_blood_vessels,2); % inhibitor only enters in certain regions

free_aPD1 = x(:,2);
conc_diffs = free_aPD1-free_aPD1'; % the first index varying varies the positive one (which is the one concentration is flowing from)
conc_diffs = conc_diffs .* M.checkpoint.neighbor_region_props;
diffusion_aPD1 = M.checkpoint.diffusion_factor_apd1*squeeze(sum(conc_diffs,1))';

% pd1
dx_tme(:,1) = pd1_apd1_reaction;

% free anti-pd1
dx_tme(:,2) = M.checkpoint.pd1_proportion .* pd1_apd1_reaction + blood_exchange_apd1 + diffusion_aPD1 - M.checkpoint.degradation_apd1*x(:,2);

% pd1-apd1 complexes
dx_tme(:,3) = -pd1_apd1_reaction;

dx(1:end-1) = dx_tme(:);

% apd1 in circulation
dx(end) = -M.checkpoint.systemic_elimination_apd1*aPD1_circulation;
