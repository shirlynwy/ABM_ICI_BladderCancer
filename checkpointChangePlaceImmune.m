function M = checkpointChangePlaceImmune(M,j)

region_gains = accumarray(M.immunes(j,M.I.region_checkpoint),1,[M.checkpoint.n_regions,1]);

gain_log = region_gains>0;

M.checkpoint.PD1(gain_log) = (M.checkpoint.PD1(gain_log) .* M.checkpoint.volumes.pd1(gain_log) + M.checkpoint.pd1_on_immune * region_gains(gain_log) );
M.checkpoint.PD1aPD1(gain_log) = M.checkpoint.PD1aPD1(gain_log).*M.checkpoint.volumes.pd1(gain_log);
M.checkpoint.volumes.pd1(gain_log) = M.checkpoint.volumes.pd1(gain_log) + region_gains(gain_log);

M.checkpoint.PD1(gain_log) = M.checkpoint.PD1(gain_log)./M.checkpoint.volumes.pd1(gain_log);
M.checkpoint.PD1aPD1(gain_log) = M.checkpoint.PD1aPD1(gain_log)./M.checkpoint.volumes.pd1(gain_log);
