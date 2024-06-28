function M = saveInitialModelData(M)

grid_size = M.grid.size;
setup = M.setup;
pars = M.pars;
plot_pars = M.plot_pars;
flags = M.flags;
immune_pars = M.immune_pars;
immune_stimulatory_factor_pars = M.immune_stimulatory_factor_pars;
save_pars = M.save_pars;
fgfr3 = M.fgfr3;
checkpoint = M.checkpoint;

save(sprintf("data/%s/%s/output_constants",M.save_pars.batch_name, M.save_pars.sim_identifier),'-regexp','[^M]','-v7.3')

M = saveModelData(M);
