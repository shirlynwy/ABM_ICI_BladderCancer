function plot_pars = basePlotPars()

%% Plotting properties

plot_pars.plot_fig = false;
plot_pars.plot_location = false; % don't bother plotting where all the cells are. it takes a long time
plot_pars.make_movie = false;
plot_pars.plot_every = 1;
plot_pars.plot_offset = 0;
plot_pars.use_carveout = false; % whether or not to carveout part of the sphere
