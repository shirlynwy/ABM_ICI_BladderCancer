function setup = baseSetup()

setup.start_day = -1; % start on Sunday

setup.grid_size_microns_x = 1000;
setup.grid_size_microns_y = 1000;
setup.grid_size_microns_z = 1000;

setup.censor_date = 20;
setup.N0 = 20;
setup.NI0 = 0;

setup.c = -3.603357085551339; % see Test_findInitialization for code that got to this number
setup.e = 2.986939791722032; % see Test_findInitialization for code that got to this number

setup.prop_ha0 = .5; % initial proportion HA
setup.prop_mut0 = 0; % initial proportion with FGFR3 mutation

%% blood vessels
setup.blood_vessel_locations = "outside";