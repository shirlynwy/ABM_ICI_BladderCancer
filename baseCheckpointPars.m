function checkpoint = baseCheckpointPars()

%% region info
checkpoint.region_type = "concentric_rectangles";
checkpoint.rectangle_compartment_size_x = 1;
checkpoint.rectangle_compartment_size_y = 1;
checkpoint.rectangle_compartment_size_z = 1;
checkpoint.blood_vessels_in_separate_region = true;
checkpoint.is_pk = true;

%% dosing parameters
checkpoint.start_day = 1; % earliest day in simulation to start (may start later if this is a weekend)
checkpoint.days_between = 3; % days between fgfr3 doses
checkpoint.n_doses = 50;
checkpoint.dose_val = 1000; % initial concentration of circulating inhibitor in nM

%% pd1 parameters
checkpoint.pd1_on_immune_nmols = ((3096 ... % mean number of PD1 on a CTL (see Katie's ODE paper for this number and refs) Cheng, Xiaoxiao, et al. "Structure and interactions of the human programmed cell death 1 receptor." Journal of Biological Chemistry 288.17 (2013): 11771-11785.
             /6.02214076e23) * 1e9); % nano mols (nmols) of PD1 on a CTL

%% pdl1 parameters
checkpoint.pdl1_on_tumor_nmols = 10 * ... % assume tumors have 10x expression of pdl1 compared to ctls
             ((9282 ... % mean number of PDL1 on a CTL (see Katie's ODE paper for this number and refs) Cheng, Xiaoxiao, et al. "Structure and interactions of the human programmed cell death 1 receptor." Journal of Biological Chemistry 288.17 (2013): 11771-11785.
             /6.02214076e23) * 1e9); % nano mols (nmols) of PD1 on a CTL
             

%% apd1 parameters
checkpoint.diffusivity_apd1 = 10 * 60*60*24; % diffusion coefficient for anti-pd1 (um^2/day) (from Thurber paper)
checkpoint.degradation_apd1 = 1; % degradation rate for anti-pd1
checkpoint.systemic_elimination_apd1 = log(2)/3; % systemic elimination rate for anti-pd1 in mouse

%% apdl1 parameters
checkpoint.diffusivity_apdl1 = 10 * 60*60*24; % diffusion coefficient for anti-pd1 (um^2/day) (from Thurber paper)
checkpoint.degradation_apdl1 = 1; % degradation rate for anti-pdl1
checkpoint.systemic_elimination_apdl1 = 0; % not using apdl1 yet

%% reaction parameters

checkpoint.kf_pd1_pdl1 = 100;
checkpoint.kr_pd1_pdl1 = 8.2e5;

checkpoint.kf_pd1_apd1 = 1e3;
checkpoint.kr_pd1_apd1 = 1.45e3;

checkpoint.kf_pdl1_apdl1 = 0; % not using apdl1 yet, so no need to set this
checkpoint.kr_pdl1_apdl1 = 0; % not using apdl1 yet, so no need to set this

% (https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12143) Bajaj et al. 2016. Model?Based Population Pharmacokinetic Analysis of Nivolumab in Patients With Solid Tumors
checkpoint.influx_apd1 = log(2)/(32.5/24); % rate of influx of anti-pd1 into TME (the log(2)/(32.5/24) value is from Bajaj (above) where they say the distribution phase half life is 32.5 h.
checkpoint.efflux_apd1 = 0; % efflux of anti-pd1 out of TME (zero because pd1 crosses the capillary walls by convection, which I can safely (I believe) assume is always into the periphery

checkpoint.influx_apdl1 = 0; % not including apdl1 yet
checkpoint.efflux_apdl1 = 0; % efflux of anti-pdl1 out of TME (zero because pd1 crosses the capillary walls by convection, which I can safely (I believe) assume is always into the periphery




