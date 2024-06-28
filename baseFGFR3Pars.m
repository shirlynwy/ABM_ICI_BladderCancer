function fgfr3 = baseFGFR3Pars()

%% region info
fgfr3.region_type = "concentric_rectangles";
fgfr3.rectangle_compartment_size_x = 1;
fgfr3.rectangle_compartment_size_y = 1;
fgfr3.rectangle_compartment_size_z = 1;
fgfr3.blood_vessels_in_separate_region = true;
fgfr3.is_pk = true;

%% tumor parameters
% https://link.springer.com/content/pdf/10.1023/B:ABME.0000030231.88326.78.pdf
% for a paper that gives a range for this in #/cell (they say the number is 10^4 but then they vary between 10^3 and 10^5 to test the effect). (Filion and Popel.
% 2003. A Reaction-Diffusion Model of Basic Fibroblast Growth Factor
% Interactions with Cell Surface Receptors)
% fgfr3.RT = (8.3056e-9/fgfr3.V)*1e-3; % in nM (per cell); Kamaldeen's ref says 1e4 # per cell which would be 1.66e-11 nmol instead of 8.3056e-9nmol; 
fgfr3.RT_num = ((1e4/6.02214e23)*1e9)*1e15; % in nM * um^3 (per cell); Kamaldeen's ref says 1e4 # per cell which would be 1.66e-11 nmol instead of 8.3056e-9nmol; 
fgfr3.RT_sigma = 0; % SD of RT on tumor cells
fgfr3.tum_prolif_up = 0.8571428571; % pars.prolif_rate/4; % (twice the) max increase to prolif probability for tumor cells with FGFR3 signaling; twice because 0<=phiD=D_A/RT<=0.5 (per day)
fgfr3.gammaT = 1/6; % ec50 for phiD (on [0 .5]) inhibiting apoptosis

%% dosing parameters
fgfr3.start_day = Inf; % earliest day in simulation to start (may start later if this is a weekend)
fgfr3.days_between = 1; % days between fgfr3 doses
fgfr3.n_doses = 0;
fgfr3.dose_val = 0; % initial concentration of circulating inhibitor in nM (based on Grunewald and the max plasma concentration with 75mg/kg given to a mouse; I have eyeballed the number, extrapolating back to t=0)

%% immune effect parameters
fgfr3.imm_evasion_ec50 = 1/6;

%% reaction parameters

fgfr3.non_mut_RT_factor = .5; % non-mutants have fewer receptors than mutants. this is the the factor to multiply by

% according to Kamaldeen's source (Bing Zhao et al. 2010. Endothelial Cell
% Capture of Heparin-Binding Growth Factors under Flow) this value should
% be 3.2e8 per M per min which would be 460.8 per nM per day
% fgfr3.kf = 1e11; % in nmol^(-1) day^(-1)
fgfr3.kf = 460.8; % in per nM per day

% the dissociation used by Zhao (see ref for kf) would be 403.2 days^-1,
% which is in the range Kamaldeen gives. If we use the value used by Filion
% (see ref for RT), we would get 69.12 days^-1.
% fgfr3.kr = 864; % dissociation rate in days^-1
fgfr3.kr = 400; % dissociation rate in days^-1

% exactly matches both Zhao and Filion (see above refs) when converting
% their value of .078 min^-1
fgfr3.kp = 112.32; % recycling rate of dimers in days^-1

% anti-FGFR3 related
% Grunewald says Kd is 7.8nM for rogaratinib with FGFR3. I will just make
% sure these have this ratio and that it goes quickly so as to achieve
% equilibrium within a few minutes. Will reference a paper by Tassa 2010
% (SMI_reactionrates.pdf in my Papers folder) for typical SMI
% association/dissociation rates.
fgfr3.k_on_R = 1e2; % inhibitor binding to monomer (in per nM per day); originally 1.28e5 from Kamaldeen
fgfr3.k_off_R = 7.8e2; % unbinding of inhibitor and monomer (in per day); orginally 95.04 from Kamaldeen
fgfr3.k_on_D = 1e2; % inhibitor binding to active dimer (in per nM per day); originally 1.28e5 from Kamaldeen
fgfr3.k_off_D = 7.8e2; % unbinding of inhibitor and active dimer (in per day); orginally 95.04 from Kamaldeen

% fgfr3.aFGFR3_influx = 576.74; % rate of inhibitor diffusing into tumor; from Kamaldeen's K12 with a 75mg/kg dose
% fgfr3.aFGFR3_degradation = 1.2477; % combo rate of anti-FGFR3 diffusing out of tumor + rate of degradation within tumor; from Kamaldeen's K21 with a 75mg/kg dose

% fgfr3.aFGFR3_influx = .48; % rate of inhibitor diffusing across blood vessel walls into and out of tumor (per day) (based on UM_ABM/FGFR_ODE/Test_PK_2Compartment.m
% fgfr3.aFGFR3_influx = 2; % rate of inhibitor diffusing across blood vessel walls into and out of tumor (per day) (assuming per
fgfr3.aFGFR3_influx = 10; % the rate that concentration is accumulated in periphery (in per day); this assumes that the volume of the periphery is exactly 1L

fgfr3.aFGFR3_degradation = 14.4; % combo rate of anti-FGFR3 diffusing out of tumor + rate of degradation within tumor (per day); formerly 1.2477 from Kamaldeen's K21 with a 75mg/kg dose
% fgfr3.aFGFR3_degradation = 1.44; % combo rate of anti-FGFR3 diffusing out of tumor + rate of degradation within tumor (per day); this goes with  of 1e-7 to get the same decay length
% this value is chosen in conjunction with the diffusion to get a length
% scale of sqrt(D/lambda) ~ 300\micro\meter which is what Macklin used as
% the length scale of VEGF (I don't know how these two length scales should
% compare, so I'm assuming they should be roughly equivalent).

fgfr3.aFGFR3_sysdecay = log(2)/.125; % decay rate of drug from system

% diffusion rate (m^2/day) (cf. https://www.sciencedirect.com/topics/nursing-and-health-professions/diffusivity#:~:text=Generally%2C%20the%20value%20of%20diffusion,2%2Fs%20is%20most%20common. 
                        % where they say 1e-8cm^2/s is common for intermediate
                        % molecular weight drugs (150-400 Da) through
                        % flexible polymers (rogaratinib is 466 Da and BGJ
                        % is 560 Da).
fgfr3.aFGFR3_diffusion = 1e6; % diffusion coefficient for fgfr3 inhibitor in TME (in um^2/d); this is 10x the value used by Macklin for diffusion of VEGF (which is 100x the molecular weight of rogaratinib); see below for possible ranges for similarly sized molecules; also used in the global ode;
% fgfr3.aFGFR3_diffusion = 1e5; % diffusion coefficient for fgfr3 inhibitor in TME (in um^2/d); this is smaller than the above value to help make pde more stable. degradation also changed to maintain decay length

% (cf. https://www.sciencedirect.com/science/article/pii/S0079656519300585
%  they show diffusion coefficients of three molecules between 130-320 Da.
%  These seem to span 1 order of magnitude (4e2-12e2 in um^2/s). These
%  correspond to >3.5e7 um^2/d (and <1e8 um^2/d) though they are smaller molecules than our
%  FGFR3 inhibitors. This also means that when varying the diffusion
%  coefficient, 2 orders of magnitude in either direction should be plenty
