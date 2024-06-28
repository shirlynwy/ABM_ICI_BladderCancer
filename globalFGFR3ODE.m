function R = globalFGFR3ODE(R,szR,M,region_props)

if isempty(region_props)
    region_props = zeros([szR(2:3),szR(2:3)]);
end

options = odeset('Nonnegative',1:numel(R));
sol = ode15s(@(t,x) globalODE(t,x,szR,M,region_props),[0,M.dt],R(:),options);
R = reshape(sol.y(:,end),szR); % store the end value to output and reshape


function dx = globalODE(t,x,szR,M,region_props)

% x = [ RF;DA;C;R_F_C;D_C_C;D_A_C;ambient drug;RF on non-mutants;C near non-mutants;R_F_C on non-mutants ] molecules all given by average concentrations

% D1: various receptors/complexes; D2: nonmut, mut, non-tumor (free or immune); D3: Regions as determined by blood vessels;
x = reshape(x,szR); % reshape x so that the dimensions have their above meaning
dx = zeros(szR); % initially comput dx/dt in this array form and then string it out at the end

dimerization = [0,1,0].*(-M.fgfr3.kf*x(1,:,:).^2 + (M.fgfr3.kr+M.fgfr3.kp)*x(2,:,:)); % dimerization and backwards reaction
monomer_inhibitor = - M.fgfr3.k_on_R*x(1,:,:).*x(3,:,:) + M.fgfr3.k_off_R*x(4,:,:); % FGFR3 reacting with aFGFR3
complex_dimerization = [0,1,0].*(- M.fgfr3.kf*x(4,:,:).^2 + (M.fgfr3.kr+M.fgfr3.kp)*x(5,:,:)); % dimerization of these complexes (and the backwards reaction)
dimer_inhibitor = - M.fgfr3.k_on_D*x(2,:,:).*x(3,:,:) + M.fgfr3.k_off_D*x(6,:,:); % aFGFR3 reacting with dimers

free_inhibitor = x(3,:,:);
conc_diffs = free_inhibitor(:)-free_inhibitor(:)'; % the first index varying varies the positive one (which is the one concentration is flowing from)
conc_diffs = reshape(conc_diffs,[szR(2:3),szR(2:3)]);
diffusion = reshape(sum(region_props.*conc_diffs,1:2),size(free_inhibitor)); % weight the concentration differences by the proportion of neighbors of this (type, region) combo that are of the source (type, region) combo; then add these all up to determine amount of diffusion into this (type, region) center combo
bv_interface = zeros(size(free_inhibitor)); % set up where drug can influx
bv_interface(1,:,M.fgfr3.regions_with_blood_vessels) = 1; % inhibitor only enters in certain regions

% free FGFR3
dx(1,:,:) = 2*dimerization + monomer_inhibitor; 

% dimerized FGFR3
dx(2,:,:) = -dimerization + dimer_inhibitor;
    
% free aFGFR3
dx(3,:,:) = M.fgfr3.aFGFR3_influx*(M.fgfr3.circ*exp(-M.fgfr3.aFGFR3_sysdecay*t)-x(3,:,:)).*bv_interface... % aFGFR3 in circulation flowing into TME
    - M.fgfr3.aFGFR3_degradation.*x(3,:,:)... % aFGFR3 leaving the TME
    + monomer_inhibitor ...% aFGFR3-FGFR3-dimer reactions
    + dimer_inhibitor ... % aFGFR3-FGFR3-dimer reactions
    + 6*M.fgfr3.aFGFR3_diffusion*diffusion/(M.pars.cell_width^2); % drug moving within TME

% aFGFR3-FGFR3 complex
dx(4,:,:) = -monomer_inhibitor + 2*complex_dimerization; 

% aFGFR3-FGFR3 dimer
dx(5,:,:) = -complex_dimerization; % dimerization of aFGFR3-FGFR3 complexes (and the backwards reaction)

% aFGFR3-dimer complex
dx(6,:,:) = -dimer_inhibitor; % aFGFR3 reacting with dimers

dx = dx(:);
