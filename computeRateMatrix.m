function rate_matrix = computeRateMatrix(M)

rate_matrix = zeros(M.NT,2);

rate_matrix(:,1) = M.pars.prolif_rate + M.fgfr3.tum_prolif_up * M.phiD;
rate_matrix(:,2) = M.pars.apop_rate ./ (1+M.phiD/M.fgfr3.gammaT);
