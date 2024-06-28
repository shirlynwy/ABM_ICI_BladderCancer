function M = initializeImmuneStimulatoryFactor(M)

M.immune_stimulatory_factor.concentration = zeros(M.grid.size);

for i = 1:M.NT

    [xx,yy,zz,xind,yind,zind] = computeISFIndices(M,M.tumors(i,M.I.subs));
    
    if M.tumors(i,M.I.type)==0
        M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) + M.immune_pars.low_antigen_isf_factor * M.immune_stimulatory_factor.stamp(xind,yind,zind);
    else
        M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) + M.immune_stimulatory_factor.stamp(xind,yind,zind);
    end
end