function M = removeApoptotic(M)

apoptosis_log = M.tumors(:,M.I.event)==2;
apoptosis_ind = find(apoptosis_log);

if ~isempty(apoptosis_ind)

    for i = 1:length(apoptosis_ind)
        [xx,yy,zz,xind,yind,zind] = computeISFIndices(M,M.tumors(apoptosis_ind(i),M.I.subs));
        if M.tumors(apoptosis_ind(i),M.I.type)==0
            M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) - M.immune_pars.low_antigen_isf_factor * M.immune_stimulatory_factor.stamp(xind,yind,zind);
        else
            M.immune_stimulatory_factor.concentration(xx,yy,zz) = M.immune_stimulatory_factor.concentration(xx,yy,zz) - M.immune_stimulatory_factor.stamp(xind,yind,zind);
        end
    end

    region_empty_volume_gains = accumarray(M.tumors(apoptosis_ind,M.I.region_fgfr3),1,[M.fgfr3.n_regions,1]);
    region_free_fgfr3_gains = accumarray(M.tumors(apoptosis_ind,M.I.region_fgfr3),M.tumors(apoptosis_ind,M.I.tumor_receptors(3)),[M.fgfr3.n_regions,1]);

    empty_region_volume = accumarray(M.fgfr3.regions(:),M.L==0,[M.fgfr3.n_regions,1]);
    M.fgfr3.concentration = (empty_region_volume .* M.fgfr3.concentration + region_free_fgfr3_gains) ./ (empty_region_volume + region_empty_volume_gains);
    M.fgfr3.concentration(empty_region_volume==0 & region_empty_volume_gains==0) = 0;
    
    M.L(M.tumors(apoptosis_ind,M.I.ind)) = 0;

    M.tumors(apoptosis_ind,:)=[]; % get rid of dead cells
end
