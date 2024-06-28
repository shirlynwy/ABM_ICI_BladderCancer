function M = updateFGFR3(M)

if M.fgfr3.circ==0 && all(M.fgfr3.concentration==0) && M.fgfr3.RT_sigma==0 % then new cells can be updated predictably; this cannot be done once the M.tumors have random amounts of receptors
    M.tumors(M.mut_log,M.I.tumor_receptors) = repmat(M.fgfr3.all_receptors_without_inhibitor,sum(M.mut_log),1);
    M.tumors(~M.mut_log,M.I.tumor_receptors) = repmat([M.fgfr3.RT*M.fgfr3.non_mut_RT_factor,0,0,0,0,0],sum(~M.mut_log),1);
else
    %% global method for update
    if M.fgfr3.circ==0 && all(M.fgfr3.concentration==0) && all(M.tumors(:,M.I.tumor_receptors(3:end))==0,'all')
        region_props = [];
    else
        region_props = determinePropNeighbors(M);
    end

    if M.NT > 0 % to allow for just letting the immune system burn itself out
        R = zeros(length(M.I.tumor_receptors),3,M.fgfr3.n_regions); % D1: various receptors/complexes; D3: Regions as determined by blood vessels; D2: nonmut, mut, non-tumor (free or immune)
        nonmut_assign = true(M.fgfr3.n_regions,1);
        mut_assign = true(M.fgfr3.n_regions,1);
        for ri = M.fgfr3.n_regions:-1:1

            nonmut_recs = M.tumors(~M.mut_log & M.tumors(:,M.I.region_fgfr3)==ri,M.I.tumor_receptors); % receptor concentrations of non-mutants in region ri
            if isempty(nonmut_recs) % then there weren't any nonmuts in this region
                nonmut_recs = zeros(1,length(M.I.tumor_receptors));
                nonmut_assign(ri) = false; % don't bother trying to assigning this later
            end
            RT_nonmut = nonmut_recs*[1;2;0;1;2;2]; % total number of FGFR on each mutant tumor cell
            [z_nonmut{ri},mean_nonmut(ri),sigma_nonmut(ri)] = zscore(RT_nonmut,1,1);

            mut_recs = M.tumors(M.mut_log & M.tumors(:,M.I.region_fgfr3)==ri,M.I.tumor_receptors); % receptor concentrations of mutants in region ri
            if isempty(mut_recs) % then there weren't any muts in this region
                mut_recs = zeros(1,length(M.I.tumor_receptors));
                mut_assign(ri) = false; % don't bother trying to assigning this later
            end
            RT_mut = mut_recs*[1;2;0;1;2;2]; % total number of FGFR on each mutant tumor cell
            [z_mut{ri},mean_mut(ri),sigma_mut(ri)] = zscore(RT_mut,1,1);

            R(:,1,ri) = mean(nonmut_recs,1); % take average value in this region to feed into global method
            R(:,2,ri) = mean(mut_recs,1); % take average value in this region to feed into global method
        end
        R(3,3,:) = M.fgfr3.concentration; % in non-tumor sites (position 3 in dimension 2), only state variable that exists (and thus matters) is the free anti-fgfr3 (Position 3 in dimension 1), which is in the M.fgfr3.concentration vector
        R(isnan(R)) = 0; % for any time there was the mean of an empty array taken

        Rnew = globalFGFR3ODE(R(:),size(R),M,region_props);
        M.fgfr3.concentration = reshape(Rnew(3,3,:),[],1); % update free anti-fgfr3
        M.fgfr3.circ = M.fgfr3.circ * exp(-M.fgfr3.aFGFR3_sysdecay*M.dt); % update circulating afgfr3
        for ri = 1:M.fgfr3.n_regions
            if nonmut_assign(ri)
                M.tumors(~M.mut_log & M.tumors(:,M.I.region_fgfr3)==ri,M.I.tumor_receptors) = Rnew(:,1,ri)'.*(1+z_nonmut{ri}*sigma_nonmut(ri)/mean_nonmut(ri));
            end
            if mut_assign(ri)
                M.tumors(M.mut_log & M.tumors(:,M.I.region_fgfr3)==ri,M.I.tumor_receptors) = Rnew(:,2,ri)'.*(1+z_mut{ri}*sigma_mut(ri)/mean_mut(ri));
            end
        end
    end
end % end of how to update M.tumors before cell fate decisions

M.phiD = M.tumors(:,M.I.tumor_receptors(2))/M.fgfr3.RT; % set  phiD value for all tumor cells

M.tracked.phiD_mean(M.i) = mean(M.phiD(M.mut_log));
M.tracked.phiD_std(M.i) = std(M.phiD(M.mut_log));