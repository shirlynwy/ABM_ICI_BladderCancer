function M = initializeSubstrateRegions(M,name)

switch M.(name).region_type
    case "concentric_circles"
        if length(M.grid.size)==2
            dist_from_origin = sqrt((((1:M.grid.size(1)) - (M.grid.size(1)+1)/2)').^2 + (((1:M.grid.size(2)) - (M.grid.size(2)+1)/2)).^2);
        elseif length(M.grid.size)==3
            dist_from_origin = sqrt((((1:M.grid.size(1)) - (M.grid.size(1)+1)/2)').^2 + (((1:M.grid.size(2)) - (M.grid.size(2)+1)/2)).^2 + ((reshape(1:M.grid.size(3),1,1,[]) - (M.grid.size(3)+1)/2)).^2);
        end
        if isfield(M.(name),"radius_transform_fn")
            dist_from_origin = M.(name).radius_transform_fn(dist_from_origin);
        end
        norm_dist_from_origin = dist_from_origin/max(dist_from_origin,[],'all');
        M.(name).regions = ceil(M.(name).num_concentric_circle_regions*norm_dist_from_origin);

    case "rectangular_grid"
        assert(all(mod(M.grid.size,M.(name).rectangle_compartment_size)==0))
        x_additions = repelem((1:(M.grid.size(1)/M.(name).rectangle_compartment_size(1))),M.(name).rectangle_compartment_size(1));
        y_additions = repelem(((1:(M.grid.size(2)/M.(name).rectangle_compartment_size(2)))-1)*(M.grid.size(1)/M.(name).rectangle_compartment_size(1)),M.(name).rectangle_compartment_size(2));
        if length(M.grid.size)==2
            M.(name).regions = x_additions' + y_additions;
        elseif length(M.grid.size)==3
            z_additions = repelem(((1:(M.grid.size(3)/M.(name).rectangle_compartment_size(3)))-1)*(prod(M.grid.size(1:2))/prod(M.(name).rectangle_compartment_size(1:2))),M.(name).rectangle_compartment_size(3));
            M.(name).regions = x_additions' + y_additions + reshape(z_additions,1,1,[]);
        end

    case "distance_to_dc"
        all_inds = reshape(1:prod(M.grid.size),M.grid.size);
        xx = cell(length(M.grid.size),1);
        [xx{:}] = ind2sub(M.grid.size,all_inds);

        yy = cell(length(M.grid.size),1);
        [yy{:}] = ind2sub(M.grid.size,M.(name).dc_inds);

        distance_to_dc = zeros(prod(M.grid.size),numel(M.(name).dc_inds));
        for i = 1:length(M.grid.size)
            distance_to_dc = distance_to_dc + (xx{i}(:) - yy{i}(:)').^2;
        end
        distance_to_dc = sqrt(min(distance_to_dc,[],2));
        M.(name).regions = reshape(ceil(distance_to_dc/M.(name).distance_delta),M.grid.size);


    case "concentric_rectangles"
        M.(name).regions = zeros(M.grid.size);
        M.(name).n_regions = ceil(min(M.grid.size)/2);
        for i = 2:M.(name).n_regions
            if length(M.grid.size)==2
                M.(name).regions(i:(end-i+1),i:(end-i+1)) = M.(name).regions(i:(end-i+1),i:(end-i+1)) - 1;
            else
                M.(name).regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) = M.(name).regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) - 1;
            end
        end

    otherwise
        M.(name).regions = reshape(1:prod(M.grid.size),M.grid.size);

end

if isfield(M.(name),"blood_vessels_in_separate_region") && M.(name).blood_vessels_in_separate_region && M.(name).is_pk

    if M.(name).region_type=="rectangular_grid" && (M.setup.blood_vessel_locations=="floor" || M.setup.blood_vessel_locations=="none")
        M.(name).regions(M.blood_vessels>0) = M.(name).regions(M.blood_vessels>0)-0.5; % keep them neighboring their original region; could perhaps choose it as +0.5 to keep it tridiagonal? not sure about this
    elseif (any(M.(name).region_type==["concentric_circles","concentric_rectangles"]) && M.setup.blood_vessel_locations=="outside") || (M.(name).region_type=="concentric_circles" && any(M.setup.blood_vessel_locations==["max_shell","weighted_shell"]))
        M.(name).regions(M.blood_vessels>0) = M.(name).regions(M.blood_vessels>0)+0.5; % move BV M.(name).regions further from center
    elseif M.(name).region_type=="distance_to_dc" && M.setup.blood_vessel_locations=="none"
        % should not need to fix this because there are no blood vessels to
        % cause adjustments to M.(name).regions
    else
        error("Have not yet decided how to handle creating separate M.(name).regions when blood vessels are at %s and M.(name).regions are in %s",M.setup.blood_vessel_locations,M.(name).region_type)
    end

end

% make sure the M.(name).regions are numbered 1:M.(name).n_regions without skipping any
% numbers so that the region number is the index in the concentration
% vector
[temp,~,new_regions] = unique(M.(name).regions);
M.(name).regions = reshape(new_regions,size(M.(name).regions));
M.(name).n_regions = length(temp);

n_region_bytes = ceil(log2(M.(name).n_regions+1));
if n_region_bytes <= 8
    M.(name).regions = uint8(M.(name).regions);
elseif n_region_bytes <= 16
    M.(name).regions = uint16(M.(name).regions);
elseif n_region_bytes <= 32
    M.(name).regions = uint32(M.(name).regions);
elseif n_region_bytes <= 64
    M.(name).regions = uint64(M.(name).regions);
else
    error("too many %s regions to store!",name)
end

M.(name).neighbors_when_empty = accumarray([reshape(M.(name).regions(1:end-1,:,:),[],1),reshape(M.(name).regions(2:end,:,:),[],1);... % count neighbors to left
            reshape(M.(name).regions(2:end,:,:),[],1),reshape(M.(name).regions(1:end-1,:,:),[],1);... % count neighbors to right
            reshape(M.(name).regions(:,1:end-1,:),[],1),reshape(M.(name).regions(:,2:end,:),[],1);... % count neighbors in front
            reshape(M.(name).regions(:,2:end,:),[],1),reshape(M.(name).regions(:,1:end-1,:),[],1);... % count neighbors behind
            reshape(M.(name).regions(:,:,1:end-1),[],1),reshape(M.(name).regions(:,:,2:end),[],1);... % count neighbors below
            reshape(M.(name).regions(:,:,2:end),[],1),reshape(M.(name).regions(:,:,1:end-1),[],1);... % count neighbors above
            reshape(M.(name).regions(1,:,:),[],1),reshape(M.(name).regions(1,:,:),[],1);... % count neighbors left of left edge as same region as center
            reshape(M.(name).regions(end,:,:),[],1),reshape(M.(name).regions(end,:,:),[],1);... % count neighbors right of right edge as same region as center
            reshape(M.(name).regions(:,1,:),[],1),reshape(M.(name).regions(:,1,:),[],1);... % count neighbors in front of front edge as same region as center
            reshape(M.(name).regions(:,end,:),[],1),reshape(M.(name).regions(:,end,:),[],1);... % count neighbors behind back edge as same region as center
            reshape(M.(name).regions(:,:,1),[],1),reshape(M.(name).regions(:,:,1),[],1);... % count neighbors below bottom edge as same region as center
            reshape(M.(name).regions(:,:,end),[],1),reshape(M.(name).regions(:,:,end),[],1);... % count neighbors above top edge as same region as center
            ],1);

M.(name).regions_with_blood_vessels = unique(M.(name).regions(M.blood_vessels));

