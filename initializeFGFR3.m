function M = initializeFGFR3(M)

M.fgfr3.circ = 0;

M = initializeSubstrateRegions(M,"fgfr3");

% 
% % M.fgfr3.concentration = zeros(M.fgfr3.n_regions,1);
% % 
% % M.fgfr3.neighbors_when_empty = zeros(M.fgfr3.n_regions,'uint32'); % total number of neighbors if lattice is totally empty
% % 
% % warning("Assuming that blood vessels are everywhere for now")
% 
% % M.fgfr3.regions = ones(M.grid.size,'uint8');
% % M.fgfr3.neighbors_when_empty(1,1) = 26*M.V_tot;
% 
% 
% switch M.fgfr3.region_type
%     case "concentric_circles"
%         if length(M.grid.size)==2
%             dist_from_origin = sqrt((((1:M.grid.size(1)) - (M.grid.size(1)+1)/2)').^2 + (((1:M.grid.size(2)) - (M.grid.size(2)+1)/2)).^2);
%         elseif length(M.grid.size)==3
%             dist_from_origin = sqrt((((1:M.grid.size(1)) - (M.grid.size(1)+1)/2)').^2 + (((1:M.grid.size(2)) - (M.grid.size(2)+1)/2)).^2 + ((reshape(1:M.grid.size(3),1,1,[]) - (M.grid.size(3)+1)/2)).^2);
%         end
%         if isfield(M.fgfr3,"radius_transform_fn")
%             dist_from_origin = M.fgfr3.radius_transform_fn(dist_from_origin);
%         end
%         norm_dist_from_origin = dist_from_origin/max(dist_from_origin,[],'all');
%         M.fgfr3.regions = ceil(M.fgfr3.num_concentric_circle_M.fgfr3.regions*norm_dist_from_origin);
% 
%     case "rectangular_grid"
%         assert(all(mod(M.grid.size,M.fgfr3.rectangle_compartment_size)==0))
%         x_additions = repelem((1:(M.grid.size(1)/M.fgfr3.rectangle_compartment_size(1))),M.fgfr3.rectangle_compartment_size(1));
%         y_additions = repelem(((1:(M.grid.size(2)/M.fgfr3.rectangle_compartment_size(2)))-1)*(M.grid.size(1)/M.fgfr3.rectangle_compartment_size(1)),M.fgfr3.rectangle_compartment_size(2));
%         if length(M.grid.size)==2
%             M.fgfr3.regions = x_additions' + y_additions;
%         elseif length(M.grid.size)==3
%             z_additions = repelem(((1:(M.grid.size(3)/M.fgfr3.rectangle_compartment_size(3)))-1)*(prod(M.grid.size(1:2))/prod(M.fgfr3.rectangle_compartment_size(1:2))),M.fgfr3.rectangle_compartment_size(3));
%             M.fgfr3.regions = x_additions' + y_additions + reshape(z_additions,1,1,[]);
%         end
% 
%     case "distance_to_dc"
%         all_inds = reshape(1:prod(M.grid.size),M.grid.size);
%         xx = cell(length(M.grid.size),1);
%         [xx{:}] = ind2sub(M.grid.size,all_inds);
% 
%         yy = cell(length(M.grid.size),1);
%         [yy{:}] = ind2sub(M.grid.size,M.fgfr3.dc_inds);
% 
%         distance_to_dc = zeros(prod(M.grid.size),numel(M.fgfr3.dc_inds));
%         for i = 1:length(M.grid.size)
%             distance_to_dc = distance_to_dc + (xx{i}(:) - yy{i}(:)').^2;
%         end
%         distance_to_dc = sqrt(min(distance_to_dc,[],2));
%         M.fgfr3.regions = reshape(ceil(distance_to_dc/M.fgfr3.distance_delta),M.grid.size);
% 
% 
%     case "concentric_rectangles"
%         M.fgfr3.regions = zeros(M.grid.size);
%         M.fgfr3.n_regions = ceil(min(M.grid.size)/2);
%         for i = 2:M.fgfr3.n_regions
%             if length(M.grid.size)==2
%                 M.fgfr3.regions(i:(end-i+1),i:(end-i+1)) = M.fgfr3.regions(i:(end-i+1),i:(end-i+1)) - 1;
%             else
%                 M.fgfr3.regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) = M.fgfr3.regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) - 1;
%             end
%         end
% 
%     otherwise
%         M.fgfr3.regions = reshape(1:prod(M.grid.size),M.grid.size);
% 
% end
% 
% if isfield(M.fgfr3,"blood_vessels_in_separate_region") && M.fgfr3.blood_vessels_in_separate_region && M.fgfr3.is_pk
% 
%     if M.fgfr3.region_type=="rectangular_grid" && (M.setup.blood_vessel_locations=="floor" || M.setup.blood_vessel_locations=="none")
%         M.fgfr3.regions(M.blood_vessels>0) = M.fgfr3.regions(M.blood_vessels>0)-0.5; % keep them neighboring their original region; could perhaps choose it as +0.5 to keep it tridiagonal? not sure about this
%     elseif (any(M.fgfr3.region_type==["concentric_circles","concentric_rectangles"]) && M.setup.blood_vessel_locations=="outside") || (M.fgfr3.region_type=="concentric_circles" && any(M.setup.blood_vessel_locations==["max_shell","weighted_shell"]))
%         M.fgfr3.regions(M.blood_vessels>0) = M.fgfr3.regions(M.blood_vessels>0)+0.5; % move BV M.fgfr3.regions further from center
%     elseif M.fgfr3.region_type=="distance_to_dc" && M.setup.blood_vessel_locations=="none"
%         % should not need to fix this because there are no blood vessels to
%         % cause adjustments to M.fgfr3.regions
%     else
%         error("Have not yet decided how to handle creating separate M.fgfr3.regions when blood vessels are at %s and M.fgfr3.regions are in %s",M.setup.blood_vessel_locations,M.fgfr3.region_type)
%     end
% 
% end
% 
% % make sure the M.fgfr3.regions are numbered 1:M.fgfr3.n_regions without skipping any
% % numbers so that the region number is the index in the concentration
% % vector
% [temp,~,new_regions] = unique(M.fgfr3.regions);
% M.fgfr3.regions = reshape(new_regions,size(M.fgfr3.regions));
% M.fgfr3.n_regions = length(temp);


M.fgfr3.concentration = zeros(M.fgfr3.n_regions,1);

% M.fgfr3.neighbors_when_empty = accumarray([reshape(M.fgfr3.regions(1:end-1,:,:),[],1),reshape(M.fgfr3.regions(2:end,:,:),[],1);... % count neighbors to left
%             reshape(M.fgfr3.regions(2:end,:,:),[],1),reshape(M.fgfr3.regions(1:end-1,:,:),[],1);... % count neighbors to right
%             reshape(M.fgfr3.regions(:,1:end-1,:),[],1),reshape(M.fgfr3.regions(:,2:end,:),[],1);... % count neighbors in front
%             reshape(M.fgfr3.regions(:,2:end,:),[],1),reshape(M.fgfr3.regions(:,1:end-1,:),[],1);... % count neighbors behind
%             reshape(M.fgfr3.regions(:,:,1:end-1),[],1),reshape(M.fgfr3.regions(:,:,2:end),[],1);... % count neighbors below
%             reshape(M.fgfr3.regions(:,:,2:end),[],1),reshape(M.fgfr3.regions(:,:,1:end-1),[],1);... % count neighbors above
%             reshape(M.fgfr3.regions(1,:,:),[],1),reshape(M.fgfr3.regions(1,:,:),[],1);... % count neighbors left of left edge as same region as center
%             reshape(M.fgfr3.regions(end,:,:),[],1),reshape(M.fgfr3.regions(end,:,:),[],1);... % count neighbors right of right edge as same region as center
%             reshape(M.fgfr3.regions(:,1,:),[],1),reshape(M.fgfr3.regions(:,1,:),[],1);... % count neighbors in front of front edge as same region as center
%             reshape(M.fgfr3.regions(:,end,:),[],1),reshape(M.fgfr3.regions(:,end,:),[],1);... % count neighbors behind back edge as same region as center
%             reshape(M.fgfr3.regions(:,:,1),[],1),reshape(M.fgfr3.regions(:,:,1),[],1);... % count neighbors below bottom edge as same region as center
%             reshape(M.fgfr3.regions(:,:,end),[],1),reshape(M.fgfr3.regions(:,:,end),[],1);... % count neighbors above top edge as same region as center
%             ],1);
% 
% M.fgfr3.regions_with_blood_vessels = unique(M.fgfr3.regions(M.blood_vessels));

