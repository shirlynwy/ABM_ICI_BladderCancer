function M = initializeBloodVessels(M)

switch M.setup.blood_vessel_locations
    case "everywhere"
        M.blood_vessels = true(M.grid.size);
    case "outside"
        M.blood_vessels = true(M.grid.size);
        M.blood_vessels(2:end-1,2:end-1,2:end-1) = false;
end
