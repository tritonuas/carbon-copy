function w_fuse = weight_fuse_analytic(density_carbon_epoxy, num_plies, num_bulkheads, t_bulkhead, area_fraction_bulkhead, density_plywood, fuse_height, fuse_width, fuse_length, fudgefactor)

    bulkhead_area_density = density_plywood*num_bulkheads*area_fraction_bulkhead*t_bulkhead; % kg/m^2   
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies;
    
    w_fuse = carbon_epoxy_area_density(2*fuse_length*fuse_width+2*fuse_length*fuse_height+2*fuse_height*fuse_width) + fuse_height*fuse_width*bulkhead_area_density;
    w_fuse = fudgefactor*w_fuse;
end
