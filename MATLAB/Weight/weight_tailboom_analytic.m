%%Authors: Allyson Chen & Kevin Vo
%%Emails:  aechen@ucsd.edu    k1vo@ucsd.edu
%%Calculates the weight of the the tailboom

%%Variable Legend
%%num_plies = number of plies
%%density_carbon_epoxy = density of carbon epoxy composite
%%surface_area = surface area of cylinder(rectangle)
%%tailboom_length = length of tailboom(meters)
%%fudgefactor = fudgefactor



function w_tailboom = weight_tailboom_analytic(num_plies, density_carbon_epoxy, tailboom_length, fudge_factor)
    
    t_carbon_epoxy = 0.0003;
    tailboom_radius = 0.0762/2;
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies*t_carbon_epoxy; % kg/m^2 
    
    surface_area = tailboom_length*2*pi*tailboom_radius; 
    w_tailboom = surface_area*carbon_epoxy_area_density; 
    w_tailboom = 9.81*w_tailboom; % convert from kg to N
    w_tailboom = w_tailboom*fudge_factor;
end 