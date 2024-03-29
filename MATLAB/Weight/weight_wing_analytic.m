%%Authors: Allyson Chen & Kevin Vo
%%Emails:  aechen@ucsd.edu    k1vo@ucsd.edu
%%Calculates the weight of the the wing

%%Variable Legend
%%num_spar = number of spars
%%num_plies = number of plies
%%spar_width = width of the spar
%%S = wing area
%%wingspan = wingspan
%%density_balsa = density of balsa wood
%%density_carbon_epoxy = density of carbon epoxy composite
%%density_divinycell = density of divinycell (aerospace foam)
%%t_tip = thickness of wing tip
%%t_root = thickness of wing root
%%t_divinycell = thickness of divinycell (aerospace foam)
%%fudgefactor = fudgefactor


function w_wing = weight_wing_analytic(S, num_spar,  spar_width, density_balsa, t_divinycell, density_divinycell, num_plies, density_carbon_epoxy, t_tip, t_root, wingspan, fudge_factor)
    
    t_carbon_epoxy = 0.0003;
    carbon_epoxy_area_density = density_carbon_epoxy*num_plies*t_carbon_epoxy; % kg/m^2 
    divinycell_area_density = t_divinycell*density_divinycell; % kg/m^2
    spar_length_density = ((t_tip+t_root)/2)*num_spar*density_balsa*spar_width; % kg/m
    
    
    w_wing = 2*S*(carbon_epoxy_area_density+divinycell_area_density) + wingspan*spar_length_density; % multiply area density by 2 to account for the top and bottom wing skins
    w_wing = 9.81*w_wing; % convert from kg to N
    w_wing = w_wing*fudge_factor;
end 
