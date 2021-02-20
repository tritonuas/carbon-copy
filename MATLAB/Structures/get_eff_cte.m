%%Finds the effective coefficients of thermal expaansions for an entire
%%laminate.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%%@param material_props a vector defined as [E11 E22 G12 Nu12]
%%@param thetas a vector of the angles of each lamina (in order)
%%@param rad_or_deg a string of 'rad' or 'deg' specifying units for theta
%%@param thicknesses a vector of the thickness of each ply
%%@param cte_mat_local the mat of cte values, each ply occupying a col
function [cte_effective] = get_eff_cte(material_props, thetas, ...
    rad_or_deg, thicknesses, cte_mat_local)

delta_T = 1;

thermal_loading = get_thermal_load(thetas, rad_or_deg,...
    material_props, thicknesses, delta_T, cte_mat_local);

%calculating the ABD matrix (essentially laminate stiffness matrix)
ABD = laminate_to_ABD_planar_ortho(thetas, rad_or_deg, thicknesses, ...
    material_props);

%calculating mid-plane strains and curvatures (nondimenionalized
%displacement in the 3 translation and 3 rotation axes)
%mid_strains_and_curvatures = AlphaBetaDelta*loading;
mid_strains_and_curvatures = ABD\thermal_loading;

%normalizing by delta_T to get the effective cte
cte_effective = mid_strains_and_curvatures(1:3)/delta_T;