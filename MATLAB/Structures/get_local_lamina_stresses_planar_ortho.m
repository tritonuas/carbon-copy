%%Takes in inputs of material and laminate properties and calculates the
%%local stresses at the top and bottom of each ply (planar).
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%%@param material_props a vector defined as [E11; E22; G12; Nu12]
%%@param thetas a vector of the angles of each lamina (in order)
%%@param rad_or_deg a string of 'rad' or 'deg' specifying units for theta
%%@param thicknesses a vector of the thickness of each ply
%%@param mech_loading the column vector of loading [Nx;Ny;Nxy;Mx;My;Mxy]
%%@param delta_T the change in temperature
%%@param cte_mat_local the matrix of cte for each ply in a local frame
%%of shape [cte1;cte2;cte12] by 
%%@return local_stresses_bot the local stresses at the bottom of plies
%%@return local_stresses_top the local stresses at the top of plies
%%@return z_all a vector of the z-coords of the top and bot of all plies
%%@return z_bot a vector of the z-coordinates of the bottom of the plies
%%@return z_top a vector of the z-coordinates of the top of the plies
function [local_stresses_bot, local_stresses_top, z_all, ...
    mid_strains_and_curvatures, thermal_loading, ABD] = ...
    get_local_lamina_stresses_planar_ortho(material_props, thetas, ...
    rad_or_deg, thicknesses, mech_loading, delta_T, cte_mat_local)

%Allow for input of just a vector if same material is used
if length(material_props(1,:)) == 1
    material_props = material_props*ones(1,length(thetas));
end
if length(cte_mat_local(1,:)) == 1
    cte_mat_local = cte_mat_local*ones(1,length(thetas));
end

thermal_loading = get_thermal_load(thetas, rad_or_deg,...
    material_props, thicknesses, delta_T, cte_mat_local);
loading = mech_loading + thermal_loading;

%calculating the ABD matrix (essentially laminate stiffness matrix)
ABD = laminate_to_ABD_planar_ortho(thetas, rad_or_deg, thicknesses, ...
    material_props);
% AlphaBetaDelta = inv(ABD);

%calculating mid-plane strains and curvatures (nondimenionalized
%displacement in the 3 translation and 3 rotation axes)
%mid_strains_and_curvatures = AlphaBetaDelta*loading;
mid_strains_and_curvatures = ABD\loading;   %essentially inv(ABD)*loading
[global_strains_bot, global_strains_top] = laminate_to_lamina_strains(...
    mid_strains_and_curvatures, thetas, thicknesses);

global_thermal_strains = get_thermal_strain(delta_T, cte_mat_local, ...
    thetas, rad_or_deg);
global_mech_strains_bot = global_strains_bot - global_thermal_strains;
global_mech_strains_top = global_strains_top - global_thermal_strains;

%Preallocation
num_planar_strain_directions = 3;
global_stresses_bot = zeros(num_planar_strain_directions, length(thetas));
global_stresses_top = zeros(num_planar_strain_directions, length(thetas));
local_stresses_bot = zeros(num_planar_strain_directions, length(thetas));
local_stresses_top = zeros(num_planar_strain_directions, length(thetas));

for k = 1:length(thetas)    %for each ply...
    %Calculate the compliance matrix
    S = get_S_planar_ortho(material_props(:,k));
    Q = inv(S);     %invert to get stiffness matrix
    %Find the stiffness matrix in the gobal coordinate frame
    Q_bar = get_Q_bar_planar(Q, thetas(k), rad_or_deg);
    
    %Solve for the global stresses
    global_stresses_bot(:,k) =  Q_bar*global_mech_strains_bot(:,k);
    global_stresses_top(:,k) =  Q_bar*global_mech_strains_top(:,k);
    
    %Convert stresses from global to local
    local_stresses_bot(:,k) = ...
        global_to_local_stress_planar(global_stresses_bot(:,k), ...
        thetas(k), rad_or_deg);
    local_stresses_top(:,k) = ...
        global_to_local_stress_planar(global_stresses_top(:,k), ...
        thetas(k), rad_or_deg);
    
end

%JUST FOR REPORTING
%find height of the middle of each lamina
total_height = sum(thicknesses);
mid_plane_height = total_height/2;
%finding the z position of the middle, bottom, and top of each ply 
z = cumsum(thicknesses) - mid_plane_height - thicknesses./2;
z_top = z + thicknesses/2;
z_bot = z - thicknesses/2;
z_all = sort([z_bot;z_top]);