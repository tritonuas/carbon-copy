%%Takes in inputs of material and laminate properties and calculates the
%%local stresses at the top and bottom of each ply (planar).
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%%@param material_props a vector defined as [E11 E22 G12 Nu12]
%%@param thetas a vector of the angles of each lamina (in order)
%%@param rad_or_deg a string of 'rad' or 'deg' specifying units for theta
%%@param thicknesses a vector of the thickness of each ply
%%@param loading the column vector of loading [Nx;Ny;Nxy;Mx;My;Mxy]
function [local_stresses_bot, local_stresses_top] = ...
    get_local_lamina_stresses_planar_ortho(material_props, thetas, ...
    rad_or_deg, thicknesses, loading)

%unpackaging material properties
E1 = material_props(1);
E2 = material_props(2);
G12 = material_props(3);
Nu12 = material_props(4);

%calculating the ABD matrix (essentially laminate compliance matrix)
ABD = laminate_to_ABD_planar_ortho(thetas, rad_or_deg, thicknesses, E1, E2, ...
    G12, Nu12);
%AlphaBetaDelta = inv(ABD)

%calculating mid-plane strains and curvatures (nondimenionalized
%displacement in the 3 translation and 3 rotation axes)
%mid_strains_and_curvatures = AlphaBetaDelta*loading;
mid_strains_and_curvatures = ABD\loading;   %essentially inv(ABD)*loading
[local_strains_bot, local_strains_top] = laminate_to_lamina_strains(...
    mid_strains_and_curvatures, thetas, rad_or_deg, thicknesses);

S = get_S_planar_ortho(E1, E2, G12, Nu12);
%Q = inv(S);    MATLAB likes A\b better than inv(A)*b

%Preallocation
num_planar_strain_directions = 3;
local_stresses_bot = zeros(num_planar_strain_directions, length(thetas));
local_stresses_top = zeros(num_planar_strain_directions, length(thetas));

for k = 1:length(thetas)
    %local_stresses_bot(:,k) =  Q*local_strains_bot(:,k);
    local_stresses_bot(:,k) =  S\local_strains_bot(:,k);
    local_stresses_top(:,k) =  S\local_strains_top(:,k);
end
