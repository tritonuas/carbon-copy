%%Finds the local lamina strains based on the laminate strains
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%@param mid_strains_and_curvatures the midplane strains and curvatures
%%@param thetas a vector of the angles of each lamina (in order)
%%@param rad_or_deg a string of 'rad' or 'deg' specifying units for theta
%%@param thicknesses a vector of the thickness of each ply
%%@return local_strains_bot the local strain at the bottom of each ply
%%@return local_strains_top the local strain at the top of each ply
function [local_strains_bot, local_strains_top] = ...
    laminate_to_lamina_strains(mid_strains_and_curvatures, ...
                                thetas, rad_or_deg, thicknesses)

%find height of the middle of each lamina
total_height = sum(thicknesses);
mid_plane_height = total_height/2;

%Preallocation
num_strain_axes = 3;
strains_bot = zeros(num_strain_axes, length(thetas));
strains_top = zeros(num_strain_axes, length(thetas)); 
local_strains_bot = zeros(num_strain_axes, length(thetas));
local_strains_top = zeros(num_strain_axes, length(thetas));

for k = 1:length(thetas)    %for each ply...
    %finding the z position of the middle, bottom, and top of each ply 
    z = cumsum(thicknesses) - mid_plane_height - thicknesses./2;
    z_top = z(k) + thicknesses(k)/2;
    z_bot = z(k) - thicknesses(k)/2;
    
    %finding the strain at the top and bottom of each ply
    strains_bot(:,k) = mid_strains_and_curvatures(1:3) ...
        + z_bot*mid_strains_and_curvatures(4:6);
    strains_top(:,k) = mid_strains_and_curvatures(1:3) ...
        + z_top*mid_strains_and_curvatures(4:6);
    
    %converting the strains from a global to local coordinate system
    local_strains_bot(:,k) = global_to_local_strain_planar(...
        strains_bot(:,k), thetas(k), rad_or_deg);
    local_strains_top(:,k) = global_to_local_strain_planar(...
        strains_top(:,k), thetas(k), rad_or_deg);
end