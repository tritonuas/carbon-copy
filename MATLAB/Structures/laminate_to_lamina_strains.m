%%Finds the local lamina strains based on the laminate strains
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%@param mid_strains_and_curvatures the midplane strains and curvatures
%%@param thetas a vector of the angles of each lamina (in order)
%%@param thicknesses a vector of the thickness of each ply
%%@return strains_bot the global strain at the bottom of each ply
%%@return strains_top the global strain at the top of each ply
function [strains_bot, strains_top] = ...
    laminate_to_lamina_strains(mid_strains_and_curvatures, ...
                                thetas, thicknesses)

%find height of the middle of each lamina
total_height = sum(thicknesses);
mid_plane_height = total_height/2;
%finding the z position of the middle, bottom, and top of each ply 
z = cumsum(thicknesses) - mid_plane_height - thicknesses./2;
z_top = z + thicknesses/2;
z_bot = z - thicknesses/2;

%Preallocation
num_strain_axes = 3;
strains_bot = zeros(num_strain_axes, length(thetas));
strains_top = zeros(num_strain_axes, length(thetas)); 

for k = 1:length(thetas)    %for each ply...
    %finding the strain at the top and bottom of each ply
    strains_bot(:,k) = mid_strains_and_curvatures(1:3) ...
        + z_bot(k)*mid_strains_and_curvatures(4:6);
    strains_top(:,k) = mid_strains_and_curvatures(1:3) ...
        + z_top(k)*mid_strains_and_curvatures(4:6);
end