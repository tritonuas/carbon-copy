%%Transforms the strains from the global coordinate system to the local
%%coordinate system.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param global_strain the straines in the global coordinate system
%@param theta the angle between the global and fiber coordinate systems
%@param radOrDeg a string of 'rad' or 'deg' to specify the units of theta
%@return local_strain the straines in the fiber coordinate system
function local_strain  = global_to_local_strain_planar(global_strain, ...
                                                        theta, radOrDeg)
%get the transformation matrix
T = get_T_planar(theta, "strain", radOrDeg);
%use the transformation matrix to transfer coordinate frame
local_strain = T*global_strain;