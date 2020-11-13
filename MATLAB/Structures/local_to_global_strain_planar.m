%%Transforms the straines from the local coordinate system to the global
%%coordinate system.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param local_strain the straines in the fiber coordinate system
%@param theta the angle between the global and local coordinate systems
%@param radOrDeg a string of 'rad' or 'deg' to specify the units of theta
%@return global_strain the straines in the global coordinate system
function global_strain  = local_to_global_strain_planar(local_strain, ...
                                                        theta, radOrDeg)
%get the transformation matrix
T = get_T_planar(-theta, "strain", radOrDeg);
%use the transformation matrix to transfer coordinate frame
global_strain = T*local_strain;