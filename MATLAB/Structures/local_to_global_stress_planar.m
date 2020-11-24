%%Transforms the stresses from the local coordinate system to the global
%%coordinate system.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param local_stress the stresses in the fiber coordinate system
%@param theta the angle between the global and local coordinate systems
%@param radOrDeg a string of 'rad' or 'deg' to specify the units of theta
%@return global_stress the stresses in the global coordinate system
function global_stress  = local_to_global_stress_planar(local_stress, ...
                                                        theta, radOrDeg)
%get the transformation matrix
T = get_T_planar(-theta, "stress", radOrDeg);
%use the transformation matrix to transfer coordinate frame
global_stress = T*local_stress;