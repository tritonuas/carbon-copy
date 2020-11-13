%%Transforms the stresses from the global coordinate system to the fiber
%%coordinate system.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param global_stress the stresses in the global coordinate system
%@param theta the angle between the global and fiber coordinate systems
%@param radOrDeg a string of 'rad' or 'deg' to specify the units of theta
%@return local_stress the stresses in the fiber coordinate system
function local_stress  = global_to_local_stress_planar(global_stress, ...
                                                        theta, radOrDeg)
%get the transformation matrix
T = get_T_planar(theta, "stress", radOrDeg);
%use the transformation matrix to transfer coordinate frame
local_stress = T*global_stress;