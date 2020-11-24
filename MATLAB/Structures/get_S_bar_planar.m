%%This function gets the S_bar compliance matrix to convert the S matrix
%%from the 11,22 coordinate system to the X,Y (planar case)
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param S the planar stiffness matrix in the 11,22,12 coordinate system
%@param theta the angle of the lamina
%@param radOrDeg a string of 'rad' or 'deg' to denote units of theta
%@return S_bar the S planar compliance matrix if the X,Y coordinate system
function S_bar = get_S_bar_planar(S, theta, radOrDeg)

%finding the appropriate transformation matrix
T_sigma = get_T_planar(theta, 'stress', radOrDeg);

%Calculating Q_bar using the transformation matrix
S_bar = T_sigma'*S*T_sigma;