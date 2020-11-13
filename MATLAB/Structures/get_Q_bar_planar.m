%%This function gets the Q_bar stiffness matrix to convert the Q matrix
%%from the 11,22 coordinate system to the X,Y (planar case)
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param Q the planar stiffness matrix in the 11,22,12 coordinate system
%@param theta the angle of the lamina
%@param radOrDeg a string of 'rad' or 'deg' to denote units of theta
%@return Q_bar the Q planar stiffness matrix if the X,Y coordinate system
function Q_bar = get_Q_bar_planar(Q, theta, radOrDeg)

%finding the appropriate transformation matrix
T_epsilon = get_T_planar(theta, 'strain', radOrDeg);

%Calculating Q_bar using the transformation matrix
Q_bar = T_epsilon'*Q*T_epsilon;