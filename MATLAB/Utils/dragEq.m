%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program is a MATLAB version of L=W=1/2*rho*velocity^2 *S * Cd
%%Note for emphasis, thrust = drag for level flight

%@param hasDrag boolean to see if a drag was input
%@param drag the lift that the wings are generating
%@param hasRho boolean to see if a rho was input
%@param rho the density of the fluid (probably air, 1.225 at sea level)
%@param hasVel boolean to see if a velocity was input
%@param velocity the velocity of the airflow (or plane in our case)
%@param hasS boolean to see if an S was input
%@param s wing area
%@param hasCd boolean to see if a Cd was input
%@param cd the drag coefficient
%@return lift see param
%@return rho see param
%@return velocity see param
%@return s see param
%@return cd see param
function [drag, rho, velocity, s, cd] = dragEq(hasDrag,drag, hasRho,rho,...
    hasVel,velocity, hasS,s, hasCd,cd)

sum = hasDrag + hasRho + hasVel + hasS + hasCd;
if sum < 4
    disp('you must already have all the variables but one');
    disp('One equation, one unknown');
elseif ~hasDrag
    drag = (1/2*rho*velocity^2)*s*cd;   %q = 1/2*rho*velocity^2
elseif ~hasRho
    rho = drag*2/velocity^2/s/cd;
elseif ~hasVel
    velocity = sqrt(drag/s/cd*2/rho);
elseif ~hasS
    s = drag/(1/2*rho*velocity^2)/cd;
elseif ~hasCd
    cd = drag/(1/2*rho*velocity^2)/s;
end

