%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This prgram is for finding the Reynold's Number
%%This is the most profound program ever.

%Source for equation: Professor Anderson (MAE2/MAE155A lecture slides),
%also google and MAE101A, it's really all over the place.
%Assumptions: none
%Other notes on equation(s): it's nondimensional
%Concerns: none
%More research recommended?: no

%@param density the density of the fluid
%@param velocity the velocity of the fluid
%@param charLength the characteristic length of the body
%@param viscosity the viscosity of the fluid    NOT KINEMATIC VISCOSITY
function [Re] = getRe(density,velocity,charLength,viscosity)
%For aicraft specifically: can look at the system as the aircraft
%stationary, and the air flowing with velocity
%density is the air density (varies with temperature (and therefore alt))
%velocity is the velocity of the plane
%charLength is the chord length of the wing (for wings)
%viscosity is the viscosity of the air
Re = density*velocity*charLength/viscosity;