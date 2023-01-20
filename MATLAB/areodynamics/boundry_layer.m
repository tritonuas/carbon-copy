function Thickness = boundry_layer(kinematic_viscosity,freestream_velocity,position_on_wings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Thickness=5*sqrt(position_on_wings*kinematic_viscosity/freestream_velocity);
end