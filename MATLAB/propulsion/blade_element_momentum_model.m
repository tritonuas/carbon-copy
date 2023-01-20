function [outputArg1,outputArg2] = blade_element_momentum_model(angel,rotation_velocity,plane_velocity,radius)% 

Cl=1.5;
CD=0.5;
chord=0.05;
air_density=1.225;
for r=0:radius
V1=sqrt((r*rotation_velocity)^2+plane_velocity^2);
end

Lift=Cl*0.5*air_density.*(V1).^2*chord;
Drag=CD*0.5*air_density.*(V1).^2*chord;

Thrust=0.5*air_density*2*(Lift*cos(angel)-Drag*cos(angel))*V1.^2*chord;
drag=CD*air_density*0.5*V1*r;
lift=Cl*0.5*air_density*V1*r;

%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
