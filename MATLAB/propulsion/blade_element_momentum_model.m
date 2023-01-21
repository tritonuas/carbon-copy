function [outputArg1,outputArg2] = blade_element_momentum_model(numsection,angel,rotation_velocity,plane_velocity,radius)% 

Cl=1.5;
CD=0.5;
chord=0.05;
air_density=1.225;
section_radius=linspace(0,radius,numsection);
for i=1:numsection
V1=sqrt((section_radius(i)*rotation_velocity)^2+plane_velocity^2);
end

dr=radius/numsection;

Lift=Cl*0.5*air_density.*(V1).^2*chord*dr;
Drag=CD*0.5*air_density.*(V1).^2*chord*dr;

Thrust=sum(0.5*air_density*2*(Lift*cos(angel)-Drag*cos(angel))*V1.^2*chord);
Torque=sum(Drag*section_radius);

%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
