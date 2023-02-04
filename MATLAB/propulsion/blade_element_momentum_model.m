%function [Thrust,Torque] = blade_element_momentum_model(numsection,angle,rotation_velocity,plane_velocity,radius)% 
clear all;




    numsection = 30;
    angle=linspace(45,20,numsection);
    rotation_velocity=340.5;
 
    plane_velocity=16;
    radius=0.33;
    num_blades=2;
    
    
    
    
    Cl=0.75;
    CD=0.5;
    chord=0.05; %% Add chord distribution next
    air_density=1.225;
    section_radius=linspace(0,radius,numsection);
    
    for i=1:numsection
    V(i)=sqrt((section_radius(i)*rotation_velocity)^2+plane_velocity^2);
    phi(i)=atan2d(plane_velocity, rotation_velocity*section_radius(i)); %Fixed phi calculations.. it was calculating alpha before
    alpha(i)=(angle(i)-phi(i))*pi/180; %Gets the Alpha
    Cl(i)=(2*pi)*alpha(i)*0.8;
    end
    
    dr=radius/numsection;
    
    Lift=Cl(i).*0.5*air_density.*(V).^2*chord*dr;
    Drag=CD*0.5*air_density.*(V).^2*chord*dr;
    
    Thrust=num_blades.*sum(Lift.*cosd(phi)-Drag.*sind(phi));
    Torque=sum(Drag.*section_radius);

    
%   Detailed explanation goes here

