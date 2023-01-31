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
    chord=0.05;
    air_density=1.225;
    section_radius=linspace(0,radius,numsection);
    
    for i=2:numsection
    V(i)=sqrt((section_radius(i)*rotation_velocity)^2+plane_velocity^2);
    phi(i)=angle(i)-atand(plane_velocity/(section_radius(i)*rotation_velocity));
    end
    
    dr=radius/numsection;
    
    Lift=Cl*0.5*air_density.*(V).^2*chord*dr;
    Drag=CD*0.5*air_density.*(V).^2*chord*dr;
    
    Thrust=num_blades.*sum(Lift.*cosd(phi)-Drag.*sind(phi));
    Torque=sum(Drag.*section_radius);

    
%   Detailed explanation goes here

