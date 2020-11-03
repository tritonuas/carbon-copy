%%Authors: Allyson Chen & Kevin Vo
%%Emails:  aechen@ucsd.edu    k1vo@ucsd.edu
%%Calculates the vertical and horizontal tail area and the tail boom length of the plane.
%%This program is a MATLAB version of CHT = (xach-xcg)Sh/cS and CVT = (xach-xcg)Sv/bS

%%Variable Legend
%%x_ach = horizontal tail aero center
%%x_acv = vertical tail aero center
%%x_cg = aircraft center-of-gravity
%%S_h = horizontal tail area
%%S_v = vertical tail area
%%S = (wing) reference area
%%b = wingspan
%%c = mean geo chord
%%C_HT = horizontal tail volume coefficient
%%C_VT = vertical tail volume coefficient

function [tail_area_h, tail_area_v, tail_boom_length] = find_tail_size(wingspan, wingarea, mean_geo_chord, C_HT, C_VT, density, viscosity, velocity)
%% General

format compact
format shortg;
addpath("Utils")

%% Code

minLenTailBoom = 0;
maxLenTailBoom = 1;
numsteps = 1000;
S_v = zeros(1000,1);
S_h = zeros(1000,1);

lenTailBoom = linspace(minLenTailBoom, maxLenTailBoom, numsteps);

for i = 1:length(lenTailBoom)

    hasx_acv = 1;
    hasx_cg = 1;
    hasS_v = 0;
    hasb = 1;
    hasS = 1;
    hasC_VT = 1;
    x_cg = 0;
    lenFuse = 0.25;
    TailBoom_Length = lenTailBoom(i);
    x_acv = TailBoom_Length + lenFuse;
    x_ach = TailBoom_Length + lenFuse;

    
    [wingspan, S_v(i), wingarea, x_acv, x_cg, C_VT] = vert_tail_size_eq(x_acv, x_cg, S_v(i), wingspan, wingarea, C_VT, hasx_acv, hasx_cg, hasS_v, hasb, hasS,...
    hasC_VT);

    hasx_ach = 1;
    hasS_h = 0;
    hasc = 1;
    hasC_HT = 1; 

    [x_ach, x_cg, S_h(i), c, S, C_HT] = horz_tail_size_eq(x_ach, x_cg, S_h(i), mean_geo_chord, wingarea, C_HT, hasx_ach, hasx_cg, hasS_h, hasc, hasS,...
    hasC_HT);

    %if variable 1 do something
    %if variable 2 == -1 do something (doesnt have a value)

    sWing = 1;
    cWing = 99999999999999999999999;
    sFuse = 0;
    lenFuse = 1; 
    sNose = 0;
    lenNose = 1;
    sTail = S_v(i) + S_h(i);
    cTail = -1;
    sTailBoom = 2*pi*x_ach*0.0762; 
    TailBoom_Length = x_ach;

    cd0 = getZeroLiftDrag(density, viscosity, velocity, ...
    sWing,cWing, sFuse,lenFuse, sNose,lenNose, sTail,cTail, ...
    sTailBoom, TailBoom_Length);

    Zero_Lift_Drag_Coeff(i) = cd0;
    
end

min(Zero_Lift_Drag_Coeff);
u = find(Zero_Lift_Drag_Coeff == min(Zero_Lift_Drag_Coeff));
tail_boom_length = lenTailBoom(u)
tail_area_h = S_h(u)
tail_area_v = S_v(u)

