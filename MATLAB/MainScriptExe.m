%%Hi! My name is Andrew, and this is so you can see a push.

%%This is the main script (executable) for the plane design program
%%As of now, this takes in inputs of wing_span, weight, and initial guesses
%%for the sizing of the parts, and outputs the ideal lift coefficient 
%%and wing reference area to maximize cl/cd and therefore maximize range 
%%and minimize drag for cruise conditions.
%%There are currently 2 approaches to this. One is an analyic method that
%%is based off of the method described in Anderson's 155A lecture slides of
%%setting cdi=cd0. This has been modified such that instead of solving for
%%velocity, it takes in a range of velocity as iterative inputs, and
%%solves for the ideal wing area for each velocity.
%%is done using the minDragEq method. The other method is a purely 
%%iterative method for the purpose of checking.
%%Both processes iterate over velocity to try to find an ideal velocity
%%for the weight. I am currently trying to use the analytic method for the
%%solution.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu

%%Sources: See functions for specific sources. 

%%Concerns: a very high sensitivity to changes in cd0, and
%%the input velocity. This is a concern because variance in manufacturing
%%and variances in flight are likely to be great enough to drastically vary
%%the optimal design according to this program.
%%Also, there are discontinuities arising from the transition point for
%%some reason. This can be seen in the iterative solution.

%%Note: If you see a function being called that you aren't familiar with, I
%%probably made the function. If the function is trivial, then it should be
%%saved under the Utils folder. If not, then it should be in the main
%%folder.

%% General
clear all;
clc;

format compact
format shortg;
addpath("Utils")
addpath("Weight");
addpath("Structures");
addpath("tail");
addpath("Airfoil");
addpath("Structures");
addpath("CM_Calc");
%% Constants
g = 9.81;                   %sea level Earth gravity
density = 1.225;            %sea level air density
viscosity = 1.789*10^-5;    %sea level air viscosity

%% Configuration
%subsonic plane
sweepAngle = 0;     %In the future, I can make this a function of velocity
taperRatio = getTaper(sweepAngle);  %see function for source of eq
                                    %taper should equal 0.45 for subsonic
%% Inputs
% I might make this whole program a function in the future
weight = 135;
lift = weight;%want to design this for cruise conditions
cl = .1;
cd0 = .1;
wing_span = 3.65;
chord = 1;
AR = 1;
v = 20;


tail_span_h = .5;
tail_span_v = .5;


%Fuselage Requirements
fuse_surface_area = 1;      %1 for carbon copy
fuse_length = 1.2;   %1.2 for carbon copy

%Tail Boom Requirement
tail_boom_radius = 0.0762/2;    %3 inch diameter converted to radius in m
%Note: This is only cause that's the size of the pole we have. This can be
%made a design variable in the future

%%Preparing for the getZeroLiftDrag function by stating that we do not 
%%have these variables (guesses will be given from within the function)
%%-1 means we do not have the information

nose_surface_area = -1;
nose_length = -1;
tail_area_h = -1;
hs_chord = -1;
tail_area_v = -1;
vs_chord = -1;
tail_boom_surface_area = -1;
tail_boom_length = -1;
wing_area = 5;

%weight inputs
battery = 24.6876;
payload = 71.5719;
num_plies_tailboom = 3;
num_plies_vtail = 2;
num_plies_htail = 2;
num_plies_wing = 2;
num_plies_fuse = 2;
num_spar_wing = 2;
num_spar_vtail = 1.5;
num_sphs_ARtail = 1.5;
spar_width_wing = 0.0127;
spar_width_htail = spar_width_wing;
spar_width_vtail = spar_width_wing;
t_divinycell = 0.003175;
t_tip = 0.12*(2*chord*taperRatio)/(1+taperRatio);
t_root = t_tip/taperRatio;
t_htail_root = (2*(tail_area_h/tail_span_h)*taperRatio)/(1+taperRatio);
t_htail_tip = t_htail_root/taperRatio;
t_vtail_root = (2*(tail_area_v/tail_span_v)*taperRatio)/(1+taperRatio);
t_vtail_tip = t_vtail_root/taperRatio;
t_bulkhead = 0.00635;
density_divinycell = 80;
density_carbon_epoxy = 1600;
density_balsa = 200;
density_plywood = 680;
density_blue_foam = 80;
num_bulkheads = 5;
area_fraction_bulkhead = 0.2;
fuse_height = 0.25;
fuse_width = 0.21;
fuse_length = 1.1;
tailboom_fudge_factor = 1; 
fudge_factor = 1;
wing_fudge_factor  = 1;
htail_fudge_factor = 1;
vtail_fudge_factor = 1;
fuse_fudge_factor  = 1;

%structures inputs
 E1 = 135e9; % Pa
E2 = 10e9; % Pa
G12 = 5e9; % Pa
Nu12 = 0.30;
thickness_per_ply = 0.005;
cte1 = -5e-7;
cte2 = 1.5e-5;
sigma_1T = 1500e6;
sigma_1C = 1200e6;
sigma_2T = 50e6;
sigma_2C = 250e6;
sigma_12 = 70e6;
 mat_props = [E1;E2;G12;Nu12];
cte_vec = [cte1; cte2; 0];
mat_strengths_t = [sigma_1T;sigma_2T;sigma_12];
mat_strengths_c = [sigma_1C;sigma_2C;sigma_12];

fail_crit = "max_stress";
print_output = true;
SF = 2;

Nx = 100;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
delta_T = 0;
t_airfoil = .1;
thetas = 90;
rad_or_deg = "deg";
thickness_per_ply = 0.0003;
thicknesses = ones(length(thetas),1)*thickness_per_ply;
Rm = -lift*(wing_span/4); % reaction moment at root
% t_airfoil = airfoil thickness
Fx = -Rm/t_airfoil; %Force
l = 0.8*chord - 0.2*chord; % length of wing box
Nx = Fx/l; % in-plane stress
%get induced drag so we can later get cd
k = getK(AR, taperRatio, sweepAngle);
cdi = cl^2*k;
%get cd so we can later get cl/cd
cd = cd0 + cdi;
drag = 1/2*density*v^2*cd*wing_area;
Rm_y = drag*(wing_span/4); % reaction moment at root in the y direction
wing_box_distance  = 0.8*chord - 0.2*chord; % distance between two spars 
Nxy = Rm_y/(wing_box_distance*(wing_span/2));
mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];

%tail area inputs
sTail = tail_area_h + tail_area_v;
tail_boom_surface_area = tail_boom_length*2*pi*tail_boom_radius;
hs_AR=4;
vs_AR=1;
tail_span_h=sqrt(tail_area_h*hs_AR);
tail_span_v=sqrt(tail_area_v*vs_AR);
tail_chord_h=tail_area_h/tail_span_h;
tail_chord_v=tail_area_v/tail_span_v;

%weight vec
battery = 1;
wing_weight = 1;
fuse_weight = 1;
vs_weight = 1;
hs_weight = 1;
tail_boom_weight = 1;
weight_vec = [battery,wing_weight,fuse_weight,vs_weight,hs_weight,tail_boom_weight];

%position vec
position_vec = [0.05,0,0;
    (chord/2)+(.2*fuse_length),0,0.05;
    fuse_length/2,0,0;
    (tail_chord_v/2)+tail_boom_length+fuse_length,0,tail_span_v/2;
    (tail_chord_h/2)+tail_boom_length+fuse_length,0,0;
    (fuse_length+tail_boom_length/2),0,0];

%% Desired Parameters/Constraints

%%CONSTRAINT: MAX CL CRUISE
maxClCruise = 1.8;  %see iterative solution for implementation
maxLoadFactorStall = 1.8;     %input max load factor for takeoff/landing
maxLoadFactorTurns = 1.667; %input max load factor for turns

% %%CONSTRIANT: ROOT CHORD LENGTH (STRUCTURAL)
% %This is for structural reasons
% rootChordReq = 0.3048;  %carbon copy is 0.3048 this is 12 inches
%   %That value is based on intuition/history and should be iterated on
% tipChordReq = rootChordReq*taperRatio;
% avgChordReq = (rootChordReq+tipChordReq)/2;
% sMinStruct = avgChordReq*wing_span;
% if sMinStruct > sMin
%     sMin = sMinStruct;
%     S = linspace(sMin, sMax, numSteps);
%     SMaster = S;
%     AR = wing_span^2./S;  %Getting ARs for induced drag calculations
%     ARMaster = AR;
% end

%%CONSTRAINT: STALL SPEED
% stallSpeed = 30;
% % L = 1/2*density*stallSpeed^2*S*clMax    nMax = clMax/cl    
% %stallSpeed = velocity/sqrt(nMax)   nMax = (velocity/stallSpeed)^2
% nStall = (velocityMax/stallSpeed)^2;
% if nStall > maxLoadFactorStall
%    nStall = maxLoadFactorStall;
%    velocityMax = stallSpeed*sqrt(nStall);
%    velocity = linspace(velocityMin, velocityMax, numSteps);
% end

%%CONSTRAINT: TURNING RADIUS
% radius = 99;    %want this radius max 33.98m
% %radius = v^2/(g*sqrt(n^2-1))
% n = sqrt((velocityMax^2/radius/g)^2 + 1);%load factor req for desired turns
% if n > maxLoadFactorTurns
%     n = maxLoadFactorTurns;
%     velocityMax = sqrt(radius*g*sqrt(n^2-1));  %velocity req for turns
%     velocity = linspace(velocityMin, velocityMax, numSteps);
% end


%%Add climb rate constraint TODO
%%PROBLEM: This requires a P/W, but it's hard to have information on thrust
% climbAngleReq = 10;     %in degrees
% G = sind(climbAngleReq);
% %Continue to work on this
% 
% %%CONSTRAINT: MISSION TIME
% %Mission should be completed within 20 minutes.
% timeLimit = 9999;
% SFTime = 1.05;   %safety factor
% timeLimit = timeLimit/SFTime;
% timeLimit = timeLimit*60;   %time limit in seconds
% %Consult with Garrett to determine whether we should take the whole time or
% %cut out time for making second passes.
% totalTravelDist = 17000;
% additionalWaypointDist = 10000;
% totalTravelDist = totalTravelDist + additionalWaypointDist;
% velocityReq = totalTravelDist/(timeLimit);
% if velocityReq > velocityMin
%     velocityMin = velocityReq;
%     velocity = linspace(velocityMin, velocityMax, numSteps);
% end
% 
% %%CONTRAINT: TAIL SIZING
% %This is based off of Anderson's MAE 155A slides(see stability, last slide)
C_HT = 0.8;       %Horizontal tail volume coefficient
C_VT = 0.04;     %Vertical tail volume coefficient

radius = 30;
%% Gradient based optimization

%TODO still need to finish validating and cleaning up stuff that is not
%needed, and finish making the output print statements
%TODO think about this from a controller standpoint, and how to add Kd and
%Ki terms if I need to.
%declaring that I have these variables for use of the lift equation
hasLift = 1;
hasDensity = 1;
hasVel = 1;

hasS = 1;       %declaring that I have S and want to calculate for Cl
hasCl = 0;      %This is in preperation for the lift equation

x = [wing_area];
Kp_S = 1;
iterNum = 0;
derivative_cl_over_cd = 1;   %to run the loop
delta_S = 0.01; %first step
while abs(derivative_cl_over_cd) > 0.001
    iterNum = iterNum + 1;
    wing_area = x(1);
    [AR,chord] = geometric_outputs(wing_area,wing_span);
    
    weight = 135; % initialize weight       
    v_guess = 20;
    v = 10;
    while abs(v_guess - v) > 0.1
        v_guess = v;
        
        wing_surface_area = wing_area*2;
%         saTail = sTail*2;
        hs_surface_area = tail_area_h*2;
        vs_surface_area = tail_area_v*2;
        %get the zero-lift drag coeff for this S and velocity
        cd0 = getZeroLiftDrag(density, viscosity, v, ...
                   wing_area, wing_surface_area, chord, fuse_surface_area,fuse_length, nose_surface_area,nose_length,...
                   hs_surface_area,hs_chord, vs_surface_area,vs_chord, tail_boom_surface_area,tail_boom_length);

        % calculate velocity and lift coefficient
        v = -1;
        [wing_area, v, cl] = minDragEq(v,...
          density, wing_area, cd0, wing_span, weight, taperRatio, sweepAngle);
        wing_surface_area = wing_area*2;
      
        % calculate tail areas and length and tail boom length (and therefore area)
        [tail_area_h, tail_area_v, tail_boom_length] = find_tail_size(wing_span, ...
        wing_area, wing_surface_area, chord, C_HT, C_VT, density, viscosity, v(end));
        
        
        [stresses_bot, stresses_top, z_all, ...
        mid_strains_and_curvatures, thermal_loading, ABD] = ...
        get_local_lamina_stresses_planar_ortho(mat_props, thetas, ...
        rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec);

        [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
        = report_ply_margins(stresses_bot, stresses_top, z_all, ...
        fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output)
        


        [weight,wing_weight,hs_weight,vs_weight,fuse_weight,tail_boom_weight] = compute_weight_analytic(battery, payload, wing_area, num_spar_wing,  spar_width_wing, density_balsa, t_divinycell, density_divinycell,...
            num_plies_wing, density_carbon_epoxy, t_tip, t_root, wing_span, fudge_factor,...
            num_sphs_ARtail, spar_width_htail, num_plies_htail, tail_area_h, density_blue_foam, t_htail_root, t_htail_tip, tail_span_h,...
            num_spar_vtail, spar_width_vtail, num_plies_vtail, tail_area_v, t_vtail_root, t_vtail_tip, tail_span_v,...
            num_plies_fuse, num_bulkheads, t_bulkhead, area_fraction_bulkhead, density_plywood, fuse_height, fuse_width, fuse_length,...
            num_plies_tailboom, tail_boom_length, wing_fudge_factor, htail_fudge_factor, vtail_fudge_factor, fuse_fudge_factor, tailboom_fudge_factor);

    end
  
    weight_vec = [battery,wing_weight,fuse_weight,vs_weight,hs_weight,tail_boom_weight];
    position_vec = [0.05,0,0;
        (chord/2)+(.2*fuse_length),0,0.05;
        fuse_length/2,0,0;
        (tail_chord_v/2)+tail_boom_length+fuse_length,0,tail_span_v/2;
        (tail_chord_h/2)+tail_boom_length+fuse_length,0,0;
        (fuse_length+tail_boom_length/2),0,0];
    [cg] = calc_cg(weight_vec, position_vec);
    
    mass_vec = [weight_vec/9.81]';
    height_vec = [.035,.13,.0475;
        .12*chord,chord,wing_span;
        fuse_height,fuse_length,fuse_width;
        tail_span_v,tail_chord_v,.12*tail_chord_v;
        .12*tail_chord_h,tail_chord_h,tail_span_h;
        .0254,tail_boom_length,.0254];
    width_vec = height_vec(:,[3, 1, 2]);

    
    component_moi = calc_component_moi(mass_vec,height_vec,width_vec);
    
    quarter_chord = [(chord/2)+(.2*fuse_length),0,0.05];
    central_moment = calc_central_moment(component_moi,weight_vec,position_vec,quarter_chord);
    
    q = 1/2*density*v^2;
    
    %Placing a maximum lift coefficient due to airfoil constraints
    if cl > maxClCruise
         cl = maxClCruise;
         hasS = 0;
         hasCl = 1;
         %Overriding the S to meet the new constraint with a fixed cl
         [lift,density,v,wing_area, cl] = liftEq(hasLift,lift,...
         hasDensity,density,hasVel,v,hasS,wing_area,hasCl,cl);
         %Resetting design variables for future iterations
         hasS = 1;
         hasCl = 0;
         chord = wing_area/wing_span;   %updating chord for new S
         AR = wing_span^2/wing_area;    %updating asepct ratio for new S
    end

    %get induced drag so we can later get cd
    k = getK(AR, taperRatio, sweepAngle);
    cdi = cl^2*k;
    %get cd so we can later get cl/cd
    cd = cd0 + cdi;
    %get cl/cd for performance comparison
    clOverCd = cl/cd;
    
    break

end

% figure(1);
% plot(2:iterNum, clOverCd(2:end), '-bo');
% xlabel("Iteration Number"); ylabel("CL/CD"); 
% title("Derivative-Based Iteration: Cl/CD vs. Iteration Number");
% figure(2);
% plot(2:iterNum, S(2:end), '-bo');
% xlabel("Iteration Number"); ylabel("Wing Area (m^2)"); 
% title("Derivative-Based Iteration: Wing Area vs. Iteration Number");

%update load factor
load_factor = sqrt((v^2/radius/g)^2 + 1);
%update stall speed
stall_speed = v(end)/sqrt(maxLoadFactorStall);

wing_loading = lift/wing_area;

disp("----------------------------");
disp("DBI_ANA SOLUTION");
disp("----------------------------");
disp("Max cl/cd: " + clOverCd(end));
disp("Velocity: " + v(end));
disp("cl: " + cl(end));
disp("Wing Area: " + wing_area);
disp("Wing loading: " + wing_loading);
disp("Chord: " + chord);
disp("Aspect ratio: " + AR);
disp("Load factor on turns: " + load_factor);
disp("Horizontal Tail Area: " + tail_area_h); 
disp("Vertical Tail Area: " + tail_area_v); 
disp("Tail Boom Length: " + tail_boom_length); 
if load_factor > maxLoadFactorTurns + 0.001
    disp("Invalid solution! Load factor constraint not satisfied.");
end
disp("Stall speed: " + stall_speed);
% 





