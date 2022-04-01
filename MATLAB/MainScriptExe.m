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
clear;
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
addpath("optimizers");
%% Constants
g = 9.81;                   %sea level Earth gravity
density = 1.225;            %sea level air density
viscosity = 1.789*10^-5;    %sea level air viscosity

%% Configuration
%subsonic plane
sweepAngle = 0;     %In the future, I can make this a function of velocity
taper_ratio = getTaper(sweepAngle);  %see function for source of eq
                                    %taper should equal 0.45 for subsonic
%% Inputs
% I might make this whole program a function in the future

alpha = .0001;
h = 1e-5;

%operating conditions
v = 20;
weight = 135;
lift = weight;%want to design this for cruise conditions
battery = 24.6876;
payload = 71.5719;

%Tail Boom Requirement    %3 inch diameter converted to radius in m
%Note: This is only cause that's the size of the pole we have. This can be
%made a design variable in the future

%geometric inputs

cl = .1;
cd0 = .1;
wing_span = 3.65;
wing_chord = 10;
wing_area = wing_span*wing_chord;
AR = wing_chord^2/wing_area;
wing_num_plies = 2;
wing_num_spar = 2;
wing_spar_width = 0.0127;
wing_tip_thickness = 0.12*(2*wing_chord*taper_ratio)/(1+taper_ratio);
wing_root_thickness = wing_tip_thickness/taper_ratio;

nose_surface_area = .05;
nose_length = .1;

hs_span = .5;
hs_chord = .15;
hs_area = hs_span*hs_chord;
hs_AR=hs_chord^2/hs_area;
hs_num_plies = 2;
hs_num_spar = 1.5; %num_sphs_ARtail
hs_spar_width = wing_spar_width/2;
hs_root_thickness = (2*(hs_area/hs_span)*taper_ratio)/(1+taper_ratio);
hs_tip_thickness = hs_root_thickness/taper_ratio;

vs_span = .5;
vs_chord = .15;
vs_area = vs_span*vs_chord;
vs_AR=vs_chord^2/vs_area;
vs_num_plies = 2;
vs_num_spar = 1.5;
vs_spar_width = wing_spar_width/2;
vs_root_thickness = (2*(vs_area/vs_span)*taper_ratio)/(1+taper_ratio);
vs_tip_thickness = vs_root_thickness/taper_ratio;

tail_boom_length = .9;
tail_boom_radius = 0.0762/2;
tail_boom_surface_area = tail_boom_length*2*pi*tail_boom_radius;
tail_boom_num_plies = 3;

fuse_surface_area = 1;      %1 for carbon copy
fuse_length = 1.2;
fuse_height = 0.25;
fuse_width = 0.21;
fuse_num_plies = 2;

bulkhead_thickness = 0.00635;
bulkheads_num = 5;
bulkhead_area_fraction = 0.2;

%weight inputs 

fudge_factor = 1;
fuse_fudge_factor  = 1;
wing_fudge_factor  = 1;
hs_fudge_factor = 1;
vs_fudge_factor = 1;
tailboom_fudge_factor = 1;

%material properties
divinycell_density = 80;
divinycell_thickness = 0.003175;
carbon_epoxy_density = 1600;
balsa_density = 200;
plywood_density = 680;
blue_foam_density = 80;

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
print_output = false;
SF = 2;

Nx = 100;
Ny = 0;
Nxy = 0;
Mx = 0;
My = 0;
Mxy = 0;
delta_T = 0;
airfoil_thickness = .1;
thetas = 90;
rad_or_deg = "deg";
thicknesses = ones(length(thetas),1)*thickness_per_ply;
Rm = -lift*(wing_span/4); % reaction moment at root
% t_airfoil = airfoil thickness
Fx = -Rm/airfoil_thickness; %Force
l = 0.8*wing_chord - 0.2*wing_chord; % length of wing box
Nx = Fx/l; % in-plane stress
%get induced drag so we can later get cd
k = getK(AR, taper_ratio, sweepAngle);
cdi = cl^2*k;
%get cd so we can later get cl/cd
cd = cd0 + cdi;
drag = 1/2*density*v^2*cd*wing_area;
Rm_y = drag*(wing_span/4); % reaction moment at root in the y direction
wing_box_distance  = 0.8*wing_chord - 0.2*wing_chord; % distance between two spars 
Nxy = Rm_y/(wing_box_distance*(wing_span/2));
mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];

%tail area inputs
tail_area_total = hs_area + vs_area;
hs_span=sqrt(hs_area*hs_AR);
vs_span=sqrt(vs_area*vs_AR);

%weight vec
wing_weight = 1;
fuse_weight = 1;
vs_weight = 1;
hs_weight = 1;
tail_boom_weight = 1;
weight_vec = [battery,wing_weight,fuse_weight,vs_weight,hs_weight,tail_boom_weight];

%position vec
position_vec = [0.05,0,0;
    (wing_chord/2)+(.2*fuse_length),0,0.05;
    fuse_length/2,0,0;
    (vs_chord/2)+tail_boom_length+fuse_length,0,vs_span/2;
    (hs_chord/2)+tail_boom_length+fuse_length,0,0;
    (fuse_length+tail_boom_length/2),0,0];
    
%mass vector  
mass_vec = (weight_vec/9.81)';
height_vec = [.035,.13,.0475;
    .12*wing_chord,wing_chord,wing_span;
    fuse_height,fuse_length,fuse_width;
    vs_span,vs_chord,.12*vs_chord;
    .12*hs_chord,hs_chord,hs_span;
    .0254,tail_boom_length,.0254];
    width_vec = height_vec(:,[3, 1, 2]);

%% Desired Parameters/Constraints

%%CONSTRAINT: MAX CL CRUISE
maxClCruise = 1.8;  %see iterative solution for implementation
maxLoadFactorStall = 1.8;     %input max load factor for takeoff/landing
maxLoadFactorTurns = 1.667; %input max load factor for turns

% %%CONSTRIANT: ROOT wing_chord LENGTH (STRUCTURAL)
% %This is for structural reasons
% rootChordReq = 0.3048;  %carbon copy is 0.3048 this is 12 inches
%   %That value is based on intuition/history and should be iterated on
% tipChordReq = rootChordReq*taper_ratio;
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

x = [wing_chord];
iter_num = 0;
gradient = ones(length(x),1);   %to run the loop
while norm(gradient) > 0.001
    iter_num = iter_num + 1;
    for i = 1:(length(x)+1)
    if i ~= 1
    x_step = zeros(length(x),1);
    x_step(i-1) = h;
    x = x + x_step;
    end
    wing_chord = x(1);
        
    [AR,wing_area,wing_tip_thickness,wing_root_thickness] = geometric_outputs(taper_ratio,wing_chord,wing_span);
    
    [weight,wing_weight,hs_weight,vs_weight,fuse_weight,tail_boom_weight] = compute_weight_analytic(battery, payload, wing_area, wing_num_spar,  wing_spar_width, balsa_density, divinycell_thickness, divinycell_density,...
            wing_num_plies, carbon_epoxy_density, wing_tip_thickness, wing_root_thickness, wing_span, fudge_factor,...
            hs_num_spar, hs_spar_width, hs_num_plies, hs_area, blue_foam_density, hs_root_thickness, hs_tip_thickness, hs_span,...
            vs_num_spar, vs_spar_width, vs_num_plies, vs_area, vs_root_thickness, vs_tip_thickness, vs_span,...
            fuse_num_plies, bulkheads_num, bulkhead_thickness, bulkhead_area_fraction, plywood_density, fuse_height, fuse_width, fuse_length,...
            tail_boom_num_plies, tail_boom_length, wing_fudge_factor, hs_fudge_factor, vs_fudge_factor, fuse_fudge_factor, tailboom_fudge_factor);
           
    v_guess = 20;
    v = 10;
    while abs(v_guess - v) > 0.1
        v_guess = v;
        
        wing_surface_area = wing_area*2;
%         saTail = tail_area_total*2;
        hs_surface_area = hs_area*2;
        vs_surface_area = vs_area*2;
        %get the zero-lift drag coeff for this S and velocity
        cd0 = getZeroLiftDrag(density, viscosity, v, ...
                   wing_area, wing_surface_area, wing_chord, fuse_surface_area,fuse_length, nose_surface_area,nose_length,...
                   hs_surface_area,hs_chord, vs_surface_area,vs_chord, tail_boom_surface_area,tail_boom_length);

        % calculate velocity and lift coefficient
        v = -1;
        [wing_area, v, cl] = minDragEq(v,...
          density, wing_area, cd0, wing_span, weight, taper_ratio, sweepAngle);
        wing_surface_area = wing_area*2;
      
        % calculate tail areas and length and tail boom length (and therefore area)
%         [hs_area, vs_area, tail_boom_length] = find_tail_size(wing_span, ...
%         wing_area, wing_surface_area, wing_chord, C_HT, C_VT, density, viscosity, v(end));
        

    end

    
    component_moi = calc_component_moi(mass_vec,height_vec,width_vec);
    
    quarter_chord = [(wing_chord/2)+(.2*fuse_length),0,0.05];
    central_moment = calc_central_moment(component_moi,weight_vec,position_vec,quarter_chord);
    
    %get induced drag so we can later get cd
    k = getK(AR, taper_ratio, sweepAngle);
    cdi = cl^2*k;
    %get cd so we can later get cl/cd
    cd = cd0 + cdi;
    drag = 1/2*density*v^2*cd*wing_area;
    %get cl/cd for performance comparison
    clOverCd = cl/cd;
    
    Nx = 100;
    Ny = 0;
    Nxy = 0;
    Mx = 0;
    My = 0;
    Mxy = 0;
    delta_T = 0;
    airfoil_thickness = .1;
    thetas = 90;
    rad_or_deg = "deg";
    thicknesses = ones(length(thetas),1)*thickness_per_ply;
    Rm = -lift*(wing_span/4); % reaction moment at root
    % t_airfoil = airfoil thickness
    Fx = -Rm/airfoil_thickness; %Force
    l = 0.8*wing_chord - 0.2*wing_chord; % length of wing box
    Nx = Fx/l; % in-plane stress
    Rm_y = drag*(wing_span/4); % reaction moment at root in the y direction
    wing_box_distance  = 0.8*wing_chord - 0.2*wing_chord; % distance between two spars 
    Nxy = Rm_y/(wing_box_distance*(wing_span/2));
    mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];
    
    [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
     thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
     fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
     rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
     mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output);
 
    min_ms = min(MS(:));
    
    ms_constraint = -min_ms;
    ms_constraint_rho = 10;
    
    constraint_vec = [ms_constraint]; %make sure to make vertical
    constraint_rho_vec = [ms_constraint_rho];
    active_constraint_vec = constraint_vec(constraint_vec > 0);
    active_constraint_rho_vec = constraint_rho_vec(constraint_vec > 0);
    penalty_scaling_factors = diag(active_constraint_rho_vec);
    objective_penalty = 1/2 * active_constraint_vec' * penalty_scaling_factors * active_constraint_vec;
    if isempty(objective_penalty)
        objective_penalty = 0;
    end
    
    if i == 1
        f_x = drag + objective_penalty;
    else
        f_x_plus_h = drag + objective_penalty;
        gradient(i-1) = (f_x_plus_h - f_x)/h;
        x = x - x_step;
    end
    end     
    x_history(:,iter_num) = x;
    objective_history(iter_num) = f_x;
    gradient_norm(iter_num) = gradient;
    x = gradient_descent_optimizer(gradient,alpha,x)
end

figure(1);
plot((1:iter_num),x_history(1,:))
xlabel("Iteration Number"); ylabel("Wing Chord"); 
title("Gradient-Based Iteration: Wing Chord vs. Iteration Number");
figure(2);
plot((1:iter_num),gradient_norm)
xlabel("Iteration Number"); ylabel("Gradient Norm"); 
title("Gradient-Based Iteration: Gradient Norm vs. Iteration Number");
figure(3);
plot((1:iter_num),objective_history)
xlabel("Iteration Number"); ylabel("Objective"); 
title("Gradient-Based Iteration: Objective vs. Iteration Number");

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
disp("Max cl/cd: " + clOverCd);
disp("Velocity: " + v);
disp("cl: " + cl);
disp("Wing Area: " + wing_area);
disp("Wing loading: " + wing_loading);
disp("wing_chord: " + wing_chord);
disp("Aspect ratio: " + AR);
disp("Load factor on turns: " + load_factor);
disp("Horizontal Tail Area: " + hs_area); 
disp("Vertical Tail Area: " + vs_area); 
disp("Tail Boom Length: " + tail_boom_length); 
if load_factor > maxLoadFactorTurns + 0.001
    disp("Invalid solution! Load factor constraint not satisfied.");
end
disp("Stall speed: " + stall_speed);

%% debugging brute force method
% wing_chords = linspace(11,.015,1000);
% iter_num = 0;
% gradient_norm = 0;
% objective_history = 0;
% for j = 1:length(wing_chords)
%     x(1) = wing_chords(j);
%     iter_num = iter_num + 1;
%     for i = 1:(length(x)+1)
%     if i ~= 1
%     x_step = zeros(length(x),1);
%     x_step(i-1) = h;
%     x = x + x_step;
%     end
%     wing_chord = x(1);
%         
%     [AR,wing_area,wing_tip_thickness,wing_root_thickness] = geometric_outputs(taper_ratio,wing_chord,wing_span);
%     
%     [weight,wing_weight,hs_weight,vs_weight,fuse_weight,tail_boom_weight] = compute_weight_analytic(battery, payload, wing_area, wing_num_spar,  wing_spar_width, balsa_density, divinycell_thickness, divinycell_density,...
%             wing_num_plies, carbon_epoxy_density, wing_tip_thickness, wing_root_thickness, wing_span, fudge_factor,...
%             hs_num_spar, hs_spar_width, hs_num_plies, hs_area, blue_foam_density, hs_root_thickness, hs_tip_thickness, hs_span,...
%             vs_num_spar, vs_spar_width, vs_num_plies, vs_area, vs_root_thickness, vs_tip_thickness, vs_span,...
%             fuse_num_plies, bulkheads_num, bulkhead_thickness, bulkhead_area_fraction, plywood_density, fuse_height, fuse_width, fuse_length,...
%             tail_boom_num_plies, tail_boom_length, wing_fudge_factor, hs_fudge_factor, vs_fudge_factor, fuse_fudge_factor, tailboom_fudge_factor);
%            
%     v_guess = 20;
%     v = 10;
%     while abs(v_guess - v) > 0.1
%         v_guess = v;
%         
%         wing_surface_area = wing_area*2;
% %         saTail = tail_area_total*2;
%         hs_surface_area = hs_area*2;
%         vs_surface_area = vs_area*2;
%         %get the zero-lift drag coeff for this S and velocity
%         cd0 = getZeroLiftDrag(density, viscosity, v, ...
%                    wing_area, wing_surface_area, wing_chord, fuse_surface_area,fuse_length, nose_surface_area,nose_length,...
%                    hs_surface_area,hs_chord, vs_surface_area,vs_chord, tail_boom_surface_area,tail_boom_length);
% 
%         % calculate velocity and lift coefficient
%         v = -1;
%         [wing_area, v, cl] = minDragEq(v,...
%           density, wing_area, cd0, wing_span, weight, taper_ratio, sweepAngle);
%         wing_surface_area = wing_area*2;
%       
%         % calculate tail areas and length and tail boom length (and therefore area)
% %         [hs_area, vs_area, tail_boom_length] = find_tail_size(wing_span, ...
% %         wing_area, wing_surface_area, wing_chord, C_HT, C_VT, density, viscosity, v(end));
%         
% 
%     end
% 
%     
%     component_moi = calc_component_moi(mass_vec,height_vec,width_vec);
%     
%     quarter_chord = [(wing_chord/2)+(.2*fuse_length),0,0.05];
%     central_moment = calc_central_moment(component_moi,weight_vec,position_vec,quarter_chord);
%     
%     %get induced drag so we can later get cd
%     k = getK(AR, taper_ratio, sweepAngle);
%     cdi = cl^2*k;
%     %get cd so we can later get cl/cd
%     cd = cd0 + cdi;
%     drag = 1/2*density*v^2*cd*wing_area;
%     %get cl/cd for performance comparison
%     clOverCd = cl/cd;
%     
%     Nx = 100;
%     Ny = 0;
%     Nxy = 0;
%     Mx = 0;
%     My = 0;
%     Mxy = 0;
%     delta_T = 0;
%     airfoil_thickness = .1;
%     thetas = 90;
%     rad_or_deg = "deg";
%     thicknesses = ones(length(thetas),1)*thickness_per_ply;
%     Rm = -lift*(wing_span/4); % reaction moment at root
%     % t_airfoil = airfoil thickness
%     Fx = -Rm/airfoil_thickness; %Force
%     l = 0.8*wing_chord - 0.2*wing_chord; % length of wing box
%     Nx = Fx/l; % in-plane stress
%     Rm_y = drag*(wing_span/4); % reaction moment at root in the y direction
%     wing_box_distance  = 0.8*wing_chord - 0.2*wing_chord; % distance between two spars 
%     Nxy = Rm_y/(wing_box_distance*(wing_span/2));
%     mech_loading = [Nx;Ny;Nxy;Mx;My;Mxy];
%     
%     [stresses_bot, stresses_top,z_all,mid_strains_and_curvatures,...
%      thermal_loading, ABD,MS, failed_plies, failed_side, failed_z,...
%      fail_mode, fail_tcs] = structural_model(mat_props, thetas, ...
%      rad_or_deg, thicknesses, mech_loading, delta_T, cte_vec,...
%      mat_strengths_t,fail_crit,mat_strengths_c, SF, print_output);
%  
%     min_ms = min(MS(:));
%     
%     ms_constraint = -min_ms;
%     ms_constraint_rho = 10;
%     
%     constraint_vec = [ms_constraint]; %make sure to make vertical
%     constraint_rho_vec = [ms_constraint_rho];
%     active_constraint_vec = constraint_vec(constraint_vec > 0);
%     active_constraint_rho_vec = constraint_rho_vec(constraint_vec > 0);
%     penalty_scaling_factors = diag(active_constraint_rho_vec);
%     objective_penalty = 1/2 * active_constraint_vec' * penalty_scaling_factors * active_constraint_vec;
%     if isempty(objective_penalty)
%         objective_penalty = 0;
%     end
% 
%     if i == 1
%         f_x = drag + objective_penalty;
%     else
%         f_x_plus_h = drag + objective_penalty;
%         gradient(i-1) = (f_x_plus_h - f_x)/h;
%         x = x - x_step;
%     end
%     end     
%     x_history(:,iter_num) = x;
%     objective_history(iter_num) = f_x;
%     gradient_norm(iter_num) = gradient;
% end
% figure(4);
% plot(wing_chords,objective_history)
% xlabel("Wing Chord"); ylabel("Objective"); 
% title("Gradient-Based Iteration: Objective vs. Wing Chord");
% figure(5);
% plot(wing_chords,gradient_norm)
% xlabel("Wing Chord"); ylabel("Gradient Norm"); 
% title("Gradient-Based Iteration: Gradient Norm vs. Wing Chord");




