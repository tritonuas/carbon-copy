%%Checks and prints failure in the each of the plies in a laminate
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%@param stresses_top the stresses at the top of the plies
%@param stresses_bot the stresses at the bot of the plies
%@param z_all the sorted z-coordinate of the top and bottom of each ply
%@param fail_crit a string for the failure criterion to use
%@param mat_strengths_t the vector or matrix of material tensile strengths
%@param mat_strengths_c the vector or matrix of mat compression strengths
%@param SF the safety factor
%@param will_print a boolean for whether the formatted output will print
%@return MS the matrix of margins of safety with each col as a ply
%@return failed_plies a boolean vector reporting which plies have failed
function [MS, failed_plies, failed_side, failed_z, fail_mode, fail_tcs] ...
    = report_ply_margins(stresses_bot, stresses_top, z_all, ...
    fail_crit, mat_strengths_t, mat_strengths_c, SF, print_output)

%allowing the input of just a vector rather than a copied matrix
if length(mat_strengths_t(1,:)) == 1
    mat_strengths_t = mat_strengths_t*ones(1,length(stresses_top(1,:)));
end
if length(mat_strengths_c(1,:)) == 1
    mat_strengths_c = mat_strengths_c*ones(1,length(stresses_top(1,:)));
end

%usage of faiilure criterion
if fail_crit == "max_stress"
    [MS, stress_signs] = max_stress_fail_crit(stresses_bot, ...
    stresses_top, mat_strengths_t, mat_strengths_c, SF);
elseif fail_crit == "tsai-wu"
    disp("Tsai-Wu has yet to be implemented.")
else
    disp("Not a valid failure criterion (either not real or not"...
        +"implemented).")
end

%output
if(print_output)
disp("------------Laminate Report------------")
end
%preallocation
num_plies = 1:length(mat_strengths_t(1,:));     %numbered plies (1-num)
%
ply_num = sort([num_plies num_plies]);  
failed_plies = [];
failed_z = [];
failed_side = [];
fail_mode = [];
fail_tcs = [];
for i = 1:length(MS(1,:))
    for j = 1:3
        if MS(j,i) < 0
            if isempty(failed_plies)
                failed_plies = ply_num(i);
                failed_z = z_all(i);
                if mod(i,2) == 0
                    failed_side = "Top";
                else
                    failed_side = "Bottom";
                end
                if j == 1
                    fail_mode = "Fiber";
                else
                    fail_mode = "Matrix";
                end
                if j == 3
                    fail_tcs = "Shear";
                elseif stress_signs(j,1) == 1
                    fail_tcs = "Tension";
                else
                    fail_tcs = "Compression";
                end
            else
                failed_plies(end+1) = ply_num(i);
                failed_z(end+1) = z_all(i);
                if mod(i,2) == 0
                    failed_side(end+1) = "Top";
                else
                    failed_side(end+1) = "Bottom";
                end
                if j == 1
                    fail_mode(end+1) = "Fiber";
                else
                    fail_mode(end+1) = "Matrix";
                end
                if j == 3
                    fail_tcs(end+1) = "Shear";
                elseif stress_signs(j,1) == 1
                    fail_tcs(end+1) = "Tension";
                else
                    fail_tcs(end+1) = "Compression";
                end
            end
        end
    end
    %Fill this out with formatted table output
end
end