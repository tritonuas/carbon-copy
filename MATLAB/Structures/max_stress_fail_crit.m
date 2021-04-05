%%Checks for failutre using the max stress criterion
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%@param stresses_bot the stresses at the bot of the plies
%@param stresses_top the stresses at the top of the plies
%@param mat_strengths_t the vector or mat of material properties in tension
%@param mat_strengths_c the vector or mat of material props in compression
%@param SF the safety factor
%@return MS the matrix of margins of safety with each col as a ply
function [MS, stress_signs] = max_stress_fail_crit(stresses_bot, ...
    stresses_top, mat_strengths_t, mat_strengths_c, SF)

%checking to make sure the transpose is not input
if length(mat_strengths_t(:,1)) ~= 3
    mat_strengths_t = mat_strengths_t';
end
if length(mat_strengths_c(:,1)) ~= 3
    mat_strengths_c = mat_strengths_c';
end
%checking to make sure input is proper
if length(mat_strengths_t(:,1)) ~= 3
    disp("Please input a material strengths vector or matrix that has" ...
        +"3 rows corresponding to max sigma_11;sigma_22;sigma_21");
end

%allowing the input of just a vector rather than a copied matrix
if length(mat_strengths_t(1,:)) == 1
    mat_strengths_t = mat_strengths_t*ones(1,length(stresses_top));
end
if length(mat_strengths_c(1,:)) == 1
    mat_strengths_c = mat_strengths_c*ones(1,length(stresses_top));
end

max_allow_t = mat_strengths_t./SF;
max_allow_c = mat_strengths_c./SF;

for i = 1:3
    counter = 0;
    for j = 1:length(mat_strengths_t(1,:))
        counter = counter + 1;
        if stresses_bot(i,j) > 0
            MS(i,counter) = abs(max_allow_t(i,j)/stresses_bot(i,j)) - 1;
        else
            MS(i,counter) = abs(max_allow_c(i,j)/stresses_bot(i,j)) - 1;
        end
        stress_signs(i,counter) = sign(stresses_bot(i,j));
        
        counter = counter + 1;
        if stresses_top(i,j) > 0
            MS(i,counter) = abs(max_allow_t(i,j)/stresses_top(i,j)) - 1;
        else
            MS(i,counter) = abs(max_allow_c(i,j)/stresses_top(i,j)) - 1;
        end
        stress_signs(i,counter) = sign(stresses_top(i,j));
    end
end