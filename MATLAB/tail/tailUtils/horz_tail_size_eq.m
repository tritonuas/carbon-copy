function [C_HT] = horz_tail_size_eq(x_ach, x_cg, S_h, c, S, C_HT, hasx_ach, hasx_cg, hasS_h, hasc, hasS,...
    hasC_HT)

if ~hasx_ach
    x_ach = C_HT*(c*S)/S_h+x_cg;

elseif ~hasx_cg
    x_cg = -C_HT*(c*S)/S_h+x_ach;
    
elseif ~hasS_h
    S_h = C_HT*(c*S)/(x_ach-x_cg);
    
elseif ~hasc
    c = (x_ach - x_cg)*S_h/(S*C_HT);

elseif ~hasS
    S = (x_ach - x_cg)*S_h/(c*C_HT);
    
elseif ~hasC_HT
    C_HT = (x_ach-x_cg)*S_h/(c*S);
    
end