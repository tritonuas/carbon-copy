%@param  

function [wingspan, S_v, wingarea, x_acv, x_cg, C_VT] = vert_tail_size_eq(x_acv, x_cg, S_v, wingspan, wingarea, C_VT, hasx_acv, hasx_cg, hasS_v, haswingspan, hasS,...
    hasC_VT)

if ~hasx_acv
    x_acv = C_VT*wingspan*wingarea/S_v+x_cg;
    
elseif ~hasx_cg
    x_cg = -C_VT*wingspan*wingarea/S_v+x_acv;
    
elseif ~hasS_v
    S_v = C_VT*wingspan*wingarea/(x_acv-x_cg);
    
elseif ~haswingspan
    wingspan = (x_acv-x_cg)*S_v/C_VT*wingarea;
    
elseif ~hasS
    wingarea = (x_acv-x_cg)*S_v/C_VT*wingspan;
    
elseif ~hasC_VT
    C_VT = (x_acv-x_cg)*S_v/(wingarea*wingspan);
   
end
    





