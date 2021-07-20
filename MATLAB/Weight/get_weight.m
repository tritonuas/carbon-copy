% get weight 

mti_weight = 0.224809; % metric (N) to imperial (poundforce) 
mti_area = 10.7639; % m^2 to ft^2
mti_press = 0.0208854;  % pascals to pound force per square foot



function [W_wing] = Weight_wing(S, battery, AR, mti_press, taperRatio, sweepAngle, load_fact_ult, gross)

    W_wing = 0.036*(S*mti_area)^(0.758)*(battery*mti_weight)^(0.0035)*...
    (AR/np.cos(sweepAngle)^2)^(0.6)*(q*mti_press)^(0.006)*taperRatio^(0.04)...
    *(100*wing_t_over_c/np.cos(sweepAngle))^(-0.3)*(load_fact_ult*gross*mti_weight)^(0.49)*1/mti_weight;

end

function [W_htail] = Weight_htail(S_h, htail_ar, htail_taper, htail_sweep, load_fact_ult, gross, htail_t_over_c) 

    W_htail = 0.016*(load_fact_ult*gross*mti_weight)^(0.414)*...
    (q*mti_press)^(0.168)*(S_h*mti_area)^(0.896)*...
    (100*htail_t_over_c/np.cos(htail_sweep))^(-0.12)*...
    (htail_ar/np.cos(htail_sweep)^2)^(0.043)*htail_taper^(-0.02)*...
    1/mti_weight;
end


function [W_vtail] = Weight_vtail(load_fact_ult,gross, S_v, vtail_t_over_c, vtail_sweep, vtail_ar, vtail_taper)

    W_vtail = 0.073*(load_fact_ult*gross*mti_weight)^(0.376)*...
    (q*mti_press)^(0.122)*(S_v*mti_area)^(0.873)*...
    (100*vtail_t_over_c/np.cos(vtail_sweep))^(-0.49)*...
    (vtail_ar/np.cos(vtail_sweep)^2)^(0.357)*vtail_taper^(0.039)*...
    1/mti_weight;
end


function [W_fuse] = Weight_fuse(fus_wet, tail_len, fuse_len, fus_str_depth, load_fact_ult) 

    W_fuse = 0.052*(fus_wet*mti_area)^(1.086)*...
    (load_fact_ult*gross*mti_weight)^(0.177)*...
	(tail_len*unit.mti_len)^(-0.051)*(fuse_len/fus_str_depth)^(-0.072)*...
	(q*mti_press)^(0.241)*1/mti_weight;
end 


function [gross] = compute_weight(W_wing, W_htail, W_vtail, W_fuse, payload, battery, misc_weight_fudge)
empty = empty_struct_fudge*(W_wing + W_htail + W_vtail + W_fuse)
new_gross = (1 + misc_weight_fudge) * (payload + empty + battery) 
   
while (1)
    if abs(new_gross - gross) < 0.1
        gross = new_gross;
        break
        gross = new_gross;
    end
end
end
