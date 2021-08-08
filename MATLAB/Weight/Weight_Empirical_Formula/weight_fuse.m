

function [w_fuse] = weight_fuse(fus_wet, tail_len, fuse_len, fus_str_depth, load_fact_ult, gross, q) 

    mti_weight = 0.224809; % metric (N) to imperial (poundforce) 
    mti_area = 10.7639; % m^2 to ft^2
    mti_press = 0.0208854;  % pascals to pound force per square foot
    mti_len = 3.28084; % m to ft

    w_fuse = 0.052*(fus_wet*mti_area)^(1.086)*...
    (load_fact_ult*gross*mti_weight)^(0.177)*...
	(tail_len*mti_len)^(-0.051)*(fuse_len/fus_str_depth)^(-0.072)*...
	(q*mti_press)^(0.241)*1/mti_weight;
end 