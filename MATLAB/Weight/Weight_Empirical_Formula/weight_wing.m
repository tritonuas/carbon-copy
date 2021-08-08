


function [w_wing] = weight_wing(S, battery, AR, taperRatio, sweepAngle, load_fact_ult, gross, q, wing_t_over_c)

    mti_weight = 0.224809; % metric (N) to imperial (poundforce) 
    mti_area = 10.7639; % m^2 to ft^2
    mti_press = 0.0208854;  % pascals to pound force per square foot

    w_wing = 0.036*(S*mti_area)^(0.758)*(battery*mti_weight)^(0.0035)*...
    (AR/cos(sweepAngle)^2)^(0.6)*(q*mti_press)^(0.006)*taperRatio^(0.04)...
    *(100*wing_t_over_c/cos(sweepAngle))^(-0.3)*(load_fact_ult*gross*mti_weight)^(0.49)*1/mti_weight;

end




