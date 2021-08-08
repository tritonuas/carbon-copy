
function [gross] = compute_weight(S, battery, AR, taperRatio, sweepAngle, load_fact_ult,...
                 S_h, htail_ar, htail_taper, htail_sweep, htail_t_over_c,...
                 S_v, vtail_t_over_c, vtail_sweep, vtail_ar, vtail_taper, ...
                 fus_wet, tail_len, fuse_len, fus_str_depth, q, wing_t_over_c, empty_struct_fudge,...
                 misc_weight_fudge, payload)

    gross = 30;
    while (1)
      w_wing = weight_wing(S, battery, AR, taperRatio, sweepAngle, load_fact_ult, gross, q, wing_t_over_c);
      w_htail = weight_htail(S_h, htail_ar, htail_taper, htail_sweep, load_fact_ult, gross, htail_t_over_c, q);
      w_vtail = weight_vtail(load_fact_ult, gross, S_v, vtail_t_over_c, vtail_sweep, vtail_ar, vtail_taper, q);
      w_fuse = weight_fuse(fus_wet, tail_len, fuse_len, fus_str_depth, load_fact_ult, gross, q); 
      empty = empty_struct_fudge*(w_wing + w_htail + w_vtail + w_fuse);
      new_gross = (1 + misc_weight_fudge) * (payload + empty + battery); 
      if abs(new_gross - gross) < 0.1
        break
      end
      gross = new_gross;
    end
    gross = gross;
end

