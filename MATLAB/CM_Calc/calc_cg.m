
function cg=calc_cg(weight_vec,position_vec)
    sums = (weight_vec*position_vec);
    cg = sums/(sum(weight_vec));
end