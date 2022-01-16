function design_variables = gradient_descent_optimizer(gradient,alpha,x_0)
delta_x = -gradient*alpha;
design_variables = x_0 + delta_x;
end