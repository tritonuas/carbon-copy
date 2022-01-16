objective = @(x)x'*x
x = [15;20]
alpha = 0.1
h = 1e-5
for i = 1:100
    gradient = finite_difference(objective,x,h);
    x = gradient_descent_optimizer(gradient,alpha,x)
end