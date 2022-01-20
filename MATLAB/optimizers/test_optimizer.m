objective = @(x)x'*x
alpha = 0.1
h = 1e-5
gradient_norm = zeros(100,1)
x = zeros(2,100)
x(:,1) = [15;20]
for i = 1:100
    gradient = finite_difference(objective,x(:,i),h);
    gradient_norm(i) = norm(gradient)
    x(:,i+1) = gradient_descent_optimizer(gradient,alpha,x(:,i))
end
iterations = 1:100
plot(iterations,gradient_norm)
plot(iterations,x(1,1:length(x)-1),'b',iterations,x(2,1:length(x)-1),'r')