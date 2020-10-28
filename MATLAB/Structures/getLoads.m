function c = getLoads(Wtotal,Wwing,L,n)

% length = 3.65; 
% taper = 0.45;
% lift = 365;
% rchord = 0.3;

load = @(x,t) 2*Wtotal*n*(sqrt(L^2 - x^2))/(L^2*pi()) - (Wwing*n*(x/L))/L

span = linspace(0,L);

%y2 = @(x) x;

shear = arrayfun(@(span) integral(load, 0, span, 'ArrayValued', true), span);
%moment = arrayfun(@(span) integral2(load, 0, y2, 0, span, 'ArrayValued', true), span);            

figure(1)
plot(span,shear)
max(shear)

%figure(2)
%plot(span,moment)


figure(3)
fplot(load, [0 L])

end
