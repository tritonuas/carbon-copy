%% testing get cf for new plane design
velocity = 20;
viscosity = 1.825*10^-5;
density = 1.225;
chord = 1;

RE_CRIT_TRANS = 1*10^5;
RE_CRIT_TURB = 1*10^5;

xCritTrans = RE_CRIT_TRANS*viscosity/(velocity*density);
xCritTurb = RE_CRIT_TURB*viscosity/(velocity*density);

lamFunc = @(x)(1.328/sqrt(getRe(density,velocity,x,viscosity)));
transFunc = @(x)(0.455/(log10(density*velocity*x/viscosity)^2.58) ...
    - 1700/(density*velocity*x/viscosity));
turbFunc = @(x)(0.455/(log10(density*velocity*x/viscosity)^2.58));

fullLamComponent = lamFunc(xCritTrans);
fullTransComponent = transFunc(xCritTurb) - transFunc(xCritTrans);

l = linspace(0.01,chord,500);
cfValues = [1:10;1:10];

% for i = 1:length(l)
%     if l(i) > xCritTurb
%         cf = (fullLamComponent*xCritTrans + ...
%             turbFunc(l(i))*(l(i)-xCritTrans))...
%             /l(i);
%         cfValues(2,i) = 3;
% %     elseif l(i) > xCritTrans
% %         cf = (fullLamComponent*xCritTrans + ...
% %             transFunc(l(i))*(l(i)-xCritTrans))...
% %             /chord;
% %         cfValues(2,i) = 2;
%     else
%         cf = lamFunc(l(i));
%         cfValues(2,i) = 1;
%     end
%     cfValues(1,i) = cf;
% end

for i = 1:length(l)
    if l(i) > xCritTurb
        cf = (lamFunc(xCritTurb)*xCritTurb + ...
            turbFunc(l(i))*(l(i)-xCritTurb))...
            /l(i);
        cfValues(2,i) = 3;
%     elseif l(i) > xCritTrans
%         cf = (fullLamComponent*xCritTrans + ...
%             transFunc(l(i))*(l(i)-xCritTrans))...
%             /chord;
%         cfValues(2,i) = 2;
    else
        cf = lamFunc(l(i));
        cfValues(2,i) = 1;
    end
    cfValues(1,i) = cf;
end
drag = l.*cfValues(1,:);


figure(1);
plot(l,cfValues(1,:))
figure(2)
plot(l,drag)