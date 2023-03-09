%viral load on mother side
initial_values = [1e-3, 0, 0,5e2,0]; % (V, E, RI, RS, RL) 
 [t_m,y] = ode15s(@mother,[0 166],initial_values); 
% xlabel('Time t');
% ylabel('Solution y');
% set(gca, 'YScale', 'log')
% legend('V')
% x = 0:.25:5;
% y = x.^2;
y_m=y(:,1);
epsilonS='No immune suppression';
V = @(z) interp1(t_m,y(:,1)*1e3,z);  % Lookup table function

t2 = 0:1:280;  % Just to compare func to data.
figure (1)
plot(t_m,y(:,1)*1e3,'sb',t2,V(t2),'*r')

%%
l=1; %width of the placenta
T=166;
D=0.0069;
nx=2000;
nt=1000;
p=0.5/D;
x = linspace(0,l,nx); %2–2.5 cm (0.8–1 inch) in thickness, 
t = linspace(0,T,nt); %if we run simulation for 20 days

m = 0;
s1=t_m;s2=y_m;


options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7) 
sol = pdepe(0,@unsatpde,@unsatic,@unsatbc,x,t,options,s1,s2,p,D); 
u= sol(:,1:1000,1);
figure (2)
x=x(1:1000);
plot(x,u(end,:))
%%
placenta_growth=132434./(1+ exp(-(0.028*(t-90))));


[ux,ut]=gradient(u,l/nx,T/nt);
Q = ux(:,1).* D.*transpose(placenta_growth);
virus_total=trapz(t,Q)
avg=trapz(t,placenta_growth)/T;
t=transpose(t);
Table1=table(t,Q);
Table1.Properties.VariableNames{1}='days';
Table1.Properties.VariableNames{2}=epsilonS;
writetable(Table1, strcat(epsilonS,'.csv'))

    figure3=figure (3)
    plot(t,Q,'b')
    %hold on
    %plot(t, ux(:,1).* D*avg,'g')
    %hold off
    ylim([0 inf])
    xlim([0 166])

    xlabel('Time t (days)')
    ylabel('CMV DNA Viral Load copies/(microliter*mm*day)')


virus_numeric=trapz(ux(:,1) .*transpose(placenta_growth) )*T/nt*D
%%

figure4= figure (4)
    surf(x,t,u,'FaceAlpha',0.5,'EdgeColor','none')
    hold on 
    
    xlabel('Distance x (mm)')
    ylabel('Time t (days)')
    ylim([0,166])
    zlim([0,inf])
    zlabel('CMV DNA Viral Load copies per microliter')
    ax=gca;
    ax.YDir="reverse"
%%

%%
function [c,f,s] = unsatpde(x,t,u,DuDx,s1,s2,p,D) 
c =1/D;% 1/(1.8*10^-3) ; %inverse of the diffusion coefficent
f = DuDx;
s =-p *u; %death of the cmv %death rate is 1, 2 virus day per day
end 
% ------------------------------------------------------------------------- 
function u0 = unsatic(x,s1,s2,p,D) 
u0 = 0; %initial is 0
end 
% ------------------------------------------------------------------------- 
function [pl,ql,pr,qr] = unsatbc(xl,ul,xr,ur,t,s1,s2,p,D) 
pl = ul; %u(l,t)=1
ql = 0; 
pr = ur-interp1(s1,s2,t,'linear','extrap'); %u(r,t)=0
qr = 0;
end 