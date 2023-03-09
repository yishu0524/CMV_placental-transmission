filename = 'healthy populationflow.csv';

%read data
    plasma = readtable(filename);column_id=16;
%%
    monkey_data = plasma(:,[1,column_id]);
    
    monkey_data.Properties.VariableNames{1}='days';
    monkey_data.Properties.VariableNames{2}='viralLoad';
   %%
    T =166;
        
    l=1.2; %width of the placenta is half of l
    D=0.0069;
    nx=1000;
    nt=1000;
    p=1/D;
    x = linspace(0,l,nx); %2–2.5 cm (0.8–1 inch) in thickness, 
    t = linspace(0,T,nt); %if we run simulation for 20 days
    m = 0;

    
    %%
    %% pass the monkey data into pde
    %%
    s1=monkey_data.days;
    s2=monkey_data.viralLoad;
    


    options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7) 
    sol = pdepe(0,@unsatpde,@unsatic,@unsatbc,x,t,options,s1,s2,p,D); 
    u= sol(:,500:1000,1);
    x=x(500:1000)-x(500);
    %%piecewise(x < 0,-1,x > 0,1)
    syms placenta_growth(time)
    placenta_growth(time)=piecewise( time<=40,666.93*time,time>40,132434./(1+ exp(-(0.028*(time-90)))));
    placenta_growth_vals = placenta_growth(t);
    %%
    [ux,ut]=gradient(u,l/nx,T/nt);

    virus_numeric=trapz(ux(:,1) .*transpose(placenta_growth_vals) )*T/nt*D;

    %%
    figure2=figure(2)
    plot(t,placenta_growth_vals)
    ylim([0 11000])
    xlim([0 166])
    title('Placental growth over pregancy')
    xlabel('Time t (days)')
    ylabel('Placental Surface Area (mm^2)')
    legend('Placental Growth Function')
    saveas(figure2, 'placenta_growth.png')
    %%
    Q = double(subs(ux(:,1).* D.*transpose(placenta_growth_vals)));
    %virus_total=trapz(t,Q);
    %avg=trapz(t,placenta_growth_vals)/T;
    %%
    figure3=figure (3)
    plot(t,Q,'b')
    %hold on
    %plot(t, ux(:,1).* D*avg,'g')
    %hold off
    ylim([0 inf])
    xlim([0 166])
    title(monkey_ID)
    xlabel('Time t (days)')
    ylabel('CMV DNA Viral Load copies/(microliter*mm*day)')
    saveas(figure3, strcat(monkey_ID,'flux.png'))
    %%
    figure4= figure (4)
    surf(x,t,u,'FaceAlpha',0.5,'EdgeColor','none')
    hold on 
    scatter3(repelem(l/2,size(monkey_data,1)),monkey_data.days,monkey_data.viralLoad,'or','filled')
    xlabel('Distance x (mm)')
    ylabel('Time t (days)')
    ylim([0,166])
    zlim([0,inf])
    xlim([0,l/2])
    zlabel('CMV DNA Viral Load copies per microliter')
    ax=gca;
    ax.YDir="reverse"
    hold off
    title(monkey_ID)
    saveas(figure4, strcat(monkey_ID,'flow.png'))
    %%
    t=transpose(t);
    Table=table(Q,t);
    Table.Properties.VariableNames{2}='days';
    %%
    Table.Properties.VariableNames{1}=monkey_ID;
    writetable(Table, strcat(monkey_ID,'flow.csv'))
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