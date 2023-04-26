%%%%%%%%% SIR-Macro Plots %%%%%%%%

%plotting
ia=2;ib=2;fsize=8;
horz=HH;
time=0:1:horz-1;
 
figure;
subplot(ia,ib,1)
plot(time,Ib0(1:horz),'b-','LineWidth',2);
box off;
title('Infected, I','FontSize',fsize);
ylabel('Share of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);
 
subplot(ia,ib,2)
plot(time,S(1:horz),'b-','LineWidth',2);
box off;
title('Susceptibles, S','FontSize',fsize);
set(gca,'FontSize',fsize);
 
subplot(ia,ib,3)
plot(time,R1(1:horz),'b-','LineWidth',2);
box off;
title('Recovered, R','FontSize',fsize);
ylabel('Share of Initial Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
 
subplot(ia,ib,4)
plot(time,D(1:horz),'b-','LineWidth',2);
box off;
title('Deaths, D','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
 
suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig1
 
 
figure;
subplot(ia,ib,3)
plot(time,100*(cs(1:horz)-cr1ss)/cr1ss,'b-','LineWidth',2); hold on
plot(time,100*(cib0(1:horz)-cr1ss)/cr1ss,'r--','LineWidth',2); hold on
plot(time,100*(cr1(1:horz)-cr1ss)/cr1ss,'k-.','LineWidth',2); hold on
box off;
%axis([0 horz-1 -(1-phi+0.02)*100 2]);
title('Consumption by Type','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Cons. Susceptibles','Cons. Infected','Cons. Recovered','Location','best');
legend boxoff
 
subplot(ia,ib,4)
plot(time,100*(ns(1:horz)-nr1ss)/nr1ss,'b-','LineWidth',2); hold on
plot(time,100*(nib0(1:horz)-nr1ss)/nr1ss,'r--','LineWidth',2); hold on
plot(time,100*(nr1(1:horz)-nr1ss)/nr1ss,'k-.','LineWidth',2); hold on
box off;
title('Hours by Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Hours Susceptibles','Hours Infected','Hours Recovered','Location','best');
legend boxoff
 
subplot(ia,ib,1)
plot(time,0*100*(aggC(1:horz)-cr1ss)/cr1ss,'m:','LineWidth',1.5);hold on
plot(time,100*(aggC(1:horz)-cr1ss)/cr1ss,'b-','LineWidth',2);
box off;
title('Aggregate Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
set(gca,'FontSize',fsize);
 
subplot(ia,ib,2)
plot(time,0*100*(aggC(1:horz)-cr1ss)/cr1ss,'m:','LineWidth',1.5);hold on
plot(time,100*(aggH(1:horz)-nr1ss)/nr1ss,'b-','LineWidth',2);
box off;
title('Aggregate Hours, H','FontSize',fsize);
set(gca,'FontSize',fsize);
 
suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig2
 
 
 
figure;
subplot(3,1,1)
plot(time,Us(1:HH),'b-','LineWidth',2); hold on;
plot(time,Uib0(1:HH),'r--','LineWidth',2);hold on;
plot(time,Ur0(1:HH),'k-.','LineWidth',2);hold on;
plot(time,U(1:HH),'m:','LineWidth',2);hold on
set(gca,'FontSize',fsize);
legend('Susceptibles','Infected','Recovered','Total','Location','best');
set(gca,'FontSize',fsize);
title('Present Value Utility','FontSize',fsize);
legend boxoff;
 
subplot(3,1,2)
plot(time,0*muc1,'k:','LineWidth',1.5); hold on
plot(time,muc1,'b-','LineWidth',2); hold on
set(gca,'FontSize',fsize);
title('Containment Policy, \mu_c','FontSize',fsize);
ylabel('%','FontSize',fsize)

suptitle('The Evolution of an Epidemic');
orient portrait
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig3
 
 