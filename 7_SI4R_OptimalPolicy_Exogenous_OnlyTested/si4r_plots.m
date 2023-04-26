%%%%%%%% SI4R-Macro Plots %%%%%%%%

%plotting
ia=4;ib=4;fsize=8;
horz=HH;
time=0:1:horz-1;
 
figure;
subplot(ia,ib,1)
plot(time,I(1:horz),'b-','LineWidth',2);
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
plot(time,R(1:horz),'b-','LineWidth',2);
box off;
title('Recovered, R','FontSize',fsize);
set(gca,'FontSize',fsize);
 
subplot(ia,ib,4)
plot(time,D(1:horz),'b-','LineWidth',2);
box off;
title('Deaths, D','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time,Ia0(1:horz),'b-','LineWidth',2);
box off;
title('Infected Asymptomatic Untested , Ia-','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,6)
plot(time,Ia1(1:horz),'b-','LineWidth',2);
box off;
title('Infected Asymptomatic Tested , Ia+','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,7)
plot(time,Ib0(1:horz),'b-','LineWidth',2);
box off;
title('Infected Symptomatic Untested , Ib-','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,8)
plot(time,Ib1(1:horz),'b-','LineWidth',2);
box off;
title('Infected Symptomatic Tested , Ib+','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,9)
plot(time,R0(1:horz),'b-','LineWidth',2);
box off;
title('Recovered 0 , R-','FontSize',fsize);
ylabel('Share of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,10)
plot(time,R1(1:horz),'b-','LineWidth',2);
box off;
title('Recovered 1 , R+','FontSize',fsize);
set(gca,'FontSize',fsize);

suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SI4Rmacro_epidemic_simulation_fig1

 
 

%plotting
ia=4;ib=2;fsize=8;
horz=HH;
time=0:1:horz-1;

figure;
 
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

subplot(ia,ib,3)
plot(time,100*(cs(1:horz)-cr1ss)/cr1ss,'b-','LineWidth',2); hold on
plot(time,100*(cia0(1:horz)-cr1ss)/cr1ss,'r--','LineWidth',2); hold on
plot(time,100*(cr0(1:horz)-cr1ss)/cr1ss,'k-.','LineWidth',2); hold on
box off;
axis([0 horz-1 -(1-phii+0.02)*100 2]);
title('Consumption DOUBTERS','FontSize',fsize);

xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Cons. Susceptibles','Cons. Asymptomatic Infected Untested','Cons. Recovered Unaware','Location','best');
legend boxoff
 
subplot(ia,ib,4)
plot(time,100*(ns(1:horz)-nr1ss)/nr1ss,'b-','LineWidth',2); hold on
plot(time,100*(nia0(1:horz)-nr1ss)/nr1ss,'r--','LineWidth',2); hold on
plot(time,100*(nr0(1:horz)-nr1ss)/nr1ss,'k-.','LineWidth',2); hold on
box off;
title('Hours by DOUBTERS Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Hours Susceptibles','Hours Asymptomatic Infected Untested','Hours Recovered Unaware','Location','best');
legend boxoff

subplot(ia,ib,5)
plot(time,100*(cib0(1:horz)-cr1ss)/cr1ss,'b-','LineWidth',2); hold on
plot(time,100*(cib1(1:horz)-cr1ss)/cr1ss,'r--','LineWidth',2); hold on
plot(time,100*(cia1(1:horz)-cr1ss)/cr1ss,'k-.','LineWidth',2); hold on
box off;
axis([0 horz-1 -(1-phii+0.02)*100 2]);
title('Consumption by INFECTED Type','FontSize',fsize);

xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Cons. Symptomatic Infected Untested','Cons. Symptomatic Infected Tested','Asymptomatic Infected Tested','Location','best');
legend boxoff

subplot(ia,ib,6)
plot(time,100*(nib0(1:horz)-nr1ss)/nr1ss,'b-','LineWidth',2); hold on
plot(time,100*(nib1(1:horz)-nr1ss)/nr1ss,'r--','LineWidth',2); hold on
plot(time,100*(nia1(1:horz)-nr1ss)/nr1ss,'k-.','LineWidth',2); hold on
box off;
title('Hours by INFECTED Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Hours Symptomatic Infected Untested','Hours Symptomatic Infected Tested','Hours Asymptomatic Infected Tested','Location','best');
legend boxoff

subplot(ia,ib,7)
plot(time,100*(cr0(1:horz)-cr1ss)/cr1ss,'b-','LineWidth',2); hold on
plot(time,100*(cr1(1:horz)-cr1ss)/cr1ss,'r--','LineWidth',2); hold on
box off;
axis([0 horz-1 -(1-phii+0.02)*100 2]);
title('Consumption by RECOVERED Type','FontSize',fsize);

xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Cons. Recovered Unaware','Cons. Recovered Aware','Location','best');
legend boxoff

subplot(ia,ib,8)
plot(time,100*(nr0(1:horz)-nr1ss)/nr1ss,'b-','LineWidth',2); hold on
plot(time,100*(nr1(1:horz)-nr1ss)/nr1ss,'r--','LineWidth',2); hold on
box off;
title('Hours by RECOVERED Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Hours Recovered Unaware','Hours Recovered Aware','Location','best');
legend boxoff

suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SI4Rmacro_epidemic_simulation_fig2
 
 
 
figure;
subplot(6,1,1)
plot(time,Us(1:HH),'b-','LineWidth',2); hold on;
plot(time,Uia0(1:HH),'r--','LineWidth',2);hold on;
plot(time,Ur0(1:HH),'k-.','LineWidth',2);hold on;
plot(time,U(1:HH),'m:','LineWidth',2);hold on
set(gca,'FontSize',fsize);
legend('Susceptibles','Infected Asymptomatic Untested','Recovered Unaware','Total','Location','best');
set(gca,'FontSize',fsize);
title('Present Value Utility','FontSize',fsize);
legend boxoff;

subplot(6,1,2)
plot(time,Uib0(1:HH),'b-','LineWidth',2); hold on;
plot(time,Uib1(1:HH),'r--','LineWidth',2);hold on;
plot(time,Uia1(1:HH),'k-.','LineWidth',2);hold on;
plot(time,U(1:HH),'m:','LineWidth',2);hold on
set(gca,'FontSize',fsize);
legend('Symptomatic Infected Untested','Symptomatic Infected Tested','Asymptomatic Infected Tested','Total','Location','best');
set(gca,'FontSize',fsize);
title('Present Value Utility','FontSize',fsize);
legend boxoff;

subplot(6,1,3)
plot(time,Ur0(1:HH),'b-','LineWidth',2); hold on;
plot(time,Ur1(1:HH),'r--','LineWidth',2);hold on;
plot(time,U(1:HH),'m:','LineWidth',2);hold on
set(gca,'FontSize',fsize);
legend('Recovered Unaware','Recovered Aware','Total','Location','best');
set(gca,'FontSize',fsize);
title('Present Value Utility','FontSize',fsize);
legend boxoff;
 
subplot(6,1,4)
plot(time,0*muc0,'k:','LineWidth',1.5); hold on
plot(time,muc0,'b-','LineWidth',2); hold on
set(gca,'FontSize',fsize);
title('Containment Policy, \mu0_c','FontSize',fsize);
ylabel('%','FontSize',fsize)

subplot(6,1,5)
plot(time,0*muc1,'k:','LineWidth',1.5); hold on
plot(time,muc1,'b-','LineWidth',2); hold on
set(gca,'FontSize',fsize);
title('Containment Policy, \mu1_c','FontSize',fsize);
ylabel('%','FontSize',fsize)

suptitle('The Evolution of an Epidemic');
orient portrait
print -dpdf -fillpage SI4Rmacro_epidemic_simulation_fig3