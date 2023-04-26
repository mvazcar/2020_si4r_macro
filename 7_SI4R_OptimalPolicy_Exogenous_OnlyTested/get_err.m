function [err,I,Ia0,Ia1,Ib0,Ib1,S,R,R0,R1,D,T,Pop,RnotSIRmacro,aggC,aggH,ns,nia0,nia1,nib0,nib1,nr0,nr1,cs,cia0,cia1,cib0,cib1,cr0,cr1,Us,Uia0,Uia1,Uib0,Uib1,Ur0,Ur1,U,probSD,probA0D,probR0D] = get_err(guess,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,C,N,Ur1ss,muc0,muc1,phii,pia,pib,pit)
                                                                     
%back out guesses for ns,nia0,nia1,nib0,nib1,nr0,nr1

ns=guess(1:HH);
nia0=guess(HH+1:2*HH);
nia1=guess(2*HH+1:3*HH);
nib0=guess(3*HH+1:4*HH);
nib1=guess(4*HH+1:5*HH);
nr0=guess(5*HH+1:6*HH);
nr1=guess(6*HH+1:7*HH);


%equilibrium equations

%Recovered Aware People (AWARE)
lambr1=(theta*nr1)./A;
cr1=((1+muc0).*lambr1).^(-1);
ur1=log(cr1)-theta/2*nr1.^2;
Ur1=NaN*ones(HH+1,1);
Ur1(HH+1)=Ur1ss;
for tt=HH:-1:1    
    Ur1(tt,1)=ur1(tt)+betta*Ur1(tt+1,1);    
end
Gamma=(1+muc1).*cr1-A.*nr1;

%Infected Symptomatic Untested People (AWARE)
lambib0=(theta*nib0)./(phii*A);
cib0=((1+muc0).*lambib0).^(-1);
uib0=log(cib0)-theta/2*nib0.^2;

%Infected Symptomatic Tested People (AWARE)
lambib1=(theta*nib1)./(phii*A);
cib1=((1+muc1).*lambib1).^(-1);
uib1=log(cib1)-theta/2*nib1.^2;

%Infected Asymptomatic Tested People (AWARE)
lambia1=(theta*nia1)./(A);
cia1=((1+muc1).*lambia1).^(-1);
uia1=log(cia1)-theta/2*nia1.^2;

% Doubters - Susceptible People (UNAWARE)
cs=1./(1+muc0).*(A.*ns+Gamma);
us=log(cs)-theta/2*ns.^2;

% Doubters - Infected Untested People (UNAWARE)
cia0=1./(1+muc0).*(A.*nia0+Gamma);
uia0=log(cia0)-theta/2*nia0.^2;

% Doubters - Recovered 0 (UNAWARE)
cr0=1./(1+muc0).*(A.*nr0+Gamma);
ur0=log(cr0)-theta/2*nr0.^2;

%pre-allocate
Pop=NaN*ones(HH+1,1);

Ta0=NaN*ones(HH,1);
Tb0=NaN*ones(HH,1);
T=NaN*ones(HH,1);

Ia0=NaN*ones(HH+1,1);
Ia1=NaN*ones(HH+1,1);
Ib0=NaN*ones(HH+1,1);
Ib1=NaN*ones(HH+1,1);
I=NaN*ones(HH+1,1);

S=NaN*ones(HH+1,1);
D=NaN*ones(HH+1,1);

R0=NaN*ones(HH+1,1);
R1=NaN*ones(HH+1,1);
R=NaN*ones(HH+1,1);

%initial conditions
Pop(1)=pop_ini;

Ia0(1)=i_ini*pia;
Ia1(1)=0;
Ib0(1)=i_ini*pib;
Ib1(1)=0;
I(1)=i_ini;

S(1)=Pop(1)-I(1);
D(1)=0;

R0(1)=0;
R1(1)=0;
R(1)=0;

%iterate on SI4R model equations
for j=1:1:HH
    Ta0(j,1)=pia*(pis1*S(j)*cs(j)*(Ia0(j)*cia0(j)+Ia1(j)*cia1(j)+Ib0(j)*cib0(j)+Ib1(j)*cib1(j))+pis2*S(j)*ns(j)*(Ia0(j)*nia0(j)+Ia1(j)*nia1(j)+Ib0(j)*nib0(j)+Ib1(j)*nib1(j))+pis3*S(j)*(Ia0(j)+Ia1(j)+Ib0(j)+Ib1(j)));
    Tb0(j,1)=pib*(pis1*S(j)*cs(j)*(Ia0(j)*cia0(j)+Ia1(j)*cia1(j)+Ib0(j)*cib0(j)+Ib1(j)*cib1(j))+pis2*S(j)*ns(j)*(Ia0(j)*nia0(j)+Ia1(j)*nia1(j)+Ib0(j)*nib0(j)+Ib1(j)*nib1(j))+pis3*S(j)*(Ia0(j)+Ia1(j)+Ib0(j)+Ib1(j)));
    T(j,1)=Ta0(j)+Tb0(j);
     
    S(j+1,1)=S(j)-Ta0(j)-Tb0(j);
    
    Ia0(j+1,1)=Ia0(j)+Ta0(j)-(pit+pir)*Ia0(j);
    Ia1(j+1,1)=Ia1(j)+pit*Ia0(j)-pir*Ia1(j);
    Ib0(j+1,1)=Ib0(j)+Tb0(j)-(pit+pir+pid)*Ib0(j);
    Ib1(j+1,1)=Ib1(j)+pit*Ib0(j)-(pir+pid)*Ib1(j);
    I(j+1,1)=Ia0(j+1)+Ia1(j+1)+Ib0(j+1)+Ib1(j+1);

    R0(j+1,1)=R0(j)+pir*Ia0(j);
    R1(j+1,1)=R1(j)+pir*(Ia1(j)+Ib0(j)+Ib1(j));    
    R(j+1,1)=R0(j+1)+R1(j+1);

    D(j+1,1)=D(j)+pid*(Ib0(j)+Ib1(j));

    Pop(j+1,1)=Pop(j,1)-pid*(Ib0(j)+Ib1(j));

end

%Infected Symptomatic Tested People (continued)
Uib1=NaN*ones(HH+1,1); 
%Uib1(HH+1)=(1-0)^HH*Uib1ss+(1-(1-0)^HH)*Ur1ss;%terminal condition 
Uib1(HH+1)=Uib1ss;% Terminal Condition
%Ui(HH+1)=Ur1ss;%
for tt=HH:-1:1    
    %Ui(tt,1)=ui(tt)+(1-0)*betta*((1-pir-pid)*Ui(tt+1,1)+pir*Ur(tt+1,1))+0*betta*Ur(tt+1,1);    
    Uib1(tt,1)=uib1(tt)+betta*((1-pir-pid)*Uib1(tt+1,1)+pir*Ur1(tt+1,1));    
end

%Infected Symptomatic Untested People (continued)
Uib0=NaN*ones(HH+1,1); 
%Uib1(HH+1)=(1-0)^HH*Uib1ss+(1-(1-0)^HH)*Ur1ss;%terminal condition 
Uib0(HH+1)=Uib1ss;% Terminal Condition
%Ui(HH+1)=Ur1ss;%
for tt=HH:-1:1    
    %Ui(tt,1)=ui(tt)+(1-0)*betta*((1-pir-pid)*Ui(tt+1,1)+pir*Ur(tt+1,1))+0*betta*Ur(tt+1,1);    
    Uib0(tt,1)=uib0(tt)+betta*((1-pir-pid-pit)*Uib0(tt+1,1)+pir*Ur1(tt+1,1)+pit*Uib1(tt+1,1));    
end

%Infected Asymptomatic Tested People - Continued
Uia1=NaN*ones(HH+1,1); 
Uia1(HH+1)=Uia1ss;% Terminal Condition
for tt=HH:-1:1    
    Uia1(tt,1)=uia1(tt)+betta*((1-pir-pid)*Uia1(tt+1,1)+pir*Ur1(tt+1,1));    
end

% DOUBTERS - Continued 

% Pre-allocate Conditional Probabilites

probSD=NaN*ones(HH+1,1);
probA0D=NaN*ones(HH+1,1);
probR0D=NaN*ones(HH+1,1);

% Iterate for conditional probabilities of Doubters

for j=1:1:HH+1
    probSD(j,1)=S(j,1)/(S(j,1)+Ia0(j,1)+R0(j,1));
    probA0D(j,1)=Ia0(j,1)/(S(j,1)+Ia0(j,1)+R0(j,1));
    probR0D(j,1)=R0(j,1)/(S(j,1)+Ia0(j,1)+R0(j,1));
end

% DOUBTERS - Susceptible People (continued)
Us=NaN*ones(HH+1,1);
Usss=Ur1ss; %PV utility of susceptibles same as recovered in steady state
Us(HH+1)=Usss;%terminal condition
for tt=HH:-1:1    
    Us(tt,1)=us(tt)+probSD(tt).*betta*((1-((T(tt)*pib)/S(tt))).*Us(tt+1,1)+((T(tt)*pib)/S(tt)).*Uib0(tt+1,1))+probA0D(tt).*betta*((1-pit).*Us(tt+1,1)+pit*Uib0(tt+1,1))+probR0D(tt).*betta*Us(tt+1,1);
end

%Lagrange multipliers doubters - susceptibles
lamtaus=betta*pib*probSD(2:HH+1).*(Uib0(2:HH+1)-Us(2:HH+1));
lambs=(cs.^(-1)+lamtaus*pis1.*I(1:HH).*cib0)./(1+muc0);

% DOUBTERS - Asymptomatically Infected Untested (continued)
Uia0=NaN*ones(HH+1,1);
Uia0ss=Ur1ss; %PV utility of susceptibles same as recovered in steady state
Uia0(HH+1)=Uia0ss;%terminal condition
for tt=HH:-1:1    
    Uia0(tt,1)=uia0(tt)+probSD(tt).*betta*((1-((T(tt)*pib)/S(tt))).*Uia0(tt+1,1)+((T(tt)*pib)/S(tt)).*Uib0(tt+1,1))+probA0D(tt).*betta*((1-pit).*Uia0(tt+1,1)+pit*Uib0(tt+1,1))+probR0D(tt).*betta*Uia0(tt+1,1);
end

%Lagrange multipliers doubters - Asymptomatically Infected Untested
lamtauia0=betta*pib*probSD(2:HH+1).*(Uib0(2:HH+1)-Uia0(2:HH+1));
lambia0=(cia0.^(-1)+lamtauia0*pis1.*I(1:HH).*cib0)./(1+muc0);

% DOUBTERS - Recovered UNAWARE People (continued)
Ur0=NaN*ones(HH+1,1);
Ur0ss=Ur1ss; %PV utility of susceptibles same as recovered in steady state
Ur0(HH+1)=Ur0ss;%terminal condition
for tt=HH:-1:1    
    Ur0(tt,1)=ur0(tt)+probSD(tt).*betta*((1-((T(tt)*pib)/S(tt))).*Ur0(tt+1,1)+((T(tt)*pib)/S(tt)).*Uib0(tt+1,1))+probA0D(tt).*betta*((1-pit).*Ur0(tt+1,1)+pit*Uib0(tt+1,1))+probR0D(tt).*betta*Ur0(tt+1,1);
end

%Lagrange multipliers doubters - Recovered UNAWARE 
lamtaur0=betta*pib*probSD(2:HH+1).*(Uib0(2:HH+1)-Ur0(2:HH+1));
lambr0=(cr0.^(-1)+lamtaur0*pis1.*I(1:HH).*cib0)./(1+muc0);

%equation residuals - Gradient Adjusting Method Extended for SI4R.
err=ones(6*HH,1);
err(1:HH)=(1+muc0).*cia0-A.*nia0-Gamma;
err(HH+1:2*HH)=(1+muc1).*cia1-A.*nia1-Gamma;
err(2*HH+1:3*HH)=(1+muc0).*cib0-phii*A.*nib0-Gamma;
err(3*HH+1:4*HH)=(1+muc1).*cib1-phii*A.*nib1-Gamma;
err(4*HH+1:5*HH)=muc0.*(S(1:HH).*cs+Ia0(1:HH).*cia0+Ib0(1:HH).*cib0+R1(1:HH).*cr0+R1(1:HH).*cr1)+muc1.*(Ia1(1:HH).*cia1+Ib1(1:HH).*cib1)-Gamma.*(S(1:HH)+Ia0(1:HH)+Ia1(1:HH)+Ib0(1:HH)+Ib1(1:HH)+R0(1:HH)+R1(1:HH));
%err(4*HH+1:5*HH)=muc1.*(S(1:HH).*cs+(Ia0(1:HH).*cia0+Ia1(1:HH).*cia1+Ib0(1:HH).*cib0+Ib1(1:HH).*cib1)+R1(1:HH).*cr0+R0(1:HH).*cr1)-Gamma.*(S(1:HH)+Ia0(1:HH)+Ia1(1:HH)+Ib0(1:HH)+Ib1(1:HH)+R0(1:HH)+R1(1:HH));
err(5*HH+1:6*HH)=-theta*ns+A.*lambs+lamtaus*pis2.*I(1:HH).*nib1;

%Aggregate consumption and hours
aggC=S(1:HH).*cs+Ia0(1:HH).*cia0+Ia1(1:HH).*cia1+Ib0(1:HH).*cib0+Ib1(1:HH).*cib1+R0(1:HH).*cr0+R1(1:HH).*cr1;
aggH=S(1:HH).*ns+Ia0(1:HH).*nia0+Ia1(1:HH).*nia1+Ib0(1:HH).*nib0*phii+Ib1(1:HH).*nib1*phii+R0(1:HH).*nr0+R1(1:HH).*nr1;

%Present value of society utility 
U=S(1:HH).*Us(1:HH)+Ia0(1:HH).*Uia0(1:HH)+Ia1(1:HH).*Uia1(1:HH)+Ib0(1:HH).*Uib0(1:HH)+Ib1(1:HH).*Uib1(1:HH)+R0(1:HH).*Ur0(1:HH)+R1(1:HH).*Ur1(1:HH);

RnotSIRmacro=(Ta0(1)+Tb0(1))/(Ia0(1)+Ia1(1)+Ib0(1)+Ib1(1))/(pir+pid);

end