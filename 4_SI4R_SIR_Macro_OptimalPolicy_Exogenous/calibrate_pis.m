function [err,pis1,pis2,pis3,RnotSIR,S,D,Ta0,Tb0,T,Ia0,Ia1,Ib0,Ib1,I,R0,R1,R,Pop] =calibrate_pis(pis_guess,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,C,N,scale1,scale2,pia,pib,pit)

%back out initial guesses
pis1=pis_guess(1)/scale1;
pis2=pis_guess(2)/scale2;
pis3=pis_guess(3);

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
    Ta0(j,1)=pia*(pis1*S(j)*C^2*I(j)+pis2*S(j)*N^2*I(j)+pis3*S(j)*I(j));
    Tb0(j,1)=pib*(pis1*S(j)*C^2*I(j)+pis2*S(j)*N^2*I(j)+pis3*S(j)*I(j));
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

%calculate equation residuals for target equations
err(1)=pis1_shr_target-(pis1*C^2)/(pis1*C^2+pis2*N^2+pis3);
err(2)=pis2_shr_target-(pis2*N^2)/(pis1*C^2+pis2*N^2+pis3);
err(3)=RplusD_target-(R(end)+D(end));


RnotSIR=(Ta0(1)+Tb0(1))/(Ia0(1)+Ia1(1)+Ib0(1)+Ib1(1))/(pir+pid);

end