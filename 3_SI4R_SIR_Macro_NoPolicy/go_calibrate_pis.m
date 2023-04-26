%calibrate pis1, pis2 and pis3 using SIR model
scale1=1000000;scale2=1000;%scale pis for numerical solver
[sol,fval,exitflag]=fsolve(@calibrate_pis,[0.2;0.2;0.2],opts_fsolve,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,cr1ss,nr1ss,scale1,scale2,pia,pib,pit);

if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else
    [err,pis1,pis2,pis3,RnotSIR,S,D,Ta0,Tb0,T,Ia0,Ia1,Ib0,Ib1,I,R0,R1,R,Pop] =calibrate_pis(sol,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,cr1ss,nr1ss,scale1,scale2,pia,pib,pit);
    
    disp(['Max. abs. error in calibration targets:',num2str(max(abs(err)))]);
    disp([' ']);
    pis1=sol(1)/scale1
    pis2=sol(2)/scale2
    pis3=sol(3)
    RnotSIR
end
%  
ia=4;
ib=4;
figure;
subplot(ia,ib,1);
plot(S);axis tight;
title('S');
subplot(ia,ib,2);
plot(I);axis tight;
title('I');
subplot(ia,ib,3);
plot(D);axis tight;
title('D');
subplot(ia,ib,4);
plot(R);axis tight;
title('R');
subplot(ia,ib,5);
plot(T);axis tight;
title('T');
subplot(ia,ib,6);
plot(Ta0);axis tight;
title('Ta-');
subplot(ia,ib,7);
plot(Tb0);axis tight;
title('Tb-');
subplot(ia,ib,9);
plot(Ia0);axis tight;
title('Ia-');
subplot(ia,ib,10);
plot(Ia1);axis tight;
title('Ia+');
subplot(ia,ib,11);
plot(Ib0);axis tight;
title('Ib-');
subplot(ia,ib,12);
plot(Ib1);axis tight;
title('Ib+');
subplot(ia,ib,13);
plot(R0);axis tight;
title('R-');
subplot(ia,ib,14);
plot(R1);axis tight;
title('R+');

suptitle('SI4R Model, Calibration')

orient landscape
print -dpdf -fillpage SIR_calibration_pis