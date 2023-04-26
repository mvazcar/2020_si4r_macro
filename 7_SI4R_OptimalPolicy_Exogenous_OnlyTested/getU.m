function [U0,err,I,S,R,D,T,Pop,cs,ns,Us,Rnot,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U]=getU(muc1_guess,n_vec_guess,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,phii,pia,pib,pit)

%append guess with zeros so that muc1 has length HH (in case you want to
%find optimal muc1 with horizon less than HH)
muc1=zeros(HH,1);
muc1(1:numel(muc1_guess))=muc1_guess;
    
%create persistent variable, i.e. dont start
%solution of nonlinear equilibrium equations from scratch
%every time in the optimization of optimal containment policy

persistent guess_n n_vec
if isempty(guess_n)==1
    guess_n=n_vec_guess;
else
    guess_n=n_vec;
end;


%get allocations for ns,nia0,nia1,nib0,nib1,nr0,nr1 given muc1

try
    [n_vec,fval,exitflag]=fsolve(@get_err,guess_n,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,muc0,muc1,phii,pia,pib,pit);

    if exitflag==1
        display('Equation solved. First Order Optimality is Small')
    elseif exitflag==2
        display('Equation solved. The change is smaller than the tolerance specified')
    elseif exitflag==3
        display('Equation solved. The change in the residuals is less than specified tolerance')
    elseif exitflag==4
        display('Equation solved. The magnitude of the search direction is less than the specified tolerance')
    else
        error('Fsolve could not solve the model')
    end    

    [err,I,Ia0,Ia1,Ib0,Ib1,S,R,R0,R1,D,T,Pop,RnotSIRmacro,aggC,aggH,ns,nia0,nia1,nib0,nib1,nr0,nr1,cs,cia0,cia1,cib0,cib1,cr0,cr1,Us,Uia0,Uia1,Uib0,Uib1,Ur0,Ur1,U,probSD,probA0D,probR0D] = get_err(n_vec,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,muc0,muc1,phii,pia,pib,pit);

catch
    U=-1e10;
    guess_n=[];
end

%PV Utility at time, t=0
%Use minimizer (fmincon) to find optimal policy muc1 such that
%U1 is maximal -- hence flip sign
U0=-U(1);
end




