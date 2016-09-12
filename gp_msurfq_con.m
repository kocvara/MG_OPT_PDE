function [x,k]=gp_membr(x0,A,intnode,vh,minit,tol,xl,xu,level)
%--------------------------------------------------------------------------
% Matlab software related to the paper 
%
% M. Kocvara and S. Mohammed. 
% A first-order multigrid method for bound-constrained convex optimization. 
% Optimization Methods and Software 31.3 (2016): 622-644.
%
% For the use please refer to the README file in this directory
%
% Developed and coded by Michal Kocvara, m.kocvara@bham.ac.uk
% June 2016
% This is academic testing software coming with no guarantees!
%--------------------------------------------------------------------------
%
% gradient projection method for box constrained minimal surface problem

global gamma_s_old

k=0;
x0=x0(:);
x=x0;
x = max(xl,x);
x = min(xu,x);

[obj0,grad] = msurfq_vh(x,A,intnode,vh,level); obj = obj0;
compl = 1e20;
gradact=grad; gradact(x==xl)=0;gradact(x==xu)=0;

while ( compl > tol && k < minit)
    gamma = 0.5; gamma_s = gamma_s_old; gaga = 2; 
    bbb1 = 0.49*grad; aaa=1e40; bbb=-1e40; iiter=1; 
    
    xn = x - gamma_s*grad; xn = max(xl,xn); xn = min(xu,xn);
    [objn,gradn] = msurfq_vh(xn,A,intnode,vh,level);
    gradnact=gradn; gradnact(xn==xl)=0;gradnact(xn==xu)=0;
    gg0 = gradnact'*grad;gg=gg0;
    if gg > 0
        while ( gg > 0 && iiter<20)
            gamma_s = gamma_s*gaga;
            xn = x - gamma_s*grad; xn = max(xl,xn); xn = min(xu,xn);
            [objn,gradn] = msurfq_vh(xn,A,intnode,vh,level);
            gradnact=gradn; gradnact(xn==xl)=0;gradnact(xn==xu)=0;
            gg=gradnact'*grad;
            gradn1=gradn;
            iiter = iiter+1;
        end
    else
        while ( gg < 0 && iiter<20)
            gamma_s = gamma_s/gaga;
            xn = x - gamma_s*grad; xn = max(xl,xn); xn = min(xu,xn);
            [objn,gradn] = msurfq_vh(xn,A,intnode,vh,level);
            gradnact=gradn; gradnact(xn==xl)=0;gradnact(xn==xu)=0;
            gg=gradnact'*grad;    
            iiter = iiter+1;
       end
    end
        if gg0>0
            x=x - (gamma_s/(gaga))*grad;
            gamma_s_old=(gamma_s/(gaga));
            grad=gradn1;
        else
            x=x - (gamma_s*(1))*grad;
            gamma_s_old=(gamma_s/(1));
            grad=gradn;
        end
    x = max(xl,x); x = min(xu,x);
    gradact=grad; gradact(x==xl)=0; gradact(x==xu)=0;
    compl = norm(grad.*(x-xl).*(xu-x));

    k=k+1;
end


