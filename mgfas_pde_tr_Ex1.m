function [xn]=mgfas_pde_tr(levels,level,Aglob,fl,intnodel,cf,Pro,vh,x,xl,xu,iiii)
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
% Truncated multigrid for Example 1

global rho gp_count xstar  nnodesg nxxx nyyy
if level==levels, rho=min(2*rho,1e9); end

gpiter=1;

if level==levels
    Aco = Aglob{level};Aco=Aco(intnodel{level},intnodel{level});
    ff=fl{level};
else
    Aco=Aglob;
    ff=fl;
end

l=level;
if level==1
%% level-one exact solution
    
     [xn,k]=gp_pde1_con_Ex1(x,Aco,ff,intnodel{1},vh,10000,1e-9,xl,xu);
        
else
    PP=Pro{l}';
    P = PP(intnodel{l},intnodel{l-1});
%% pre-smoothing
    [x,kk]=gp_pde1_con_Ex1(x,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level);
        
    if level==levels, gp_count=gp_count+kk; end
    
%% truncation
    ndotl = find(x==xl); ndotu = find(x==xu);
    ndot = union(ndotl,ndotu);
    if level==levels
        nF = length(intnodel{l});
        eediag = ones(nF,1); eediag(ndot)=0;
        eeF=sparse(1:nF,1:nF,eediag);
        Acot=eeF*Aco*eeF;
        AcoH=P'*Acot*P;
    else
        AcoH=P'*Aco*P;
    end
    
%%
    xH=0.25.*P'*x;
    xfull = zeros(nnodesg{l},1); xfull(intnodel{l})= x;
    xHfull = xfull(cf{level}); xH = xHfull(intnodel{l-1});
    
    [dum,grad] = pde1_vh_Ex1(x,Aco,ff,intnodel{level},vh);
    grad(ndot)=0;
    vH1=P'*grad;
    
    ffH=P'*ff;
    %xfull = zeros(nnodesg{l},1); xfull(intnodel{l})= ff;
    %xHfull = xfull(cf{level}); ffH = xHfull(intnodel{l-1});
    [dum,vH2] = pde1_vh_Ex1(xH,AcoH,ffH,intnodel{level-1},zeros(size(xH)));
    
    % RHS for the coarse level
    vH = vH2 - vH1;
    
    x1=xH;
    xPPxl = xl - x ;
    xPPxu = xu - x ;
    if level==levels
        xPPxl(xl==x)=-100;
        xPPxu(xu==x)= 100;
    end
    
    xl1=zeros(length(x1),1);
    xu1=zeros(length(x1),1);
    for i=1:length(x1);
        [inbgh,idum,vdum] = find(P(:,i));
        xl1(i) = max(xPPxl(inbgh));
        xu1(i) = min(xPPxu(inbgh));
    end
    xl1 = xl1 + x1; xu1 = xu1 + x1;

    [cor]=mgfas_pde_tr_Ex1(levels,l-1,AcoH,ffH,intnodel,cf,Pro,vH,xH,xl1,xu1,iiii);
    
    xcor = cor-xH;
    xcorF = P*xcor;
    xc = x+xcorF;

%% post-smoothing
    
    [xn,kk]=gp_pde1_con_Ex1(xc,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level);
    
    if level==levels, gp_count=gp_count+kk;end
    
end

