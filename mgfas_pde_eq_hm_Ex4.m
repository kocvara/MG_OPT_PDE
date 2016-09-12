function [xn,merit]=mgfas_pde_hm_Ex4(levels,level,Aglob,fl,intnodel,cf,Pro,vh,x,xl,xu,iiii,Area)
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
% Solving Example 4 with an equality constraint by Hackbusch-Mittelmann 
% algorithm

global rho gp_count xstar  nnodesg nxxx nyyy
if level==levels, rho=min(2*rho,1e9); end

gpiter=1;

Aco = Aglob{level};Aco=Aco(intnodel{level},intnodel{level});
ff=fl{level};
Area_v=Area{level};

l=level;
if level==1
%% level-one exact solution
    if iiii==0
        vh=zeros(size(vh));
    end
    [xn,funcevals]=gp_pde1_con_eq(x,Aco,ff,intnodel{level},vh,1000,1e-11,xl,xu,level,Area_v);
else
    PP=Pro{l}';
    P = PP(intnodel{l},intnodel{l-1});
%% Pre-smoothing
        kk=0;
        [x,kk]=gp_pde1_con_eq(x,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level,Area_v);

    if level==levels, gp_count=gp_count+kk;end
    
    xH=0.25.*P'*x;
    
    [dum,grad] = pde1_vh_Ex4(x,Aco,ff,intnodel{level},vh);
    vH1=P'*grad;
    
    AcoH = Aglob{level-1};AcoH=AcoH(intnodel{level-1},intnodel{level-1});
    ffH=P'*ff;
    [dum,vH2] = pde1_vh_Ex4(xH,AcoH,fl{level-1},intnodel{level-1},zeros(size(xH)));  
    
    % RHS for the coarse level
    vH = vH2 - vH1;
    
    x1=xH;
    xPPxl = xl - x + 1.*P*x1;
    xPPxu = xu - x + 1.*P*x1;
    xl1=zeros(length(x1),1);
    xu1=zeros(length(x1),1);
    for i=1:length(x1);
        [inbgh,idum,vdum] = find(P(:,i));
        xl1(i) = max(xPPxl(inbgh));
        xu1(i) = min(xPPxu(inbgh));
    end
    
   xl1=0.25.*P'*xl; xu1=0.25.*P'*xu;
    
    for i=1:length(x1);
        [inbgh,idum,vdum] = find(P(:,i));
        if abs(max(xl(inbgh)-x(inbgh)))<1e-6
            xl1(i,1) = x1(i);xu1(i,1) = x1(i);
        end
        if abs(min(xu(inbgh)-x(inbgh)))<1e-6
            xu1(i,1) = x1(i);xl1(i,1) = x1(i);
        end
    end
    if min(xu1-xl1)<0, max(xl1),min(xu1), end
    
    [cor]=mgfas_pde_eq_hm_Ex4(levels,l-1,Aglob,fl,intnodel,cf,Pro,vH,xH,xl1,xu1,iiii,Area);
    
    xcor = cor-xH; 
    xcorp = 1.*P*xcor; 
    xc = x+xcorp;
    
    %% Post-smoothing
    
    [xn,kk]=gp_pde1_con_eq(xc,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level,Area_v);
    
    if level==levels, gp_count=gp_count+kk;end
    merit=0;
    
end

