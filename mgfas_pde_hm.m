function [xn,merit]=mgfas_pde_hm(levels,level,Aglob,fl,intnodel,cf,Pro,vh,x,xl,xu,iiii)
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
% Solving PDE1 and 2 from Wen-Goldfarb
% Hackbusch-Mittelmann algorithm

global rho gp_count xstar  nnodesg nxxx nyyy
if level==levels, rho=min(2*rho,1e9); end

gpiter=1;

Aco = Aglob{level};Aco=Aco(intnodel{level},intnodel{level});
ff=fl{level};

l=level;
if level==1
%% level-one exact solution
    
     [xn,k]=gp_pde1_con(x,Aco,ff,intnodel{1},vh,10000,1e-9,xl,xu);
        
else
    PP=Pro{l}';
    P = PP(intnodel{l},intnodel{l-1});
%% Pre-smoothing
        kk=0;
        [x,kk]=gp_pde1_con(x,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level);

    if level==levels, gp_count=gp_count+kk;end
    
    xH=0.25.*P'*x;
    %xfull = zeros(nnodesg{l},1); xfull(intnodel{l})= x;
    %xHfull = xfull(cf{level}); xH = xHfull(intnodel{l-1});
    
    [dum,grad] = pde1_vh(x,Aco,ff,intnodel{level},vh);
    vH1=P'*grad;
    
    AcoH = Aglob{level-1};AcoH=AcoH(intnodel{level-1},intnodel{level-1});
    ffH=P'*ff;
    [dum,vH2] = pde1_vh(xH,AcoH,fl{level-1},intnodel{level-1},zeros(size(xH)));
    
    
    % RHS for the coarse level
    vH = vH2 - vH1;
    
    x1=xH;
    xPPxl = xl - x + 1.*P*x1;
    xPPxu = xu - x + 1.*P*x1;
    xl1=zeros(length(x1),1);
    xu1=zeros(length(x1),1);
    for i=1:length(x1);
        [inbgh,idum,vdum] = find(P(:,i));
        xl1(i) = min(xPPxl(inbgh));
        xu1(i) = max(xPPxu(inbgh));
    end
    
    xl1=0.25.*P'*xl; xu1=0.25.*P'*xu;
    
    for i=1:length(x1);
        [inbgh,idum,vdum] = find(P(:,i));
        if max(xl(inbgh)-x(inbgh))==0
            xl1(i,1) = x1(i); %xu1(i,1) = x1(i);
        end
        if min(xu(inbgh)-x(inbgh))==0
            xu1(i,1) = x1(i); %xl1(i,1) = x1(i);
        end
    end
    if min(xu1-xl1)<0, max(xl1),min(xu1), end
    
    [cor]=mgfas_pde_hm(levels,l-1,Aglob,fl,intnodel,cf,Pro,vH,xH,xl1,xu1,iiii);
    
    xcor = cor-xH;
    
    if iiii==0
        nxl = nxxx*2^(level-2)+1; nyl = nyyy*2^(level-2)+1;
        a=zeros(nxl,1);b=zeros(nyl,1);c=zeros(nxl,1);d=zeros(nyl,1);
        xx2=zeros(nxl,nyl);
        xx2(:,1)=a; xx2(end,:)=b; xx2(:,end)=c; xx2(1,:)=d;
        xx2(2:end-1,2:end-1) = reshape(xcor,nyl-2,nxl-2);
        xx=reshape(xx2,nyl*nxl,1);
        xcorpf = PP*xx; xcorp = xcorpf(intnodel{l});
    else
        xcorp = 1.*P*xcor;
    end
    
    xc = x+xcorp;
    
%% Post-smoothing
    
    [xn,kk]=gp_pde1_con(xc,Aco,ff,intnodel{level},vh,gpiter,1e-11,xl,xu,level);
    
    if level==levels, gp_count=gp_count+kk;end
    merit=0;
    
end

