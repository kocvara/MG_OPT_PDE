function [xn,merit]=mgfas_msurfq_tr(levels,level,Avec,intnodel,cf,elc,Pro,vh,x,xl,xu,iiii)
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
% Solving the minimum surface with obstacle on a rectangular domain by FAS 
% multigrid with truncation

global rho gp_count nnodesg
if level==levels, rho=min(2*rho,1e9); end

gpiter=1;

if level==levels
    Avec_h = Avec{level};
else
    Avec_h=Avec;
end

l=level;
if level==1
%% level-one exact solution
    
    %opts.verbose=0; opts.optTol=1e-15;
    %[xn, dumdum, funcevals] = minConf_TMP(@(x) msurfq_vh(x,Avec_h,intnodel{1},vh,1),x,xl,xu,opts);
    
    [xn,kk]=gp_msurfq_con(x,Avec_h,intnodel{level},vh,10000,1e-9,xl,xu,level);
else
    PP=Pro{l}';
    P = PP(intnodel{l},intnodel{l-1});
%% Pre-smoothing
    
    [x,kk]=gp_msurfq_con(x,Avec_h,intnodel{level},vh,gpiter,1e-11,xl,xu,level);
    
    if level==levels, gp_count=gp_count+kk;end
    
%% truncation
    ndotl = find(x==xl); ndotu = find(x==xu);
    ndot = union(ndotl,ndotu);
    ndotfull=intnodel{l}(ndot);
    if level==levels
        nF = nnodesg{l}; nFF = nnodesg{l-1};
        eediag = ones(nF,1); eediag(ndotfull)=0;
        eeF=sparse(1:nF,1:nF,eediag);
        for ielH=1:size(Avec_h,1)/nF
            Ako=eeF'*Avec_h(nF*(ielH-1)+1:nF*(ielH-1)+nF,:)*eeF;
            Avec_hr(nF*(ielH-1)+1:nF*(ielH-1)+nF,:)=Ako;
        end
        for ielH=1:size(elc{level},2)
            elH=elc{level}(:,ielH);
            Ako=PP'*Avec_hr(nF*(elH(1)-1)+1:nF*(elH(1)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_hr(nF*(elH(2)-1)+1:nF*(elH(2)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_hr(nF*(elH(3)-1)+1:nF*(elH(3)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_hr(nF*(elH(4)-1)+1:nF*(elH(4)-1)+nF,:)*PP;
            Avec_H(nFF*(ielH-1)+1:nFF*(ielH-1)+nFF,:) = Ako;
        end
    else
        nF = nnodesg{l}; nFF = nnodesg{l-1};
        for ielH=1:size(elc{level},2)
            elH=elc{level}(:,ielH);
            Ako=PP'*Avec_h(nF*(elH(1)-1)+1:nF*(elH(1)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_h(nF*(elH(2)-1)+1:nF*(elH(2)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_h(nF*(elH(3)-1)+1:nF*(elH(3)-1)+nF,:)*PP;
            Ako=Ako+PP'*Avec_h(nF*(elH(4)-1)+1:nF*(elH(4)-1)+nF,:)*PP;
            Avec_H(nFF*(ielH-1)+1:nFF*(ielH-1)+nFF,:) = Ako;
        end
    end
    
%%
    xH=0.25.*P'*x;
    %xfull = zeros(nnodesg{l},1); xfull(intnodel{l})= x;
    %xHfull = xfull(cf{level}); xH = xHfull(intnodel{l-1});
    
    [dum,grad] = msurfq_vh(x,Avec_h,intnodel{level},vh,level);
    grad(ndot)=0;
    vH1=1.*P'*grad;
    
    [dum,vH2] = msurfq_vh(xH,Avec_H,intnodel{level-1},zeros(size(xH)),level-1);
    
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
    
    [cor]=mgfas_msurfq_tr(levels,l-1,Avec_H,intnodel,cf,elc,Pro,vH,xH,xl1,xu1,iiii);
    
    xcor = cor-xH;
    
    xcorp = P*xcor;
    
    xc = x+xcorp;
    kokos=1;
    
%% Post-smoothing
    
    [xn,kk]=gp_msurfq_con(xc,Avec_h,intnodel{level},vh,gpiter,1e-11,xl,xu,level);
    
    if level==levels, gp_count=gp_count+kk;end
    merit=0;
    
end

