function [x,gamma_s]=gp_pde1_con_eq(x0,A,ff,intnode,vh,minit,tol,xl,xu,level,Area)
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
%gradient projection method for box constrained PDE 1 and 2

% projected gradient method with Armijo back-tracking line-search

global gamma_s_old

ooo = optimset('LargeScale','on','TolCon',1e-13,'TolFun',1e-13,'Display','off',...
    'MaxIter',1000,'Algorithm','interior-point-convex');
k=0;
x0=x0(:);
x=x0;
m = length(x);
Ar=Area;  qAeq = Ar'; qbeq =sum(Area.*x0); 
qH = speye(m);

%%
vol=1;
[obj0,grad] = pde1_vh_Ex4(x,A,ff,intnode,vh); obj = obj0; grad0=grad;
diffc = Inf;
while ( diffc > tol && k < minit)
    gamma = 0.5; gamma_s = 1;
    aaa=1e40;
    bbb=-1e40;
    iiter=1;
    while ( aaa > bbb && iiter<20)
        xn = x - gamma_s*grad;
        qf = -xn;
        
%% QUADPROG projection        
        %         [xn,fval,exitflag,outputqp,lambda] = ...
        %             quadprog(qH,qf,[],[],qAeq,qbeq,xl,xu,xn,ooo);
        

%% CVX projection       
%                 cvx_begin
%                 cvx_solver gurobi
%                 cvx_quiet true
%                 cvx_precision([1e-13])
%                 variable xnn(size(xn))
%                 minimize( 0.5*xnn'*qH*xnn + xnn'*qf )
%                 minimize( norm(xn-xnn) )
%                 subject to
%                 xnn >= xl
%                 xnn <= xu
%                 qAeq*xnn == qbeq
%                 cvx_end
%                 xn1=xnn;

%% GUROBI projection
        clear model;
        model.Q = qH;
        model.A = sparse(qAeq);
        model.obj = 2.*qf;
        model.rhs = qbeq;
        model.sense = '=';
        model.lb = xl;
        model.ub = xu;
        params.OutputFlag = 0;
        params.FeasibilityTol = 1e-9;
        params.BarConvTol  = 1e-13;  
%        gurobi_write(model, 'qp.lp');
        results = gurobi(model,params);  
        xn=results.x;

 %%
        [obj1,grad1] = pde1_vh_Ex4(xn,A,ff,intnode,vh);
        aaa = obj1;
        bbb = obj -.1*grad'*(x-xn);
        gamma_s = gamma_s*gamma;
        iiter = iiter+1;
    end
    
    if iiter<20
        x=xn; grad=grad1;
        diffc = abs(obj1-obj);
        obj=obj1;
    else
        diffc = 0;
        fprintf(' line search failure\n');
    end
    k=k+1;
end

