function [f,g] = msurfq_vh(x,AA,intnode,vh,level)
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
% objective function and gradient for the 2D minimal surface problem with 
% additional linear term -vh'*x in the objective function
%
global fcount xstar nxxx nyyy
if length(x)==length(xstar), fcount = fcount+1; end


%nel = size(AA,2);
nxl = nxxx*2^(level-1)+1; nyl = nyyy*2^(level-1)+1;
nnod = nxl*nyl; nel=size(AA,1)/nnod;
area = 1/(nxl-1)^2;area=1;
a=zeros(nxl,1);b=zeros(nyl,1);c=zeros(nxl,1);d=zeros(nyl,1);
xxx=0:1/(nxl-1):1;

a=-sin(2*pi*xxx);b=-a;c=-a;d=a;

xx2=zeros(nxl,nyl);
xx2(:,1)=a; xx2(end,:)=b; xx2(:,end)=c; xx2(1,:)=d;
xx2(2:end-1,2:end-1) = reshape(x,nyl-2,nxl-2);
xx=reshape(xx2,nyl*nxl,1);

ss = zeros(nel,1); gg=sparse(nnod,nel); Alocxx=sparse(nnod,nel);

    AAA = sparse(AA*xx);
    Alocxx = reshape(AAA',nnod,nel);
    ssl = 1./sqrt(1+xx'*Alocxx);
    gg = Alocxx*spdiags(ssl',0,nel,nel);
    ss=1./ssl;

f = area*sum(ss) - vh'*x;
g1 = area.*sum(gg,2);
g = g1(intnode);
g=g-vh;
end
