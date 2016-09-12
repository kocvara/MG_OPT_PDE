function [Avec]=Amatq(levels,ivg,vxy)
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
% "stiffness matrix" for the minimum surface problem 
% rectangular domain, square bilinear finite elements

Alocloc = sparse([4 -1 -2 -1;-1 4 -1 -2;-2 -1 4 -1;-1 -2 -1 4])./6;

for level=1:levels
    nnod =  size(vxy{level},2);
    nelem =  size(ivg{level},2);
    Avec1 = []; 
    h = 2^-(level-1);
    for ie = 1:nelem;
        % node numbers for triangle K
        nodes = ivg{level}(:,ie);
        % element stiffness matrix
        Aloc=sparse(nnod,nnod);
        Aloc(nodes,nodes) = Alocloc;
        Avec1 = [Avec1 Aloc];
    end
    Avec{level} = Avec1';
end

