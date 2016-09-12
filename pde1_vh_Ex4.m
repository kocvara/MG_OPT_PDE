function [f,g] = pde1_vh(x,AA,b,intnode,vh)
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
% objective function and gradient for Example 4 with 
% additional linear term -vh'*x in the objective function
%
global fcount xstar

% checking if we are on the top level
if length(x)==length(xstar), fcount = fcount+1; end
% fcount = fcount+1;

h=1/length(x); lambda=1.0;

%% EXAMPLE 4
  f = 0.5*x'*AA*x - 0.3333333*h*sum(x.^3) - (b+vh)'*x;
  g = AA*x-h.*x.^2-b;  g=g-vh;

end

