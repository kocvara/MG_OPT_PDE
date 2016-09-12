% README file for the Matlab software related to the paper 
%--------------------------------------------------------------------------
% M. Kocvara and S. Mohammed. 
% A first-order multigrid method for bound-constrained convex optimization. 
% Optimization Methods and Software 31.3 (2016): 622-644.
%
% Developed and coded by Michal Kocvara, m.kocvara@bham.ac.uk
% June 2016
% This is academic testing software coming with no guarantees!
%--------------------------------------------------------------------------
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% NOTICE that Example 4 requires an "external" quadratic programming solver
% called in file 'gp_pde1_con_eq.m'. At the moment, I use Gurobi, but the
% file includes uncommented calls to CVX-Gurobi or Matlab's own 'quadprog'.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% There are four "main" Matlab functions associated with the four examples
% in the paper called as:
%
% >> mg_pde_ex1(nx,ny,levels);
% >> mg_pde_ex2(nx,ny,levels);
% >> mg_msurf_ex3(nx,ny,levels);
% >> mg_pde_eq_ex4(nx,ny,levels);
%
% where 'nx, ny' are the numbers of elements in the x and y axis on the
% coarsest level and 'levels' is the number of refinement levels (at least
% 2). I typically use nx=ny=2 and levels=2..9 .
%
% Examples 1,2,3 use the truncation multigrid algorithm described in the
% paper, while Example 4 uses the Hackbusch-Mittelmann algorithm.
%
% In all main functions the currect solution is compared with the "exact"
% solution which has been precomputed with high precision and stored in the
% mat-files called 'xstar*'. Note that the example number in the name in
% these files does not correspond to the example number in the paper. The
% relation is 
% Example 1 (paper) -> xstar_ex3_lev*.mat
% Example 2 (paper) -> xstar_ex1_lev*.mat 
% Example 3 (paper) -> xstar_ex4_lev*.mat 
% Example 4 (paper) -> xstar_ex5_lev*.mat 
% If you don't want to use these precomputed solutions (for instance, if 
% you want to change the obstacle, you can use the (uncommented) call to
% alternative solvers, found in the main functions.
%
% NOTICE that I have made some changes in the codes after the publication
% of the paper, so the numbers you will obtain will not be exactly the same 
% as those reported in the paper. 

