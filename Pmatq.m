function [Pro]=Pmatq(levels,ivg,vxy,cf,nex,ney,intnodel)
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
% prolongation operator for square elements on a square domain,
% regular mesh, regular refinement

for lev=2:levels
    %lev
    intno = intnodel{lev-1};
    nexl = nex*2^(lev-1); nx = nexl+1;
    neyl = ney*2^(lev-1); ny = neyl+1;
    nexlc = nex*2^(lev-2); nxc = nexlc+1;
    neylc = ney*2^(lev-2); nyc = neylc+1;
    %nnc = size(vxy{lev-1},2);
    pmat=sparse(nxc*nyc,nx*ny);
    
    ic = 1;
    i = cf{lev}(ic);
    pmat(ic,i) = 1;
    pmat(ic,i+1) = 0.5;
    pmat(ic,i+ny) = 0.5;
    pmat(ic,i+ny+1) = 0.25;
    
    ic = nyc;
    i = cf{lev}(ic);
    pmat(ic,i) = 1;
    pmat(ic,i-1) = 0.5;
    pmat(ic,i+ny) = 0.5;
    pmat(ic,i+ny-1) = 0.25;
    
    ic = nyc*(nxc-1)+1;
    i = cf{lev}(ic);
    pmat(ic,i) = 1;
    pmat(ic,i+1) = 0.5;
    pmat(ic,i-ny) = 0.5;
    pmat(ic,i-ny+1) = 0.25;
    
    ic = nxc*nyc;
    i = cf{lev}(ic);
    pmat(ic,i) = 1;
    pmat(ic,i-1) = 0.5;
    pmat(ic,i-ny) = 0.5;
    pmat(ic,i-ny-1) = 0.25;
    
    for ic=2:nyc-1
        i = cf{lev}(ic);
        pmat(ic,i) = 1;
        pmat(ic,i-1) = 0.5;
        pmat(ic,i+1) = 0.5;
        pmat(ic,i+ny) = 0.5;
        pmat(ic,i+ny-1) = 0.25;
        pmat(ic,i+ny+1) = 0.25;
    end
    
    for ic=nyc*(nxc-1)+2:nxc*nyc-1
        i = cf{lev}(ic);
        pmat(ic,i) = 1;
        pmat(ic,i-1) = 0.5;
        pmat(ic,i+1) = 0.5;
        pmat(ic,i-ny) = 0.5;
        pmat(ic,i-ny-1) = 0.25;
        pmat(ic,i-ny+1) = 0.25;
    end
    
    for ic=nyc+1:nyc:(nxc-1)*(nyc-1)
        i = cf{lev}(ic);
        pmat(ic,i) = 1;
        pmat(ic,i+1) = 0.5;
        pmat(ic,i-ny) = 0.5;
        pmat(ic,i+ny) = 0.5;
        pmat(ic,i-ny+1) = 0.25;
        pmat(ic,i+ny+1) = 0.25;
    end
    
    for ic=2*nyc:nyc:nxc*(nyc-1)
        i = cf{lev}(ic);
        pmat(ic,i) = 1;
        pmat(ic,i-1) = 0.5;
        pmat(ic,i-ny) = 0.5;
        pmat(ic,i+ny) = 0.5;
        pmat(ic,i-ny-1) = 0.25;
        pmat(ic,i+ny-1) = 0.25;
    end
    
    for iic=1:length(intno)
        ic = intno(iic);
        i = cf{lev}(ic);
        pmat(ic,i) = 1;
        pmat(ic,i-1) = 0.5;
        pmat(ic,i+1) = 0.5;
        pmat(ic,i-ny) = 0.5;
        pmat(ic,i+ny) = 0.5;
        pmat(ic,i-ny-1) = 0.25;
        pmat(ic,i-ny+1) = 0.25;
        pmat(ic,i+ny-1) = 0.25;
        pmat(ic,i+ny+1) = 0.25;
    end
    Pro{lev} = pmat;
    
end

