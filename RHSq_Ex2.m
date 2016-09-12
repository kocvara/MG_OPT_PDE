function [fl,intnodel]=RHSq(levels,ivg,vxy,nx,ny)
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
% Right-hand side and vector of internal nodes for Example 2
% rectangular domain descretized by square finite elements

fl=cell(1,levels);
intnodel=cell(1,levels);

for level=1:levels
    %clear p t tt Kbig Kvec
    
    %% Definition of volume and boundary forces

    % EXAMPLE 2
    f1=inline('0.01.*(9*pi^2+10*exp((x.^2-x.^3).*sin(3*pi.*y)).*(x.^2-x.^3)+6.*x-2)*sin(3*pi.*x)');

    %% Volume forces
    
    nxl = nx*2^(level-1)+1; nyl = ny*2^(level-1)+1;
    n=nxl*nyl;
    h = 1/(nxl-1);
    f = ones(n,1).*h*h;
    
    xxx=0:1/(nxl-1):1;  yyy=0:1/(nyl-1):1;  [X,Y]=meshgrid(xxx,yyy');
    fdum=f1(X,Y); f=fdum(:).*h*h;
    
    %% Internal nodes
    nx=nx;ny=ny;
    nxl = nx*2^(level-1)+1; nyl = ny*2^(level-1)+1;
    intno=[]; lum=0;
    for i=1:nxl-2
        intno = [intno [nyl+2:nyl+nyl-1]+lum];
        lum = lum+nyl;
    end
    
    intnodel{level} = intno;
    fl{level}=f(intno);
    
end

