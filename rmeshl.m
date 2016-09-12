
function [ivg,vxy,cf,elc]=rmeshl(nex,ney,levels)
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
% defining regular discretization of a rectangle 

for lev=1:levels
    nexl = nex*2^(lev-1); neyl = ney*2^(lev-1);
    nx = nexl+1; ny = neyl+1;
    nnod =  nx*ny;
    nelem =  (nx-1)*(ny-1);
    hx = 1; hy = 1;
    x = 0:hx:nx-1; y = 0:hy:ny-1;
    
    vxyl = zeros(2,nnod);
    adu = repmat(x,ny,1); vxyl(1,:) =adu(:)';
    vxyl(2,:) = repmat(y,1,nx);
    
    [ivgl] = rmesh(nx,ny,x,y);
    ivg{lev} = ivgl;
    vxy{lev} = vxyl;
    
    if lev>1
        nxc = (nx-1)/2+1; nyc = (ny-1)/2+1;
        cfl = zeros(nxc*nyc,1);
        inc = 1; infi = 1;
        for ix=1:nxc
            for iy=1:nyc
                cfl(inc)  =  infi;
                inc    =  inc + 1;
                infi   =  infi + 2;
            end
            infi    =  infi + ny-1;
        end
        cf{lev} = cfl;
    end
    
    if lev>1
        nexlH = nex*2^(lev-2); neylH = ney*2^(lev-2);
        elcl = zeros(4,nexlH*neylH);
        ielH = 1; iel = 1;
        for ix=1:nexlH
            for iy=1:neylH
                elcl(1,ielH)  =  iel;
                elcl(2,ielH)  =  iel+1;
                elcl(3,ielH)  =  iel+neyl;
                elcl(4,ielH)  =  iel+neyl+1;
                iel = iel + 2;
                ielH = ielH + 1;
            end
            iel = iel + neyl;
        end
        elc{lev} = elcl;
    end
    
end

end

function [ivg,vxy]=rmesh(nx,ny,x,y)

in    =  1;
nnod =  nx*ny;
nelem =  (nx-1)*(ny-1);

ivg = zeros(4,nelem);

ie=1;
for ix=1:nx-1
    for iy=1:ny-1
        iv(1,1)  =  in;
        iv(4,1)  =  in + 1;
        iv(3,1)  =  iv(4) + ny;
        iv(2,1)  =  iv(1) + ny;
        ivg(:,ie) = iv;
        in    =  in + 1;
        ie = ie+1;
    end
    in    =  in + 1;
end

end