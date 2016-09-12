function [x]=mg_pde_ex2(nx,ny,levels)
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
% Solving Example 2
% Discretization by square bilinear finite elements
% nx,ny ... number of elements in x and y direction

global rho gp_count xstar fcount nnodesg nxxx nyyy levlev
global gamma_s_old

gamma_s_old=1;
rho=100; gp_count=0; fcount=0; nxxx=nx;nyyy=ny; levlev=levels;

 [ivg,vxy,cf]=rmeshl(nx,ny,levels);

load Amat_file

 [ivg,vxy,cf]=rmeshl(nx,ny,levels);
 [fl,intnodel]=RHSq_Ex2(levels,ivg,vxy,nx,ny);

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;

for i=2:levels
    nnodesg{i} = size(vxy{i},2);
end

level=levels;

n = length(intnodel{levels});
x=zeros(n,1);
vh=zeros(n,1);

%% Definition of the obstacle

% EXAMPLE 2: parabolic lower bound
aup=0.5; xu(1:n,1)=aup;
xxx=0:1/(nxf-1):1;yyy=0:1/(nyf-1):1;[X,Y]=meshgrid(xxx,yyy');
XLXL=-8.*(X-7/16).^2-8.*(Y-7/16).^2+.2;    
xlxl=XLXL(:); xl=xlxl(intnodel{levels});


%% Exact solution
Aco = Aglob{levels};Aco=Aco(intnodel{levels},intnodel{levels});

%opts    = struct( 'factr', 1e-1, 'pgtol', 1e-19, 'm', 9, 'maxIts',10000000,'maxTotalIts',100000000);
%opts.printEvery     = 100;
%[xstar, dumdum, info] = lbfgsb(@(x) pde1_vh(x,Aco,fl{levels},intnodel{levels},zeros(size(x))), xl, xu, opts);

% [xstar,kk]=gp_pde1_con(x,Aco,fl{levels},intnodel{levels},zeros(size(x)),10000000,1e-11,xl,xu,levels);toc

% reading exact solution from a file
% EXAMPLE 1 = ex3; EXAMPLE 2 = ex1;
name = strcat('xstar_ex1_lev',num2str(levels));
load(name)
funcevals

error(1) = sqrt(sum((x-xstar).^2)/length(x));

tic
%% Multigrid iterations
g=ones(size(x));

niter = 30;
eee=1;i=1;ndot=[];
while eee>2e-6 & i<niter
    vh = zeros(length(x),1); %vh=fl{levels};
    
  % truncated multigrid  
    [x]=mgfas_pde_tr_Ex2(levels,level,Aglob,fl,intnodel,cf,Pro,vh,x,xl,xu,i);

    error(i+1) = sqrt(sum((x-xstar).^2)/length(x));
    
    eee = norm(error(i+1));
    i=i+1;
%
    fprintf('convergence speed, error: %3d %6f %6f\n',i,exp((log(error(i))-log(error(i-1)))),eee);
end
nitt=i-1;

%%
toc

gp_count
fcount

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;
xx2=zeros(nxf,nyf);
xx2(2:end-1,2:end-1) = reshape(x,nyf-2,nxf-2);
surf(xx2);
surfl(xx2);
shading interp
colormap(gray);
figure
plot(log10(error))

slope = exp((log(error(end))-log(error(2)))/(nitt-1));
slope3 = exp((log(error(end))-log(error(end-3)))/3);
slope5 = exp((log(error(end))-log(error(end-5)))/5);

fprintf('average convergence speed: %6f %6f %6f\n',slope,slope3,slope5);

end
