function [x,fcount]=mg_pde_eq_ex4(nx,ny,levels)
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
% Solving Example 4 with an equality constraint
% Discretization by square bilinear finite elements
% nx,ny ... number of elements in x and y direction

global rho gp_count xstar fcount nnodesg nxxx nyyy levlev
%global a b c d
global gamma_s_old

gamma_s_old=1;
rho=100; gp_count=0; fcount=0; nxxx=nx;nyyy=ny; levlev=levels;

tic
[ivg,vxy,cf]=rmeshl(nx,ny,levels);
%[Aglob,~,nelem]=Amatq(levels,ivg,vxy);
[fl,intnodel]=RHSq_Ex4(levels,ivg,vxy,nx,ny);
%fl{levels} = ones(length(intnodel{levels}),1);
Pro=Pmatq(levels,ivg,vxy,cf,nx,ny,intnodel);

load Amat_file

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;

for i=2:levels
    %pp = Pro{i};
    %Pro{i} = pp(intnodel{i-1},intnodel{i});
    nnodesg{i} = size(vxy{i},2);
end

for i=1:levels
    nxf=nx*2^(i-1);nyf=ny*2^(i-1);
    nelem{i} = nxf*nyf;
end

level=levels;

n = length(intnodel{levels});
x=zeros(n,1);
vh=zeros(n,1);
%%
for l=1:levels
    intel=intnodel{l};
    dintel=length(intel);
    Area{l}=1/nelem{l};
    A_r=Area{l};
    Area_v{l}=A_r*ones(dintel,1);
end

% RIGHT HAND SIDE FOR VOLUME (EQUALITY) CONSTRAINT
vol=1;                                 

maxnelem = length(intnodel{levels});

% Initial point must be feasible in the equality constraint
x=vol/maxnelem./Area_v{levels};

%% Definition of the obstacle
nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;
alow=-10; %alow=-0.05;
aup=10;  %aup=1.8;% aup=0.05;
% simple bounds
xl(1:n,1)=alow; xu(1:n,1)=aup;

% parabolic lower bound
xxx=0:1/(nxf-1):1;yyy=0:1/(nyf-1):1;[X,Y]=meshgrid(xxx,yyy');
XLXL=-4.*(X-0.5).^2-4.*(Y-0.5).^2+.1;
XLXL=-32.*(X-0.5).^2-32.*(Y-0.5).^2+2.5;
xlxl=XLXL(:);xl=xlxl(intnodel{levels});


%% Exact solution
Aco = Aglob{levels};Aco=Aco(intnodel{levels},intnodel{levels});
% tic,[xstar,kk]=gp_pde1_con_eq(x,Aco,fl{levels},intnodel{levels},zeros(size(x)),10000,1e-13,xl,xu,levels,Area_v{levels});toc
%kk
%fcount
name = strcat('xstar_ex5_lev',num2str(levels));
load(name)
%funcevals
error(1) = sqrt(sum((x-xstar).^2)/length(x));
%error
tic
%% Multigrid iterations
g=ones(size(x));

niter = 50;
eee=1;i=1;ndot=[];
while eee>2e-6 & i<niter
    vh = zeros(length(x),1); %vh=fl{levels};
    [x]=mgfas_pde_eq_hm_Ex4(levels,level,Aglob,fl,intnodel,cf,Pro,vh,x,xl,xu,i,Area_v);
    
    error(i+1) = sqrt(sum((x-xstar).^2)/length(x));
    
    eee = norm(error(i+1));
    i=i+1;
    %
    fprintf('convergence speed, error: %3d %6f %6f\n',i,exp((log(error(i))-log(error(i-1)))),eee);
end
nitt=i-1;

%%
toc

%gp_count
fcount

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;
xx2=zeros(nxf,nyf);
xx2(2:end-1,2:end-1) = reshape(x,nyf-2,nxf-2);
%surf(xx2);
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
