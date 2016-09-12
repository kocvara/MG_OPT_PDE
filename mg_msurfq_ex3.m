function [x]=mg_msurfq_ex3(nx,ny,levels)
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
% Solving the minimum surface problem with obstacle on a rectangular domain 
% by FAS multigrid with truncation
% nx,ny ... number of elements in x and y direction

global rho gp_count xstar fcount nnodesg nxxx nyyy levlev
global gamma_s_old

gamma_s_old=1;
rho=100; gp_count=0; fcount=0; nxxx=nx;nyyy=ny; levlev=levels;

[ivg,vxy,cf,elc]=rmeshl(nx,ny,levels);
[Avec]=Amat_surfq(levels,ivg,vxy);
[fl,intnodel]=RHSq_Ex3(levels,ivg,vxy,nx,ny);  
Pro=Pmatq(levels,ivg,vxy,cf,nx,ny,intnodel);

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;
a=zeros(nxf,1);b=zeros(nyf,1);c=zeros(nxf,1);d=zeros(nyf,1);
xxx=0:1/(nxf-1):1;
a=-sin(2*pi*xxx);b=-a;c=-a;d=a;

nnodesg{1} = (nx+1)*(ny+1);
for i=2:levels
    nnodesg{i} = size(vxy{i},2);
end

level=levels;

n = length(intnodel{levels});
x=zeros(n,1);
vh=zeros(n,1);

%% definition of the obstacle
aup=10; xu(1:n,1)=aup;

xxx=0:1/(nxf-1):1;yyy=0:1/(nyf-1):1;[X,Y]=meshgrid(xxx,yyy');
XLXL=-8.*(X-0.5).^2-8.*(Y-0.5).^2+.55;
xlxl=XLXL(:); xl=xlxl(intnodel{levels});


%% Exact solution
tic
fcount=0;xstar = zeros(size(x));
% opts    = struct( 'factr', 1e-1, 'pgtol', 1e-19, 'm', 9, 'maxIts',10000000,'maxTotalIts',100000000);
% opts.printEvery     = 100;
% [xstar, dumdum, info] = lbfgsb(@(x) msurfq_vh(x,Avec{levels},intnodel{levels},zeros(size(x)),levels), xl, xu, opts);

% exact solution from a file
name = strcat('xstar_ex4_lev',num2str(levels));
load(name)
funcevals

fcount=0;

error(1) = sqrt(sum((x-xstar).^2)/length(x));

toc
tic
%% Multigrid iterations
g=ones(size(x));

niter = 50;
eee=1;i=1;ndot=[];
while eee>4e-6 & i<niter
    vh = zeros(length(x),1);
    
    [x]=mgfas_msurfq_tr(levels,level,Avec,intnodel,cf,elc,Pro,vh,x,xl,xu,i);

    error(i+1) = sqrt(sum((x-xstar).^2)/length(x));
    
    eee = norm(error(i+1));
    i=i+1;

    fprintf('convergence speed, error: %3d %6f %6f\n',i,exp((log(error(i))-log(error(i-1)))),error(i));
end
nitt=i-1;

%%
toc
gp_count
fcount

nxf=nx*2^(levels-1)+1;nyf=ny*2^(levels-1)+1;
xx2=zeros(nxf,nyf);
xx2(:,1)=a; xx2(end,:)=b; xx2(:,end)=c; xx2(1,:)=d;
xx2(2:end-1,2:end-1) = reshape(xstar,nyf-2,nxf-2);
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
