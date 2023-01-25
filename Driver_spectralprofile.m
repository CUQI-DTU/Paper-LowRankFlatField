% Low-rank flat-field correction for artifact reduction in 
% spectral computed tomography
%
% Driver for producing Figure 7 and Table 4
%
% Katrine O. Bangsgaard, DTU, 05/05-2022

clear all, close all, clc

% Predefined colors
mucol = [31, 119, 180; ...
255, 127, 14; ...
44, 160, 44; ...
214, 39, 40; ...
148, 103, 189]/255;

% Save figures path
graphicpath = './figures/';

% Add paths
addpath('./functions')
addpath('./data')

% Load and store data
load('neutrondata.mat')
load('angstrom_bins.mat')
sino = double(permute(sino,[2 3 1]));
flat = double(permute(flat,[2 3 1]));

% System specifications
[p,r,k]    = size(sino);
pixelsize  = 0.055;
dsz        = pixelsize*r;
theta      = linspace(0,p-1,p)'*1.5;
s          = 8;
n          = r;

% Generate system matrix
A_n        = paralleltomo(n,theta,n,n);
A_n        = (dsz/n)*A_n;  

% Generate mask
[Xt,Yt]     = meshgrid(linspace(-1,1,n),linspace(-1,1,n));
mask_radius = 1;
mask        = Xt.^2 + Yt.^2 > mask_radius;
X(mask)     = 0;
clear('Xt','Yt');

% Load attenuation coefficients
data = importdata('attenuation.txt');
data2= importdata('Zn_new.txt');

% Store energies
ENG    = data.data(:,1);
ENG2   = data2(:,1);

% Store materials
Fe     = data.data(:,3);
Ni     = data.data(:,4);
Cu     = data.data(:,5);
Al     = data.data(:,9);
Zn     = data2(:,2);
Theory = [Fe,Ni,Cu,Al,Al];
%% Identify and fix sinograms with "0 measurements" 
id = zeros(k,1);
for i=1:k
    test = sino(:,:,i);
    if sum(any(test(:) == 0)) > 0    
    id(i) = 1;
    end
end

id = find(id > 0);
testsino= reshape(sino(:,:,id),p,r,length(id));
masksino= zeros(p,r,length(id));
for i=1:length(id)
    [idr, idc] = find(testsino(:,:,i) == 0);
    
    if ~isempty(idr)
    
    for j =1:length(idr)
        if idr(j) == 1
        testsino(idr(j),idc(j),i) = 1/3*(testsino(idr(j)+1,idc(j),i) ...
                +testsino(idr(j),idc(j)+1,i)+testsino(idr(j),idc(j)-1,i));    
        elseif idr(j) == p
        testsino(idr(j),idc(j),i) = 1/3*(testsino(idr(j)-1,idc(j),i) ...
                +testsino(idr(j),idc(j)+1,i)+testsino(idr(j),idc(j)-1,i));
        elseif idc(j) == 1
        testsino(idr(j),idc(j),i) = 1/3*(testsino(idr(j)+1,idc(j),i)+testsino(idr(j)-1,idc(j),i) ...
                +testsino(idr(j),idc(j)+1,i));
        elseif idc(j) == r
        testsino(idr(j),idc(j),i) = 1/3*(testsino(idr(j)+1,idc(j),i)+testsino(idr(j)-1,idc(j),i) ...
                +testsino(idr(j),idc(j)-1,i));   
        else
        testsino(idr(j),idc(j),i) = 1/4*(testsino(idr(j)+1,idc(j),i)+testsino(idr(j)-1,idc(j),i) ...
                +testsino(idr(j),idc(j)+1,i)+testsino(idr(j),idc(j)-1,i));
        end
    end
    end
end

sino(:,:,id) = testsino;
%% Compute FBP and TV reconstructions for all energies (time consuming use 
% [U,S,V]   = svd(reshape(flat(:,:,:),s*r,k),'econ');
% SS        = diag(S);
% 
% flatappr  = reshape(U(:,1)* S(1,1) * V(:,1)',s,r,k);
% 
% pixel=zeros(n*n,k,4);
% for i = 1:k
% F = squeeze(mean(flat(:,:,i),[1 3]));
% Y = squeeze(mean(sino(:,:,i),3)); 
% b = -log(Y./F)';
% x = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n);
% x(mask) = 0;    
% 
% bp = xRemoveStripesVertical(b',3,'db5',0.9)';
% xp = reshape(fbp(A_n,bp(:),theta)*(n/dsz)^2,n,n);
% xp(mask) = 0;  
% 
% bp2 = tsino_nonlocal_medfilt(b, 31);
% xp2 = reshape(fbp(A_n,bp2(:),theta)*(n/dsz)^2,n,n);
% xp2(mask) = 0;  
% 
% F  = squeeze(mean(flatappr(:,:,i),[1 3]))';
% blr = -log(Y./F')';
% xlr = reshape(fbp(A_n,blr(:),theta)*(n/dsz)^2,n,n);
% xlr(mask) = 0;
% 
% 
% pixel(:,i,1) = x(:);
% pixel(:,i,2) = xp(:);
% pixel(:,i,3) = xp2(:);
% pixel(:,i,4) = xlr(:);
% display(i);
% end
% % TV for all energies
% pixelTV=zeros(n*n,k,4);
% for i = 1:k
%     
% F = squeeze(mean(flat(:,:,i),[1 3]))';
% Flr = squeeze(mean(flatappr(:,:,i),[1 3]))';
% Y = squeeze(mean(sino(:,:,i),3))'; 
% clear options;
% 
% maxiters = 1000;    % Number of iterations
% options = struct(...
%     'u0',zeros(n*n,1),...
%     'rho',1.8,...
%     'tau',1e-3,...
%     'tolf',1e-8,... 
%     'lambda',0.005,...%5e-2,...
%     'maxiters',maxiters,...
%     'uhold',[],...
%     'mask',true,...
%     'verbose',0);
% 
% xlr = gd_wls(A_n,Flr,Y,options); % SVD 
% 
% x= gd_wls(A_n,F,Y,options); % FBP 
% 
% options.preproc = 1;
% xp  = gd_wls(A_n,F,Y,options);   % P-TV
% options.preproc = 2;
% xp2  = gd_wls(A_n,F,Y,options);   % P-TV
% 
% 
% 
% pixelTV(:,i,1) = x(:);
% pixelTV(:,i,2) = xp(:);
% pixelTV(:,i,3) = xp2(:);
% pixelTV(:,i,4) = xlr(:);
% display(i);
% end
% 
% FBPrec = pixel(:,:,1);save('FBPrec.mat','FBPrec')
% PFBPrec = pixel(:,:,2);save('PFBPrec.mat','PFBPrec')
% P2FBPrec = pixel(:,:,3);save('P2FBPrec.mat','P2FBPrec')
% LRFBPrec = pixel(:,:,4);save('LRFBPrec.mat','LRFBPrec')
% TVrec = pixelTV(:,:,1);save('TVrec.mat','TVrec')
% PTVrec = pixelTV(:,:,2);save('PTVrec.mat','PTVrec')
% P2TVrec = pixelTV(:,:,3);save('P2TVrec.mat','P2TVrec')
% LRTVrec = pixelTV(:,:,4);save('LRTVrec.mat','LRTVrec')

%% Load mat files for spectral profiles
addpath('./spectralplot/')
load('FBPrec.mat')
load('PFBPrec.mat')
load('P2FBPrec.mat')
load('LRFBPrec.mat')
load('TVrec.mat')
load('PTVrec.mat')
load('P2TVrec.mat')
load('LRTVrec.mat')

TVtest = reshape(TVrec,n,n,k)/pixelsize;
PTVtest = reshape(PTVrec,n,n,k)/pixelsize;
P2TVtest = reshape(P2TVrec,n,n,k)/pixelsize;
LRTVtest = reshape(LRTVrec,n,n,k)/pixelsize;

FBPtest = reshape(FBPrec,n,n,k)/pixelsize;
PFBPtest = reshape(PFBPrec,n,n,k)/pixelsize;
P2FBPtest = reshape(P2FBPrec,n,n,k)/pixelsize;
LRFBPtest = reshape(LRFBPrec,n,n,k)/pixelsize;

%% Figure (not in paper): Position of the chosen pixels for spectral profiles
load('angstrom_bins.mat')
roi_row = [362, 242, 120, 132, 252];
roi_col = [182, 382, 150, 312, 100];
                    
plottest = reshape(TVrec,n,n,k);
figure
imagesc(plottest(:,:,200))
axis square
axis off
hold on;
for i=1:length(roi_row)
plot(roi_col(i),roi_row(i), '.','color',mucol(i,:), 'MarkerSize', 30, 'LineWidth', 2);
end
%saveas(gcf,fullfile(graphicpath,'Experiment3POINT'),'epsc')

%% Figure 7: Spectral plot TV

coeff  = [0.65, 0.9, 0.55, 1, 0.9];
addrow = 1;

% plot settings
ystart    = [0.82, 0.64, 0.46, 0.28 0.1];
xstart    = [0.11 0.32 0.58 0.79];
dim       = 0.4;
dim2      = 0.15;
cmin      = [0.4, 1.3,0.0,0.0,0.0];
cmax      = [2.1,2.8,1.7,0.7,0.4];
cmin2     = [0.4, 1.4,0.7,0.0,0.0];
cmax2     = [2.1,2.8,1.5,0.5,0.3];
xmin      = [2.2 3.3 3.4 4.0 3.9];
xmax      = [3.2 4.3 4.4 5.0 4.9];
ylabstr   = {'Fe','Ni','Cu','Zn','Al'};
linw      = 0.5;

figure('Position',[100 100 1000 1400]) 


for i =1:5

    if i == 4
        ax1 = axes('Position',[xstart(1) ystart(i) dim dim2]);
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(1,:),'linewidth',linw)
        hold on
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(PTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(2,:),'linewidth',linw)
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(P2TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(3,:),'linewidth',linw)
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(LRTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(4,:),'linewidth',linw)

        plot(ax1,ENG2,Zn,'color','k','linewidth',linw)
        hold off    
        
        ax2 = axes('Position',[xstart(3) ystart(i) dim dim2]);
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(1,:),'linewidth',linw)
        hold on
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(PTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(2,:),'linewidth',linw)
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(P2TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(3,:),'linewidth',linw)
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(LRTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(4,:),'linewidth',linw)

        plot(ax2,ENG2,Zn,'color','k','linewidth',linw)
        hold off  
    else     
        ax1 = axes('Position',[xstart(1) ystart(i) dim dim2]);
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(1,:),'linewidth',linw)
        hold on
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(PTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(2,:),'linewidth',linw)
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(P2TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(3,:),'linewidth',linw)
        plot(ax1,angstrom,1/coeff(i)*squeeze(mean(LRTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(4,:),'linewidth',linw)

        plot(ax1,ENG,Theory(:,i),'color','k','linewidth',linw)
        hold off
        
        ax2 = axes('Position',[xstart(3) ystart(i) dim dim2]);
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(1,:),'linewidth',linw)
        hold on
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(PTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(2,:),'linewidth',linw)
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(P2TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(3,:),'linewidth',linw)
        plot(ax2,angstrom,1/coeff(i)*squeeze(mean(LRTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2])),'color',mucol(4,:),'linewidth',linw)

        plot(ax2,ENG,Theory(:,i),'color','k','linewidth',linw)
        hold off  
    
    end
 
    
xlim(ax1,[1 5])
xlim(ax2,[xmin(i) xmax(i)])

ylim(ax1,[cmin(i) cmax(i)])
ylim(ax2,[cmin2(i) cmax2(i)])

ylabel(ax1,ylabstr{i})
if  i == 3
    ylabel(ax1,{'$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]', 'Cu'})
    
elseif i == 5
    xticklabels(ax1,[1,2,3,4,5])
    xlabel(ax1,'\AA')
    xlabel(ax2,'\AA')
    g = legend(ax1,'TV','WF-TV','NLM-TV','LR-TV','location','northwest','fontsize',20);
    
end    
set(ax1,'box','off','XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
set(ax2,'box','off','XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
end
set(gcf, 'PaperPositionMode','manual')
saveas(gcf,fullfile(graphicpath,['Experiment3TV', '2']),'epsc')

%% Table 4: Mean square error for Figure 7
MSE = zeros(5,4);

for i=1:5
    Mat = zeros(k,4);
    Mat(:,1) = 1/coeff(i)*squeeze(mean(TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2]));
    Mat(:,2) = 1/coeff(i)*squeeze(mean(PTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2]));
    Mat(:,3) = 1/coeff(i)*squeeze(mean(P2TVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2]));
    Mat(:,4) = 1/coeff(i)*squeeze(mean(LRTVtest(roi_row(i):roi_row(i)+addrow,roi_col(i):roi_col(i)+addrow,:),[1 2]));
    if i == 4
    xint       = ENG2;
    yint       = Zn;
    xq         = angstrom;
    vq         = interp1(xint,yint,xq);    
    MSE(i,1) = norm(vq(2:k-1)-Mat(2:k-1,1))/norm(vq(2:k-1));
    MSE(i,2) = norm(vq(2:k-1)-Mat(2:k-1,2))/norm(vq(2:k-1));
    MSE(i,3) = norm(vq(2:k-1)-Mat(2:k-1,3))/norm(vq(2:k-1));
    MSE(i,4) = norm(vq(2:k-1)-Mat(2:k-1,4))/norm(vq(2:k-1));
    
    else
    xint       = ENG;
    yint       = Theory(:,i);
    xq         = angstrom;
    vq         = interp1(xint,yint,xq);    
    MSE(i,1) = norm(vq-Mat(:,1))/norm(vq);
    MSE(i,2) = norm(vq-Mat(:,2))/norm(vq);
    MSE(i,3) = norm(vq-Mat(:,3))/norm(vq);
    MSE(i,4) = norm(vq-Mat(:,4))/norm(vq);
    end
    
end

disp(MSE)

