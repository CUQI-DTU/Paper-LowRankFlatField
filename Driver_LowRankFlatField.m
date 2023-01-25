% Low-rank flat-field correction for artifact reduction in 
% spectral computed tomography
%
% Driver for producing all Figures except Figure 7.
%
% Katrine O. Bangsgaard, DTU, 05/05-2022

clear all, close all, clc

% Predefined colors
red     = [0.60,  0.00,  0.00];
green   = [0.4392  0.6392 0.4392];
orange  = [0.8    0.33   0.0];
blue    = [0.3725 0.6588 0.8431];

% Save figures path
graphicpath = './figures/';

% Add paths
addpath('./functions')
addpath('./data')

% Choose colormap 
cmap = 0; % standard
%cmap = 1; % inferno

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

%% Figure 1: Naive FBP rec for two energies
figure('Position',[100 100 1000 1000]) 

% Energy channel 145
i  = 145;
F  = squeeze(mean(flat(:,:,i),[1 3]));
Y  = squeeze(mean(sino(:,:,i),3)); 
b  = -log(Y./F)';
x  = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;    

ax = axes('Position',[0.05 0.55 0.3 0.3]);
imagesc(ax,x)
title(ax,['(a) ',num2str(round(angstrom(i),1)), ' \AA'])
caxis(ax,[-0.3 2.5]),
set(ax, 'XTick', [],'YTick', []);


% Energy channel 200
i =  200;
F = squeeze(mean(flat(:,:,i),[1 3]));
Y = squeeze(mean(sino(:,:,i),3)); 
b = -log(Y./F)';
xpm = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpm(mask) = 0;    

ax = axes('Position',[0.37 0.55 0.3 0.3]);
imagesc(ax,xpm)
caxis(ax,[-0.3 2.5]), 
set(ax, 'XTick', [],'YTick', []);
title(ax,['(b) ',num2str(round(angstrom(i),1)), ' \AA']) 
h = colorbar('location','eastoutside','Position',[0.69 0.55 0.02 0.3],'fontsize',20);
h.Label.Interpreter = 'latex';
h.Label.String = '$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]';

if cmap == 1
colormap inferno
end
saveas(gcf,fullfile(graphicpath,'FBPrec'),'epsc')
%% Figure 2: Flat-field analysis
% Singular value decomposition
[U,S,V]   = svd(reshape(flat(:,:,:),s*r,k),'econ');
SS = diag(S);

% colorbar specifications
cmin = 0.01;
cmax = 0.06;

figure('Position',[100 100 1000 1000]) 
ax1 = axes('Position',[0.07 0.72 0.25 0.2]);
ax2 = axes('Position',[0.4 0.72 0.25 0.2]);
imagesc(ax1,reshape(flat,s*r,k)),caxis(ax1,[cmin cmax])
title(ax1,'(a) Flat-field measurements')
xlabel(ax1,'\AA')
ylabel(ax1,'detector $\times$ measurement')
colorbar(ax1,'location','eastoutside','Position',...
      [0.33 0.72 0.02 0.2],'fontsize',20);
set(ax1,'XTick',1:84:339,'XTickLabels', [1 2 3 4 5])

if cmap == 1    
    loglog(ax2,SS,'-','linewidth',3,'marker','.','color',red,'markersize',30)
else   
    loglog(ax2,SS,'-','linewidth',3,'marker','.','color',blue,'markersize',30)
end
set(ax2,'YTick',[1e-2, 1e-1 1e-0 1e+1 1e+2]);
title(ax2,'(b) Singular values')
xlim(ax2,[0 k])
ylim([min(SS) max(SS)])
set(ax2,'XTick',[1 1e+1 1e+2],'box','off');
set(ax2,'YAxisLocation', 'right')
legend(ax2,'$\sigma_k$')
xlabel(ax2,'$k$')
grid on
if cmap == 1
    colormap inferno
end
saveas(gcf,fullfile(graphicpath,'SVDintro'),'epsc')

%% Figure 3: SVD analysis
% rank-one approximation
flatappr  = reshape(U(:,1)* S(1,1) * V(:,1)',s,r,k);
%rank five-approximation
flatappr5 = reshape(U(:,1:5)* S(1:5,1:5) * V(:,1:5)',s,r,k);

%figure settings
cmin = 0.01;
cmax = 0.06;

figure('Position',[100 100 1000 1000]) 
ax1 = axes('Position',[0.07 0.72 0.25 0.2]);
ax2 = axes('Position',[0.33 0.72 0.25 0.2]);
ax3 = axes('Position',[0.59 0.72 0.25 0.2]); 
ax4 = axes('Position',[0.07 0.43 0.25 0.2]);
ax5 = axes('Position',[0.33 0.43 0.25 0.2]);
ax6 = axes('Position',[0.59 0.43 0.25 0.2]); 

imagesc(ax1,reshape(flat,s*r,k)),caxis(ax1,[cmin cmax])
title(ax1,'(a) Spectral flat-fields')
set(ax1,'XTick',[10    85   169   253   328],'XTickLabels', [1 2 3 4 5])
xlabel(ax1,'\AA')
ylabel(ax1,'detector $\times$ measurement')

imagesc(ax2,reshape(flatappr,s*r,k)),caxis(ax2,[cmin cmax])
title(ax2,'(b) Rank-one approx.')
set(ax2, 'YTick', [],'XTick',[10    85   169   253   328],'XTickLabels', [1 2 3 4 5])
xlabel(ax2,'\AA')

imagesc(ax3,reshape(flatappr5,s*r,k)),caxis(ax2,[cmin cmax])
title(ax3,'(e) Rank-five approx.')
caxis(ax3,[cmin cmax])
set(ax3, 'YTick', [],'XTick',[10    85   169   253   328],'XTickLabels', [1 2 3 4 5])
xlabel(ax3,'\AA')

colorbar(ax3,'location','eastoutside','position',[0.850 0.720 0.0160 0.2000],'fontsize',20)

imagesc(ax5,reshape(flatappr-flat,s*r,k))
title(ax5,'(e) Rank-one difference')
set(ax5, 'YTick', [],'XTick',[10    85   169   253   328],'XTickLabels', [1 2 3 4 5])
xlabel(ax5,'\AA')
caxis(ax5,[0 0.01])


imagesc(ax6,reshape(flatappr5-flat,s*r,k)),caxis(ax6,[0 0.01])
title(ax6,'(f) Rank-five difference')
colorbar(ax6,'location','eastoutside','position',[0.850 0.4300 0.0160 0.2000],'fontsize',20)
set(ax6, 'YTick', [],'XTick',[10    85   169   253   328],'XTickLabels', [1 2 3 4 5])
xlabel(ax6,'\AA')

if cmap == 1
    loglog(ax4,SS/SS(1),'-','linewidth',3,'marker','.','color',red,'markersize',30)
else
    loglog(ax4,SS/SS(1),'-','linewidth',3,'marker','.','color',blue,'markersize',30)

end
set(ax4,'YTick',[1e-3 1e-2, 1e-1 1e-0 1e+1 1e+2]);
xlim(ax4,[0 k])
xlabel(ax4,'$k$')
legend(ax4,'$\sigma_k / \sigma_1$')
set(ax4,'XTick',[1 1e+1 1e+2],'box','off');
set(ax4,'XGrid','on','YGrid','on')
title(ax4,'(d) Relative singular value')
grid on
if cmap == 1
colormap inferno
end
saveas(gcf,fullfile(graphicpath,'SVDrank15'),'epsc')

%% Relative difference measure eq. (8) in paper 
RD=zeros(k,1);
for i = 1:k
% FBP reconstruction    
F  = squeeze(mean(flat(:,:,i),[1 3]));
Y  = squeeze(mean(sino(:,:,i),3)); 
b  = -log(Y./F)';
x  = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;    

% Low rank reconstruction
F   = squeeze(mean(flatappr(:,:,i),[1 3]))';
blr = -log(Y./F')';
xlr = reshape(fbp(A_n,blr(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xlr(mask) = 0;

% Relative difference
RD(i) = norm(x(:)-xlr(:),2)/norm(xlr(:),2);
display(i);
end

% compute statistics for relative difference
[~,minid] = min(RD);
[~,maxid] = max(RD);
y         = quantile(RD,[0.25  0.50 0.75]);
[~,q1]    = min(abs(RD-y(1)));
[~,q2]    = min(abs(RD-y(2)));
[~,q3]    = min(abs(RD-y(3)));

%% Figure 5: Different energies FBP

%CNR settings 
yb = n/2-30;
ys = 160;
xb = n/2-30;
xs = n/2+90;
yb2  = 280;
xb2  = n/2+90;

cmin =-0.3;
cmax = 2.5;
ystart    = [0.72, 0.51, 0.30, 0.09];
xstart    = [0.05 0.26 0.47 0.68];
dim = 0.2;
figure('Position',[100 100 900 900]) ;

% Minimum R3
F   = squeeze(mean(flat(:,:,minid),[1 3]));
Flr = squeeze(mean(flatappr(:,:,minid),[1 3]));
Y   = squeeze(mean(sino(:,:,minid),3)); 


b = -log(Y./F)';
x = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;    


bpv = tsino_nonlocal_medfilt(b, 31);
xpv = reshape(fbp(A_n,bpv(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpv(mask) = 0;  

bpm = xRemoveStripesVertical(b',3,'db5',0.9)';
xpm = reshape(fbp(A_n,bpm(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpm(mask) = 0;

blr   = -log(Y./Flr)';
xlr   = reshape(fbp(A_n,blr(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xlr(mask) = 0;    


sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(1) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]),  title('FBP')
set(ax1, 'XTick', [], 'YTick', []);
ylabel(ax1,[num2str(round(angstrom(minid),1)), ' \AA'])

ax2 = axes('Position',[xstart(2) ystart(1) dim dim]);
imagesc(ax2,xpm), 
caxis([cmin cmax]),  
set(ax2, 'XTick', [], 'YTick', []);
title(ax2,'WF-FBP')

ax3 = axes('Position',[xstart(3) ystart(1) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);
title(ax3,'NLM-FBP')
ax4 = axes('Position',[xstart(4) ystart(1) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);
title(ax4,'LR-FBP')

% plot CNR ROI
  hold on 
rectangle(ax4,'Position',[ys,xs,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','-')
rectangle(ax4,'Position',[yb2,xb2,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','--')

hold off 

% Median RD
F   = squeeze(mean(flat(:,:,q2),[1 3]));
Flr = squeeze(mean(flatappr(:,:,q2),[1 3]));
Y   = squeeze(mean(sino(:,:,q2),3)); 

b = -log(Y./F)';
x = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;    


bpv = tsino_nonlocal_medfilt(b, 31);
xpv = reshape(fbp(A_n,bpv(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpv(mask) = 0;  

bpm = xRemoveStripesVertical(b',3,'db5',0.9)';
xpm = reshape(fbp(A_n,bpm(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpm(mask) = 0;

blr   = -log(Y./Flr)';
xlr   = reshape(fbp(A_n,blr(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xlr(mask) = 0;    

sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(2) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]), 

ylabel(ax1,[num2str(round(angstrom(q2),1)), ' \AA'])
set(ax1, 'XTick', [], 'YTick', []);

ax2 = axes('Position',[xstart(2) ystart(2) dim dim]);
imagesc(ax2,xpm), 
caxis([cmin cmax]), 
set(ax2, 'XTick', [], 'YTick', []);

ax3 = axes('Position',[xstart(3) ystart(2) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);

ax4 = axes('Position',[xstart(4) ystart(2) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);


% Max RD
F   = squeeze(mean(flat(:,:,maxid),[1 3]));
Flr = squeeze(mean(flatappr(:,:,maxid),[1 3]));
Y   = squeeze(mean(sino(:,:,maxid),3)); 


b = -log(Y./F)';
x = reshape(fbp(A_n,b(:),theta)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;    


bpv = tsino_nonlocal_medfilt(b, 31);
xpv = reshape(fbp(A_n,bpv(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpv(mask) = 0;  

bpm = xRemoveStripesVertical(b',3,'db5',0.9)';
xpm = reshape(fbp(A_n,bpm(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xpm(mask) = 0;

blr   = -log(Y./Flr)';
xlr   = reshape(fbp(A_n,blr(:),theta)*(n/dsz)^2,n,n)/pixelsize;
xlr(mask) = 0;    

sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(3) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]), 
ylabel(ax1,[num2str(round(angstrom(maxid),1)), ' \AA'])
set(ax1, 'XTick', [], 'YTick', []);

ax2 = axes('Position',[xstart(2) ystart(3) dim dim]);
imagesc(ax2,xpm),  
caxis([cmin cmax]), 
set(ax2, 'XTick', [], 'YTick', []);

ax3 = axes('Position',[xstart(3) ystart(3) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);

ax4 = axes('Position',[xstart(4) ystart(3) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);

if cmap == 1
colormap inferno
end

h = colorbar('location','eastoutside','Position',[0.89 ystart(3) 0.02 ystart(1)-0.1],'fontsize',20);
h.Label.Interpreter = 'latex';
h.Label.String = '$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]';

saveas(gcf,fullfile(graphicpath,'Experiment1FBP'),'epsc')  


    
%% Figure 6: Different energies TV
figure('Position',[100 100 1000 1000]) ;

mucol = [31, 119, 180; ...
255, 127, 14; ...
44, 160, 44; ...
214, 39, 40; ...
148, 103, 189]/255;

tau = 1e-3;
lamb = 0.005;
% Minimum R3
F   = squeeze(mean(flat(:,:,minid),[1 3]))';
Flr = squeeze(mean(flatappr(:,:,minid),[1 3]))';
Y   = squeeze(mean(sino(:,:,minid),3))'; 

clear options;

maxiters = 1000;
options = struct(...
    'u0',zeros(n*n,1),...
    'rho',1.8,...
    'tau',tau,...
    'lambda',lamb,...
    'tolf',1e-4,... 
    'maxiters',maxiters,...
    'uhold',[],...
    'mask',true,...
    'verbose',0);

xlr = gd_wls(A_n,Flr,Y,options)/pixelsize; % LR-TV 
x   = gd_wls(A_n,F,Y,options)/pixelsize;   % TV 
options.preproc = 1;
xpm  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV
options.preproc = 2;
xpv  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV

x   = reshape(x,n,n);
xpm  = reshape(xpm,n,n);
xpv  = reshape(xpv,n,n);
xlr = reshape(xlr,n,n);  


sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(1) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]),  title('TV')
set(ax1, 'XTick', [], 'YTick', []);

ylabel(ax1,[num2str(round(angstrom(minid),1)), ' \AA'])
ax2 = axes('Position',[xstart(2) ystart(1) dim dim]);
imagesc(ax2,xpm), 
caxis([cmin cmax]),  
set(ax2, 'XTick', [], 'YTick', []);
title(ax2,'WF-TV')

ax3 = axes('Position',[xstart(3) ystart(1) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);
title(ax3,'NLM-TV')
ax4 = axes('Position',[xstart(4) ystart(1) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);
title(ax4,'LR-TV')

% plot CNR ROI
  hold on 
rectangle(ax4,'Position',[ys,xs,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','-')
rectangle(ax4,'Position',[yb2,xb2,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','--')
roi_row = [362, 242, 120, 132, 252];
roi_col = [182, 382, 150, 312, 100];
                    

for i=1:length(roi_row)
ppp=plot(roi_col(i),roi_row(i), 'o','color','white', 'MarkerSize', 8, 'LineWidth', 1);
ppp.MarkerFaceColor = 'white';%,mucol(i,:);
end
hold off 
% Median RD
F   = squeeze(mean(flat(:,:,q2),[1 3]))';
Flr = squeeze(mean(flatappr(:,:,q2),[1 3]))';
Y   = squeeze(mean(sino(:,:,q2),3))'; 
clear options;

options = struct(...
    'u0',zeros(n*n,1),...
    'rho',1.8,...
    'tau',tau,...
    'lambda',lamb,...
    'tolf',1e-4,... 
    'maxiters',maxiters,...
    'uhold',[],...
    'mask',true,...
    'verbose',0);
xlr = gd_wls(A_n,Flr,Y,options)/pixelsize; % LR-TV 
x   = gd_wls(A_n,F,Y,options)/pixelsize;   % TV 
options.preproc = 1;
xpm  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV
options.preproc = 2;
xpv  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV

x    = reshape(x,n,n);
xpm  = reshape(xpm,n,n);
xpv  = reshape(xpv,n,n);
xlr  = reshape(xlr,n,n);  


sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(2) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]), 


ylabel(ax1,[num2str(round(angstrom(q2),1)), ' \AA'])
set(ax1, 'XTick', [], 'YTick', []);

ax2 = axes('Position',[xstart(2) ystart(2) dim dim]);
imagesc(ax2,xpm), 
caxis([cmin cmax]), 
set(ax2, 'XTick', [], 'YTick', []);

ax3 = axes('Position',[xstart(3) ystart(2) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);

ax4 = axes('Position',[xstart(4) ystart(2) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);


% Max RD
F   = squeeze(mean(flat(:,:,maxid),[1 3]))';
Flr = squeeze(mean(flatappr(:,:,maxid),[1 3]))';
Y   = squeeze(mean(sino(:,:,maxid),3))'; 

clear options;

options = struct(...
    'u0',zeros(n*n,1),...
    'rho',1.8,...
    'tau',tau,...
    'lambda',lamb,...
    'tolf',1e-4,... 
    'maxiters',maxiters,...
    'uhold',[],...
    'mask',true,...
    'verbose',0);
xlr = gd_wls(A_n,Flr,Y,options)/pixelsize; % LR-TV 
x   = gd_wls(A_n,F,Y,options)/pixelsize;   % TV 
options.preproc = 1;
xpm  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV
options.preproc = 2;
xpv  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV

x   = reshape(x,n,n);
xpm  = reshape(xpm,n,n);
xpv  = reshape(xpv,n,n);
xlr = reshape(xlr,n,n);  


sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])

ax1 = axes('Position',[xstart(1) ystart(3) dim dim]);
imagesc(ax1,x), 
caxis([cmin cmax]),  

ylabel(ax1,[num2str(round(angstrom(maxid),1)), ' \AA'])
set(ax1, 'XTick', [], 'YTick', []);

ax2 = axes('Position',[xstart(2) ystart(3) dim dim]);
imagesc(ax2,xpm),  
caxis([cmin cmax]), 
set(ax2, 'XTick', [], 'YTick', []);

ax3 = axes('Position',[xstart(3) ystart(3) dim dim]);
imagesc(ax3,xpv), 
caxis([cmin cmax]), 
set(ax3, 'XTick', [], 'YTick', []);

ax4 = axes('Position',[xstart(4) ystart(3) dim dim]);
imagesc(ax4,xlr), 
caxis([cmin cmax]), 
set(ax4, 'XTick', [], 'YTick', []);
if cmap == 1
colormap inferno
end
h = colorbar('location','eastoutside','Position',[0.89 ystart(3) 0.02 ystart(1)-0.1],'fontsize',20);
h.Label.Interpreter = 'latex';
h.Label.String = '$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]';
saveas(gcf,fullfile(graphicpath,'Experiment1TV'),'epsc')  

%% Figure 8: Number of flat-fields FBP
% plot settings 
flatcount = [1 2 4 8];
ystart    = [0.72, 0.51, 0.30, 0.09];
xstart    = [0.05 0.26 0.47 0.68];
cmin      =-0.3;
cmax      = 2.5;
dim       = 0.2;
filter = 'Hann';
figure('Position',[100 100 1000 1000]) 


for i=1:length(flatcount)
% Compute SVD approximation     
[U,S,V]   = svd(reshape(flat(1:flatcount(i),:,:),flatcount(i)*r,k),'econ');
flatappr2 = reshape(U(:,1)* S(1,1) * V(:,1)',flatcount(i),r,k);

ax1 = axes('Position',[xstart(1) ystart(5-i) dim dim]);
ax2 = axes('Position',[xstart(2) ystart(5-i) dim dim]);
ax3 = axes('Position',[xstart(3) ystart(5-i) dim dim]); 
ax4 = axes('Position',[xstart(4) ystart(5-i) dim dim]); 
energy   = 300; % 30, 90
F        = squeeze(mean(flat(1:flatcount(i),:,energy),[1 3]));
Y        = squeeze(mean(sino(:,:,energy),3)); 
% FBP
b        = -log(Y./F)';
x       = reshape(fbp(A_n,b(:),theta,filter)*(n/dsz)^2,n,n)/pixelsize;
x(mask) = 0;

% Pm-FBP
bpm      = xRemoveStripesVertical(b',3,'db5',0.9)';
xpm       = reshape(fbp(A_n,bpm(:),theta,filter)*(n/dsz)^2,n,n)/pixelsize;
xpm(mask) = 0;  

% Pv-FBP
bpv       = tsino_nonlocal_medfilt(b, 31);
xpv       = reshape(fbp(A_n,bpv(:),theta,filter)*(n/dsz)^2,n,n)/pixelsize;
xpv(mask) = 0; 

%LR-FBP
F        = squeeze(mean(flatappr2(1:flatcount(i),:,energy),[1 3]))';
blr       = -log(Y./F')';
xlr       = reshape(fbp(A_n,blr(:),theta,filter)*(n/dsz)^2,n,n)/pixelsize;
xlr(mask) = 0;

% Contrast-to-noise
sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);
sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);
sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);
sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])
if i == 4
 % plot CNR ROI    
 imagesc(ax1,x)

 imagesc(ax2,xpm)
 imagesc(ax3,xpv)
 imagesc(ax4,xlr)
  hold on 
rectangle(ax4,'Position',[ys,xs,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','-')
rectangle(ax4,'Position',[yb2,xb2,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','--')
hold off 



else
 imagesc(ax1,x)
 imagesc(ax2,xpm)
 imagesc(ax3,xpv)
 imagesc(ax4,xlr)
end

% plot settings 
caxis(ax1,[cmin cmax])
set(ax1, 'XTick', [], 'YTick', []);
caxis(ax2,[cmin cmax])
set(ax2, 'XTick', [], 'YTick', []);
caxis(ax3,[cmin cmax])
set(ax3, 'XTick', [], 'YTick', []);
caxis(ax4,[cmin cmax])
set(ax4, 'XTick', [], 'YTick', []);
if i == 1
    ylabel(ax1,'1')
elseif i == 2
    ylabel(ax1,'2')
elseif i ==3
    ylabel(ax1,'4')
else
    ylabel(ax1,'8')
    title(ax1,'FBP')
    title(ax2,'WF-FBP')
    title(ax3,'NLM-FBP')
    title(ax4,'LR-FBP')
    g = colorbar('location','eastoutside','Position',[0.89 ystart(4) 0.02 0.83],'fontsize',20);
    g.Label.Interpreter = 'latex';
    g.Label.String = '$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]';
end



end
if cmap == 1
    colormap inferno
end
saveas(gcf,fullfile(graphicpath,'Experiment2FBP'),'epsc')

%% Figure 9: Number of flat-fields TV
% plot settings 
flatcount = [1 2 4 8];
ystart    = [0.72, 0.51, 0.30, 0.09];
xstart    = [0.05 0.26 0.47 0.68];
cmin      = -0.3;
cmax      = 2.5;
dim       = 0.2;
figure('Position',[100 100 1000 1000]) 


for i=1:length(flatcount)
% Compute SVD approximation     
[U,S,V]   = svd(reshape(flat(1:flatcount(i),:,:),flatcount(i)*r,k),'econ');
flatappr2 = reshape(U(:,1)* S(1,1) * V(:,1)',flatcount(i),r,k);

ax1 = axes('Position',[xstart(1) ystart(5-i) dim dim]);
ax2 = axes('Position',[xstart(2) ystart(5-i) dim dim]);
ax3 = axes('Position',[xstart(3) ystart(5-i) dim dim]); 
ax4 = axes('Position',[xstart(4) ystart(5-i) dim dim]); 
energy   = 300; % 30, 90

F        = squeeze(mean(flat(1:flatcount(i),:,energy),[1 3]))';
Y        = squeeze(mean(sino(:,:,energy),3))';
Flr      = squeeze(mean(flatappr2(:,:,energy),[1 3]))';

clear options;

options = struct(...
    'u0',zeros(n*n,1),...
    'rho',1.8,...
    'tau',tau,...
    'lambda',lamb,...
    'tolf',1e-4,... 
    'maxiters',maxiters,...
    'uhold',[],...
    'mask',true,...
    'verbose',0);


xlr = gd_wls(A_n,Flr,Y,options)/pixelsize; % LR-TV 
x   = gd_wls(A_n,F,Y,options)/pixelsize;   % TV 
options.preproc = 1;
xpm  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV
options.preproc = 2;
xpv  = gd_wls(A_n,F,Y,options)/pixelsize;   % P-TV

x   = reshape(x,n,n);
xpm  = reshape(xpm,n,n);
xpv = reshape(xpv,n,n);
xlr = reshape(xlr,n,n);  

 

sROI = x(xs:xs+60,ys:ys+60);
bROI = x(xb2:xb2+60,yb2:yb2+60);
t    = round(constrastnoise(sROI,bROI),2);

sROI = xpm(xs:xs+60,ys:ys+60);
bROI = xpm(xb2:xb2+60,yb2:yb2+60);
tpm   = round(constrastnoise(sROI,bROI),2);

sROI = xpv(xs:xs+60,ys:ys+60);
bROI = xpv(xb2:xb2+60,yb2:yb2+60);
tpv   = round(constrastnoise(sROI,bROI),2);

sROI = xlr(xs:xs+60,ys:ys+60);
bROI = xlr(xb2:xb2+60,yb2:yb2+60);
tlr  = round(constrastnoise(sROI,bROI),2);
display([t tpm tpv tlr])
if i == 4
 % plot CNR ROI    
 imagesc(ax1,x)

 imagesc(ax2,xpm)
 imagesc(ax3,xpv)
 imagesc(ax4,xlr)
  hold on 
rectangle(ax4,'Position',[ys,xs,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','-')
rectangle(ax4,'Position',[yb2,xb2,60,60],...
  'EdgeColor', 'white',...
  'LineWidth', 2,...
  'LineStyle','--')
hold off 



else
 imagesc(ax1,x)
 imagesc(ax2,xpm)
 imagesc(ax3,xpv)
 imagesc(ax4,xlr)
end

% plot settings 
caxis(ax1,[cmin cmax])
set(ax1, 'XTick', [], 'YTick', []);
caxis(ax2,[cmin cmax])
set(ax2, 'XTick', [], 'YTick', []);
caxis(ax3,[cmin cmax])
set(ax3, 'XTick', [], 'YTick', []);
caxis(ax4,[cmin cmax])
set(ax4, 'XTick', [], 'YTick', []);
if i == 1
    ylabel(ax1,'1')
elseif i == 2
    ylabel(ax1,'2')
elseif i ==3
    ylabel(ax1,'4')
else
    ylabel(ax1,'8')
    title(ax1,'TV')
    title(ax2,'WF-TV')
    title(ax3,'NLM-TV')
    title(ax4,'LR-TV')
    g = colorbar('location','eastoutside','Position',...
                [0.89 ystart(4) 0.02 0.83],'fontsize',20);
             
    g.Label.Interpreter = 'latex';
    g.Label.String = '$\Sigma_{\rm{tot}}(\lambda)$, [cm$^{-1}$]';
end

if cmap == 1
    colormap inferno
end

end

saveas(gcf,fullfile(graphicpath,'Experiment2TV'),'epsc')
