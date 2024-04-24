% Script file to run uncertainty quantification 
% Calculated parameters were [mG, sG, bx, by] 
clear all;
close all;

p=0; % exponential
Nsim = 50; % Monte carlo iterations
X1 = FEMmodel(); % Initialize the object
Lmax=100;
X1.Lx=Lmax;
X1.Ly=X1.Lx;
    
%% Choose posterior results:
C11paramSelected = [ 1.55, 0.129, 20.18, 16.6];
% C11paramSelected = [ 1.55, 0, 20.18, 16.6];

C33paramSelected = [ 0.316, 0.027, 20.86, 15.08];   
% C33paramSelected = [ 0.316, 0, 20.86, 15.08];   

% Output file
fileName = ['output_Disp_Plate_L=' num2str(sprintf('%.0f',X1.Lx)) 'x' num2str(sprintf('%.0f',X1.Ly))  '.txt'];
fileID = fopen(fileName, 'w');
  
% element size from min corr. length
minbx = min(C11paramSelected(3),C33paramSelected(3));
minby = min(C11paramSelected(4),C33paramSelected(4));
b = min(minbx, minby);
Delem = b/4;
X1.nx = round(X1.Lx/Delem)+1; % Number of nodes in X
X1.ny = X1.nx; % Number of nodes in Y
nelX = X1.nx-1; % Number of elements in x-dir
nelY = X1.ny-1; % Number of elements in y-dir
nel = nelX*nelY; % Total number of elements

%% Computation of ksix and ksiy matrices
xMtrx = repmat((0:nelX-1)*Delem,nelY,1);
yMtrx = repmat((0:nelY-1)'*Delem,1,nelX);
xMtrxTemp = xMtrx(:);
yMtrxTemp = yMtrx(:);
numMeasuremnts = length(xMtrxTemp);
% numMeasuremnts = numRVEs;
ksix=zeros(numMeasuremnts,numMeasuremnts);
ksiy=zeros(numMeasuremnts,numMeasuremnts);
% distMtrx(numMeasuremnts,numMeasuremnts)=0;

for im=1:numMeasuremnts
    ksix(im,:)=xMtrxTemp(:)-xMtrxTemp(im);
    ksiy(im,:)=yMtrxTemp(:)-yMtrxTemp(im);
end

%%%%%%%%%%%%%%% COMPUTE CORRELATION MATRICES %%%%%%%%%%%%%%%%%%%

% C11
muGC11 = C11paramSelected(1);
sigmaGC11 = C11paramSelected(2);
bxC11 = C11paramSelected(3);
byC11 = C11paramSelected(4);
% r = sqrt(ksix.^2 + ksiy.^2);
r = sqrt((ksix/bxC11).^2 + (ksiy/byC11).^2);
% corMtrx2DC11 = exp(-r); %exponential correlation matrix
% corMtrx2DC11 = sigmaGC11^2*exp(-r); %exponential correlation matrix
corMtrx2DC11 = exp(-r.^2); %square exponential correlation matrix

% C33
muGC33 = C33paramSelected(1);
sigmaGC33 = C33paramSelected(2);
bxC33 = C33paramSelected(3);
byC33 = C33paramSelected(4);
r = sqrt((ksix/bxC33).^2 + (ksiy/byC33).^2);
% corMtrx2DC33 = exp(-r); %exponential correlation matrix
% corMtrx2DC33 = sigmaGC33^2*exp(-r); %exponential correlation matrix
corMtrx2DC33 = exp(-r.^2); %square exponential correlation matrix




%%%%%%%%%%%%%%%%%%%% COVARIANCE DECOMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%
m=numMeasuremnts;

[VC11,DC11] = eigs(corMtrx2DC11,m); %diagonal matrix D containing the eigenvalues on the main diagonal, and matrix V whose columns are the corresponding eigenvectors
[VC33,DC33] = eigs(corMtrx2DC33,m);

% Thresholding negative eigenvalues to zero
DC11(DC11 < 0) = 0;
DC33(DC33 < 0) = 0;

sigma_hat_C11=  trace(DC11)/numMeasuremnts;
sigma_hat_C33=  trace(DC33)/numMeasuremnts;
%NEW_COV= VC11*DC11*VC11';          %new field is approximately homogeneous

m = size(DC11,1);
UC11 = randn(m,Nsim);
UC33 = randn(m,Nsim);
 % From eq. 3.59 Lecture Notes Iason
C11weights = sqrt(diag(DC11)).*UC11;
C33weights = sqrt(diag(DC33)).*UC33;
XC11 = VC11*C11weights; % X is the same as U(omega,x) in eq. 26 CMAME2020 the underlying Gaussian field 
XC33 = VC33*C33weights;

C11 = exp(muGC11 +sigmaGC11*XC11); % This is A(omega,x) in eq. 26 CMAME2020 the LogNormal homogeneous translation field
C33 = exp(muGC33 +sigmaGC33*XC33);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize  (Remove default loading and boundary conditions)
X1.BCedge=[]; 
X1.BCpoint=[];
X1.Fedge=[];
X1.Fpoint=[];
X1 = X1.mesh2d();

% figure
% X1.plotmesh(); drawnow;

%% plot lambda and mu
x = linspace(0, Lmax, nelX);
y = linspace(0, Lmax, nelY);

lambda_ = C11 - 2*C33;
mu_ = C33;

figure;
for i = 1:10
    subplot(2, 5, i);
    imagesc(x, y, reshape(lambda_(:,i), nelX,nelY))
    set(gca,'YDir','normal')
    title(['\lambda, sample ' num2str(i)]);
    colorbar;
    colormap(subplot(2, 5, i), 'copper');
end
lambda_mean = mean(lambda_,2);
% lambda_mean = lambda_(:,1);
lambda_std = std(lambda_,0,2);
mu_mean = mean(mu_,2);
mu_std = std(mu_,0,2);

figure;
subplot(1, 2, 1);
imagesc(x, y, reshape(lambda_mean, nelX,nelY))
set(gca,'YDir','normal')

title('Mean of \lambda');
colorbar;
subplot(1, 2, 2);
imagesc(x, y, reshape(lambda_std, nelX,nelY))
set(gca,'YDir','normal')
title('Std of \lambda');
colorbar;

colormap(subplot(1, 2, 1), 'copper'); % You can change 'parula' to any colormap you prefer
colormap(subplot(1, 2, 2), 'copper'); % You can change 'parula' to any colormap you prefer

%%
% Apply bcs
X1.BCedge  = [1 0 1 0 0 0]; % [direction value fix-x fix-y disp-x disp-y]
X1.BCedge  = [X1.BCedge;2 0 0 1 0 0]; % [direction value fix-x fix-y disp-x disp-y]
X1.BCpoint = [0 0 1 1 0 0];% [x-val y-val fix-x fix-y disp-x disp-y] 
pstress = 0.1;
Area = X1.Ly; % thickness is one
xforce = pstress*Area;
X1.Fedge   = [1  X1.Lx 1 xforce];% [direction value force-direction force-value]
monitorNode = find(X1.node(:,1)==X1.Lx & X1.node(:,2)==X1.Ly);


%% Perform the analysis
X1 = X1.bcforce();

% MCS
% MCS
ux_result = cell(Nsim,1);
uy_result = cell(Nsim,1);
exx_result = cell(Nsim,1);
eyy_result = cell(Nsim,1);
exy_result = cell(Nsim,1);
sxx_result = cell(Nsim,1);
syy_result = cell(Nsim,1);
sxy_result = cell(Nsim,1);
for isim = 1:Nsim
    
    X1 = X1.assembleKF(C11(:,isim), C33(:,isim));
    X1 = X1.solver();
    X1 = X1.writeOutput(monitorNode,isim,fileID);
    if rem(isim,100)==0
        fprintf('%d\n',isim)
    end
    ux_result{isim} = X1.u(:,1);
    uy_result{isim} = X1.u(:,2);
end

%% Plot mean and std heatmap of displacement

% Convert cell matrices to numeric arrays
ux_data = squeeze(cat(2, ux_result{:})); % Combine the 1x784 matrices into a 2D matrix
uy_data = squeeze(cat(3, uy_result{:}));

% Calculate the mean and standard deviation
mean_ux = mean(ux_data, 2); % Calculate the mean along rows
std_ux = std(ux_data, 0, 2); % Calculate the standard deviation along rows
mean_uy = mean(uy_data, 2);
std_uy = std(uy_data, 0, 2);

% Create a square domain for plotting
% Create a square domain for plotting
n = sqrt(size(mean_ux, 1)); % Assuming it's a square domain
x = linspace(0, Lmax, n);
y = linspace(0, Lmax, n);

% Plot mean and std as heatmaps
figure;
subplot(2, 2, 1);
imagesc(x, y, reshape(mean_ux, n, n));
set(gca,'YDir','normal')
title('Mean of u_x');
colorbar;

subplot(2, 2, 2);
imagesc(x, y, reshape(std_ux, n, n));
set(gca,'YDir','normal')
title('Standard Deviation of u_x');
colorbar;

subplot(2, 2, 3);
imagesc(x, y, reshape(mean_uy, n, n));
set(gca,'YDir','normal')
title('Mean of u_y');
colorbar;

subplot(2, 2, 4);
imagesc(x, y, reshape(std_uy, n, n));
set(gca,'YDir','normal')
title('Standard Deviation of u_y');
colorbar;

colormap(subplot(2, 2, 1), 'copper'); % You can change 'parula' to any colormap you prefer
colormap(subplot(2, 2, 3), 'copper'); % You can change 'parula' to any colormap you prefer
colormap(subplot(2, 2, 2), 'bone'); % You can change 'parula' to any colormap you prefer
colormap(subplot(2, 2, 4), 'bone'); % You can change 'parula' to any colormap you prefer


%% Save the data

save('Ux_samples.mat', 'ux_data');
save('Uy_samples.mat', 'uy_data');