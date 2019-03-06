%% Enhancement of the 2D coherence slice
clear 
close all
clc
%add path
addpath('coherencefilter_version5b\')
addpath('phasemap_v1.1\phasemap\')
% load input image
I = (imread('fundus2D.png'));            %fundus2D.png    1.jpg  cameraman.gif
if size(I,3) == 3
    I = rgb2gray(I);
end
load my
%I = my;
figure,imshow(I);
title('原图')
%%
%扩散滤波
IC = CoherenceFilter(I,struct('T',2,'rho',10,'Scheme','N','verbose','none'));
IC = KuwaharaFourier(IC,7);
V2 = vesselness2D(IC, 0.5:0.2:3, [1;1], 1, false);%最重要的函数！！
%计算方位角（根据Hessian矩阵特征向量的几何意义）
theta = [-90:10:90];
I = V2;
I = medfilt2(I);
[Dxx,Dxy,Dyy] = Hessian2D(I,1);
[~,~,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);
thetaI = atand(Iy./Ix);
thetaI(isnan(thetaI)) = 90;
[row,col] = size(I);
J0 = steerGaussFilterOrder2(I,0,2,false);
J90 = steerGaussFilterOrder2(I,90,2,false);
J0(J0 > 0) = 0; J0 = J0 / min(J0(:));
J90(J90 > 0) = 0; J90 = J90 / min(J90(:));
J = sqrt(J0.^2 + J90.^2);
%带方向约束，去除倾角为thetaRM的裂缝
thetaRM = -60;
J2 = steerGaussFilterOrder2(I,thetaRM,2,false); 
J2(J2 > 0) = 0; J2 = J2 / min(J2(:));
%不带方向约束
J2 = zeros(size(J2));
J = J - J2;
J (J < 0) = 0;
J = J/max(J(:));
%% 绘图
thetaIT = medfilt2(thetaI,[5,5]);
thetaIT = (thetaIT+90)/(180);
HLS(:,:,1) = thetaIT;
HLS(:,:,2) = ones(size(J));
HLS(:,:,3) = 1-J;
RGB = hsl2rgb(HLS);
figure,imshow(RGB)
title('裂缝与方位角融合显示')
colormap('hsv')
colorbar,caxis([-90,90])
% hold on
% phasemap
% phasebar('location','sw','size',0.3,'deg')
%% 画theta 矢量图
thetaIT = medfilt2(thetaI,[5,5]);
V2x = cosd(thetaIT);
V2y = sind(thetaIT);
dh = 10;
[ny, nx] = size(thetaIT);
figure,imshow(I),alpha(0.7),hold on
h=quiver(1:dh:nx,ny:-dh:1,V2x(ny:-dh:1,1:dh:nx),V2y(ny:-dh:1,1:dh:nx),'b');
view(0,90);
title('方位角矢量图')