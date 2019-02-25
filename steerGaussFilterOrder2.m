function J = steerGaussFilterOrder2(I,theta,sigma,disFlag)

%    This function implments the steerable filter of the    
%    second derivative of Gaussian function (X-Y separable version)
%    of the paper
% 
%    W. T. Freeman and E. H. Adelson, "The Design
%     and Use of Steerable Filters", IEEE PAMI, 1991.
%
%    J = steerGaussFilterOrder2(I,theta,sigma,disFlag) 
%    calculated 
%    Input:
%    1. I: input image
%    2. theta: the orientation
%    3. sigma: standard deviation of the Gaussian template   
%    Output:
%    J. The response of derivative in theta direction
%
%    Author: Jincheng Pang, Tufts University, Dec. 2013.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I: Assign algorithm parameters.

I = mean(double(I),3);

% Process input arguments (according to usage mode).
if ~exist('arg1','var') || ~isstruct(arg1)
   
    
    % Assign default filter orientation (if not provided).
    if ~exist('theta','var') || isempty(theta)
      theta = 0;
    end
    theta = -theta*(pi/180);

    % Assign default standard deviation (if not provided).
    if ~exist('sigma','var') || isempty(sigma)
       sigma = 1;
    end
    
    % Assign default visualization state (if not provided).
    if ~exist('disFlag','var') || isempty(disFlag)
       disFlag = false;
    end   
    
end % End of input pre-processing.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II: Evaluate separable filter kernels.

% Determine necessary filter support (for Gaussian).
Wx = floor((8/2)*sigma); 
if Wx < 1
  Wx = 1;
end
x = [-Wx:Wx];

[xx,yy] = meshgrid(x,x);

g0 = exp(-(xx.^2+yy.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));
% G2a = 0.9213*(2*xx.^2-1).*g0;
% G2b =  1.843*xx.*yy.*g0;
% G2c = 0.9213*(2*yy.^2-1).*g0;
G2a=(1*xx.^2/sigma^4-1/sigma^2).*g0;   %各种二阶的幅度系数
G2b=xx.*yy/sigma^4.*g0;
G2c=(1*yy.^2/sigma^4-1/sigma^2).*g0;
% H2a=(-2.25*xx/sigma^4+xx.^3/sigma^6).*g0;
% H2b=(-0.7515/sigma^4+xx.^2/sigma^6).*g0;
% H2c=(-0.7515/sigma^4+yy.^2/sigma^6).*g0;
% H2d=(-2.254*yy/sigma^4+yy.^3/sigma^6).*g0;
% kh2a = cos(theta)^3;
% kh2b=-3*cos(theta)^2*sin(theta);
% kh2c=3*cos(theta)*sin(theta)^2;
% kh2d=-sin(theta)^3;
% G4a=(0.75/sigma^4-3*xx.^2/sigma^6+xx.^4/sigma^8).*g0;
% G4b=(-1.5*xx/sigma^6+xx.^3/sigma^8).*g0;
% G4c=(xx.^2/sigma^6-0.5/sigma^4).*yy/sigma^6.*g0;
% G4d=(-1.5*yy/sigma^6+yy.^3/sigma^8).*xx.*g0;
% G4e=(0.75/sigma^4-3*yy.^2/sigma^6+yy.^4/sigma^8).*g0;
% kg4a=cos(theta)^4;
% kg4b=-4*cos(theta)^3*sin(theta);
% kg4c=6*cos(theta)^2*sin(theta)^2;
% kg4d=-4*cos(theta)*sin(theta)^3;
% kg4e=sin(theta)^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Determine oriented filter response.

% Calculate image gradients (using separability).
I2a = imfilter(I,G2a,'same','replicate');
I2b = imfilter(I,G2b,'same','replicate');
I2c = imfilter(I,G2c,'same','replicate');
% I2a = imfilter(I,H2a,'same','replicate');
% I2b = imfilter(I,H2b,'same','replicate');
% I2c = imfilter(I,H2c,'same','replicate');
% I2d = imfilter(I,H2d,'same','replicate');
% I4a = imfilter(I,G4a,'same','replicate');
% I4b = imfilter(I,G4b,'same','replicate');
% I4c = imfilter(I,G4c,'same','replicate');
% I4d = imfilter(I,G4d,'same','replicate');
% I4e = imfilter(I,G4e,'same','replicate');
% Evaluate oriented filter response.
J = (cos(theta))^2*I2a+sin(theta)^2*I2c-2*cos(theta)*sin(theta)*I2b;
%J1 = J(floor(size(I,1)/2),floor(size(I,2)/2));

% J = I2a*kh2a+I2b*kh2b+I2c*kh2c+I2d*kh2d;
% F = H2a*kh2a+H2b*kh2b+H2c*kh2c+H2d*kh2d;
% J = I4a*kg4a+I4b*kg4b+I4c*kg4c+I4d*kg4d+I4e*kg4e;
% F = G4a*kg4a+G4b*kg4b+G4c*kg4c+G4d*kg4d+G4e*kg4e;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV: Visualization

% Create figure window and display results.
% Note: Will only create output window is requested by user.
if disFlag
   fig = figure(1); clf; set(gcf,'Name','Oriented Filtering');
   subplot(1,3,1); imagesc(I); axis image; colormap(gray);
      title('Input Image');
   subplot(1,3,2); imagesc(J); axis image; colormap(gray);
      title(['Filtered Image (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
   F = (cos(theta))^2*G2a+sin(theta)^2*G2c-2*cos(theta)*sin(theta)*G2b;
   subplot(1,3,3); imagesc(F);%imagesc((cos(theta))^2*G2a+sin(theta)^2*G2c-2*cos(theta)*sin(theta)*G2b);
      axis image; colormap(gray);
      title(['Oriented Filter (\theta = ',num2str(-theta*(180/pi)),'{\circ})']);
   saveas(fig,[int2str(int8(fig.Number)),'.jpg']);
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%