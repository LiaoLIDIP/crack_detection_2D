function filtered = KuwaharaFourier(original,winsize)
    %% Edge aware Kuwahara filter for 2D and 3D images: 
    % - black and white only (scalar values) 
    % - not suitable for 2D colored image (eg. RGB)
    %% Author Job G. Bouwman: jgbouwman@hotmail.com

    %% Incorrect input handling
    error(nargchk(2, 2, nargin, 'struct'));
    % non-double data will be cast
    if ~isa(original, 'double')
        original = double(original);
    end % if
    % wrong-sized kernel is an error
    if mod(winsize,2)~=1
        error([mfilename ':IncorrectWindowSize'],'winsize must be odd');
    end % if
    %% Build the square convolution kernel
    N_or            = size(original); 
    Dimensionality  = length(N_or);        
    Nsquare         = (winsize+1)/2;
    Nkernel         = repmat(Nsquare, [1, Dimensionality]); 
    Nrem            = winsize - Nsquare;
    squareKernel    = 1/Nsquare^Dimensionality*ones(Nkernel); 
    %% Convolution-wise calculation of averages and variances
    % NB: I have accelerated the convolution, by processing the averages and
    % the variances in a single convolution. This can be done by adding the 
    % squared image as a an imaginary passenger to the normal image.
    complexConvInput = original + 1i*original.^2; 
    % Convolute... 
    if prod(Nkernel) < 6*log(prod(N_or))  % This boundary could be improved 
        % ... (in case of small kernels) in the spatial domain:
        complexConvOutput = convn(complexConvInput, squareKernel);
    else
        % ... (in case of large kernels) in the Fourier domain: 
        % NB: Zero-padding and ensuring equal size of image and kernel
        complexConvInput = padarray(complexConvInput, Nkernel, 'post');  
        squareKernel     = padarray(squareKernel    , N_or   , 'post'); 
        % Multiplication in frequency domain:  
        complexConvOutput = ...
            ifftn(fftn(complexConvInput).*fftn(squareKernel));
    end 
    

%% Separate the real and imaginary part: 
    % Real part is the running average: 
    runningAvrge = real(complexConvOutput); 
    % imag part is the stds, only squared average must be yet subtracted: 
    runningStdev = imag(complexConvOutput)-runningAvrge.^2;  

%% Now select from the four quadrants (the eight octants) the right value
% - first create stack of four (eight) images
% - select minimal stdev
% - take the corresponding avarage 

%% Differerent implementations for 2D and 3D:  
    if Dimensionality == 2
        % Patch the North-East side, North-West side, ... onto a stack: 
        avg_Stack = cat(3,...
            runningAvrge(Nrem+(1:N_or(1)), Nrem+(1:N_or(2))),...
            runningAvrge(Nrem+(1:N_or(1)),       1:N_or(2)) ,... 
            runningAvrge(      1:N_or(1) , Nrem+(1:N_or(2))),...
            runningAvrge(      1:N_or(1) ,       1:N_or(2))); 
        std_Stack = cat(3,...
            runningStdev(Nrem+(1:N_or(1)), Nrem+(1:N_or(2))),...
            runningStdev(Nrem+(1:N_or(1)),       1:N_or(2)) ,... 
            runningStdev(      1:N_or(1) , Nrem+(1:N_or(2))),...
            runningStdev(      1:N_or(1) ,       1:N_or(2))); 

        % Choice of the index with minimum variance
        [minima,indices] = min(std_Stack,[],3); %#ok<ASGLU>
        
        % Selecting the accompanying average value: 
            [y,x] = meshgrid(1:N_or(2),1:N_or(1));
            lookupIndices = x+N_or(1)*(y-1)+prod(N_or)*(indices-1);
            filtered = avg_Stack(lookupIndices);
    else
        % 3D implementation:   
        % Combine the avgs and the stdev in a complex number.
        % (complete other trick than before) 
        % - stdev = absolute value
        % - avgs  -->  phase (normalize to [-pi, + pi]
        % This allows the min function to find minimum stdev
        % Take phase of that minimum, to get corresponding average. 
        av_Min = min(runningAvrge(:));
        av_Max = max(runningAvrge(:));
        av_Rng = av_Max - av_Min;
        normalizedAverages = (runningAvrge - av_Min)./av_Rng;
        complexPacking = runningStdev.*exp(1i*2*pi*normalizedAverages); 
        % Patch the North-East side, North-West side, ... onto a stack: 
        for dim = 1:3
            circShiftVector = zeros(1,3); 
            circShiftVector(dim) = -Nrem; 
            complexPacking = cat(4, ...
                complexPacking, ...
                circshift(complexPacking, circShiftVector)); 
        end
        % chop redundant part: 
        complexPacking = complexPacking(1:N_or(1),1:N_or(2),1:N_or(3),:);  
        % find minimum, and convert its phase into avgs again: 
        filtered = mod(angle(min(complexPacking,[],4)),2*pi)*(av_Rng/(2*pi)) + av_Min;
    end
end % of function
