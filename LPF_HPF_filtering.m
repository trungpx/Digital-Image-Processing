
    %% How to run: Just click "Run" button and hit "Select folder"
    
function LPF_HPF_filtering()   
    clear variables; close all; clc;
	%Read the gray image
    imgdir = uigetdir('Test_images');
    file = fopen(fullfile(imgdir,'Gray_monarch_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    figure; imshow(gray_image,[]); title('Gray\_monarch\_512x512.raw');
    fclose(file);
    %% Transform Domain Filtering
    gray_image = double(gray_image);
    DFT_image = FFT2D(gray_image);
    [M,N] = size(gray_image);
    disp('Start...');
	tic;
    %..LPF.. filtering
    disp('(1) No noise, LPF, DFT...');
    LPF = LPF_add(M,N);
    variable = LPF.*DFT_image;
    filter_image = real(IFFT2D(double(variable)));
    figure; imshow(filter_image,[]); title('No noise, LPF, DFT');
    %% ---- DFT start ----
    %..HPF.. filtering
    disp('(2) No noise, HPF, DFT...');
    HPF = HPF_add(M,N);
    variable = HPF.*DFT_image;
    filter_image = real(IFFT2D(double(variable)));
    figure; imshow(filter_image,[]); title('No noise, HPF, DFT');
    %% Noise added
    noise_image = addGaussian(gray_image);
    noise_image_DFT = FFT2D(noise_image);
    %LPF
    disp('(3) Gaussian noise, LPF, DFT...');
    variable = LPF.*noise_image_DFT;
    filter_image = real(IFFT2D(double(variable)));
    figure; imshow(filter_image,[]); title('Gaussian noise, LPF, DFT');
    %..HPF
    disp('(4) Gaussian noise, HPF, DFT...');
    variable = HPF.*noise_image_DFT;
    filter_image = real(IFFT2D(double(variable)));
    figure; imshow(filter_image,[]); title('Gaussian noise, HPF, DFT');
    %% DFT finished

    %% DCT start
    % LPF part
    disp('(5) No noise, LPF, DCT...');
    DCT_image = myDCT2D(gray_image);
    variable = LPF.*DCT_image;
    filter_image = myIDCT2D(double(variable));
    figure; imshow(filter_image,[]); title('No noise, LPF, DCT');
    %...HPF
    disp('(6) No noise, HPF, DCT...');
    variable = HPF.*DCT_image;
    filter_image = myIDCT2D(double(variable));
    figure; imshow(filter_image,[]); title('No noise, HPF, DCT');
    %% Gaussian noise added
    noise_image_DCT = myDCT2D(noise_image);
    %..LPF
    disp('(7) Gaussian noise, LPF, DCT...');
    variable = LPF.*noise_image_DCT;
    filter_image = myIDCT2D(double(variable));
    figure; imshow(filter_image,[]); title('Gaussian noise, LPF, DCT');
    %..HPF
    disp('(8) Gaussian noise, HPF, DFT...');
    variable = HPF.*noise_image_DCT;
    filter_image = myIDCT2D(double(variable));
    figure; imshow(filter_image,[]); title('Gaussian noise, HPF, DCT');
	toc;
    disp('Completed!');

end

%%  Low pass filter
function y = LPF_add(M,N) % M,N is size of filter
        F = 80; %The cut off frequency initialization 
        a = 0:(M-1);
        b = 0:(N-1);
        index = find(a>M/2);
        a(index) = a(index)-M;
        indey = find(b>N/2);
        b(indey) = b(indey)-N;
        [B,A] = meshgrid(b,a);
        T = sqrt(A.^2+B.^2);
        y = double(T<=F); % Cut-off frequency comparison
end

function y = HPF_add(M,N) % arguments are size of filter
        LPF = LPF_add(M,N);
        y = 1 - LPF;
end

function y = addGaussian (x)
    SNR_dB = 5;
    SNR = 10^(SNR_dB/10); % Change SNR to linear scale
    [M,L] = size(x);
    rng('default');% set the random generator seed (for comparison)
    y = zeros(M,L);
    for i=1:M
        Esym = sum(abs(x(i,:)).^2)/(L); %Calculate actual symbol energy
        N0 = Esym/SNR; %Find the noise spectral density
        sigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
        n = sigma*randn(1,L);%compute noise
        y(i,:) = x(i,:) + n; %output y
    end
end

%% FFT2D functions
% Function for FFT with 2D dimension
function y = FFT2D(x)
y = FFT_matrix(FFT_matrix(x).').'; % perform FFT for the vector that is column of the matrix x
end

function y = FFT_matrix(x)
    y = zeros(size(x,1),1);
    for i = 1:size(x,2)
        y = [y myFFT(x(:,i))];
    end
    y = y(:,2:end);
end

%code for incurve FFT algorithm
function y = myFFT(x)
    x = x(:);
    N = length(x);
    even = x(1:2:end); % get the even part
    odd = x(2:2:end); % get the odd part
    if (N >= 8) %perform the incurve algorithm
        x_even = myFFT(even);
        x_odd = myFFT(odd);
        factor = exp(-1i*2*pi*((0:N/2-1)')/N);
        z = factor .* x_odd;
        y = [(x_even + z);(x_even -z)];
    else
        switch N
            case 2
                y = [1 1;1 -1]*x; % if N = 2 then calculate instantly
            case 4
                y = [1 0 1 0;0 1 0 -1i;1 0 -1 0;0 1 0 1i]*[1 0 1 0;1 0 -1 0;0 1 0 1;0 1 0 -1]*x; % if N = 4 then calculate instantly
            otherwise
                error('The number of N must be 2^q');
        end
    end
end

function y = IFFT2D(x)
    [M,N] = size(x);
     y = conj(FFT2D(conj(x)))/(M*N); % using the formulation of inverse FFT in slide
end

function y = IFFT1D(x)
    M = size(x);
     y = conj(myFFT(conj(x)))/M; % using the formulation of inverse FFT in slide
end
%% function to calculate 2D DCT of an image
function y = myDCT2D(x)
    signal = double(x);
    N = size(signal,1);
    y = zeros(N);
    for k=1:N     %calculate 1D DCT of each row of image
        y(k,:) = myDCT1D(signal(k,:));  
    end
    for k=1:N     %calculate 1D DCT of each column of image
        y(:,k) = myDCT1D(y(:,k));
    end
end

%% function to calculate DCT of a 1D signal
function y = myDCT1D(signal)
    N = length(signal);
    signal = signal(:).';
    alpha1 = sqrt(2/N);
    for i=1:N
        if i==1   %for signal index of 1, alpha is 1/sqrt(l)
            alpha = sqrt(1/N);
        else
            alpha = alpha1;
        %for signal index of greater than 1
        end
        j=1:N;
        % summation calculates single entry of output by applying the  
        summ = sum(signal.*cos((pi*(2*(j-1)+1).*(i-1))/(2*N))); % formula of DCT on the signal
        y(i) = alpha*summ;
    end
end

%% function to calculate 2D IDCT of an image
function y = myIDCT2D(x)
y = IDCT_matrix(IDCT_matrix(x).').';
end
%function to calculate the IDCT column of a matrix
function y = IDCT_matrix(x)
    y = zeros(size(x,1),1);
    for i = 1:size(x,2)
        y = [y myIDCT1D(x(:,i)')]; 
    end
    y = y(:,2:end);
end

%% function to calculate IDCT of a 1D signal
function y = myIDCT1D(x)
    N = length(x);
    CN = zeros(N);
    alpha = sqrt(2/N);
    alpha1 = sqrt(1/N);
    for n=0:N-1
        for k=0:N-1
            if k==0
                CN(k+1,n+1)=alpha1;
            else
                CN(k+1,n+1)=alpha*cos(pi*(n+0.5)*k/N); %using the formula of IDCT
            end
        end
    end
    y = CN'*x'; %sum all result to get the output
end
