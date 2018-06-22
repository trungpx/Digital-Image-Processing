
    %% How to run: Just click "Run" button and hit "Select folder"
    
function FFT2D_DCT2D_DHT_WAVELET_transform()
    clc;
    clear;
    close all;
    imgdir = uigetdir('Test_images');
    %% read gray noise image
    file = fopen(fullfile(imgdir,'\Gray_monarch_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);
    %% Display gray noise image
    figure; imshow(gray_image,[]); title('Problem 1) Origional');
    gray_image_1 = double(gray_image);
    %% Calculate forward FFT
    fprintf('Calculating for FFT forward...\n');
    tic % calculate time elapsed
    result = FFT2D(gray_image_1); %using the FFT2D as defined below
    % "log" method for magnitude
    mag = log(1+abs(result));
    phase = angle(result);
    % Display FFT image
    figure; imshow(mag,[]); title('Magnitude FFT');
    figure; imshow(phase,[]); title('Phase FFT');
    toc; %finished time calculation
    %% Calculate inverse FFT
    fprintf('Calculating for IFFT...\n');
    tic;
    inverse = real(IFFT2D(result));
    % Display inverse FFT image
    figure; imshow(inverse,[]); title('IFFT');
    toc;
    %% FFT for only magnitude of transformed result
    fprintf('Calculating for...\n');
    tic;
    inverse_mag = real(IFFT2D(abs(result)));
    figure; imshow(inverse_mag,[]); title('IFFT for maginitude');
    toc;
    %% FFT for only phase of transformed result
    fprintf('Calculating for...\n');
    tic;
    inverse_phase = real(IFFT2D(angle(result)));
    figure; imshow(inverse_phase,[]); title('IFFT for phase');
    toc;
    %% DCT forward
    fprintf('Calculating for DCT...\n');
    tic;
    DCT_result = myDCT2D(gray_image_1);
    figure; imshow(log(1+abs(DCT_result)),[]); title('DCT');
    toc;
    %% DCT inverse
    fprintf('Calculating for IDCT...\n');
    tic;
    DCT_inv = myIDCT2D(DCT_result);
    figure; imshow(DCT_inv,[]); title('IDCT');
    toc;
    %% DCT 16x16 blocks
    fprintf('Calculating for 16x16 block DCT...\n');
    tic;
    DCT_16 = DCT16(gray_image_1);
    figure; imshow(log(1+abs(DCT_16)),[]); title('16x16 block DCT');
    toc;
    %% IDCT 16x16 blocks
    fprintf('Calculating for 16x16 block IDCT...\n');
    tic;
    DCT_16_inv = IDCT16(DCT_16);
    figure; imshow(DCT_16_inv,[]); title('16x16 block IDCT');
    toc;
    %% Discrete Hadamard transform (DHT)
    fprintf('Calculating for DHT...\n');
    tic;
    DHT = myDHT2D(gray_image_1,'T');
    figure; imshow(log(1+abs(DHT)),[]); title('DHT');
    toc;
    %% Discrete Hadamard transform inverse (IDHT)
    fprintf('Calculating for IDHT...\n');
    tic;
    DHT_inv = myDHT2D(DHT,'i');
    figure; imshow(DHT_inv,[]); title('IDHT');
    toc;
    %% Wavelet Transform
    fprintf('Calculating for Wavelet Transform...\n');
    tic;
    DWT = myDWT2D(gray_image_1);
    figure; imshow(log(1+abs(DWT)),[]); title('Wavelet Transform');
    toc;
end

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
% function to get exactly the output is row or column
function y = myFFT(x)
    y = myFFT1(x);    
    if (isrow(x)~= 0)
        y = y';
    end
end
%code for incursive FFT algorithm
function y = myFFT1(x)
    x = x(:);
    N = length(x);
    even = x(1:2:end); % get the even part
    odd = x(2:2:end); % get the odd part
    if (N >= 8) %perform the incurve algorithm
        x_even = myFFT1(even);
        x_odd = myFFT1(odd);
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
     M = length(x);
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
%% function to calculate DCT16x16
function y = DCT16(x)
    [M,N] = size(x);
    y = zeros(M,N);
    for i=1:16:M
        for j=1:16:N %jump to each block 16x16
            block = myDCT2D(x(i:i+15,j:j+15));
            y(i:i+15,j:j+15) = block; % give the output the value when apply DCT for each block
        end
    end
end
%% function to calculate DCT16x16
function y = IDCT16(x)
    [M,N] = size(x);
    y = zeros(M,N);
    for i=1:16:M
        for j=1:16:N %jump to each block 16x16
            block = myIDCT2D(x(i:i+15,j:j+15));
            y(i:i+15,j:j+15) = block; % give the output the value when apply DCT for each block
        end
    end
end

%% function to calculate DHT & IDHT
function y = myDHT2D(D,Type)
    % D: 2D input matrix 
    % Type: 'T' = Transform , or 'i' = Inverse Trasnform
    n = size(D,1);
    H = Hadamard_matrix(n);
    %Forward Transfomation -----------
    if (Type=='T')
        y = H*D;
        y = y*H;
    end
    %Inverse Transfomation -----------
    if (Type=='i')
        y = D*H^-1;
        y = H^-1*y;
    end
    %--------------------------------------
end

%% Build H matrix
function H = Hadamard_matrix(n)
    [f,e] = log2([n n/12 n/20]);
    k = find(f==1/2 & e>0);
    e = e(k)-1;
if k == 1 %N = 1 * 2^e;
    H = ones(); % initiate H with 1 value
elseif k == 2 %N = 12 * 2^e;
    H = [ones(1,12); ones(11,1) toeplitz([-1 -1 1 -1 -1 -1 1 1 1 -1 1],[-1 1 -1 1 1 1 -1 -1 -1 1 -1])]; %create toeplitz matrix
elseif k == 3 %N = 20 * 2^e;
    H = [ones(1,20); ones(19,1) hankel([-1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1 1], [1 -1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1])];
end
%  Kronecker creation
for i = 1:e
    H = [H  H; H -H];
end
end

%% Wavelet transform
function y = myDWT2D(x)
    % Inputs: x: 2D data array (e.g., image) 
    h_matrix=[0.683 1.183 0.317 -0.183];
    K = length(h_matrix);
    [M,N]=size(x); 
    h_matrix=h_matrix/norm(h_matrix);% normalize h
    h0_matrix=zeros(1,N); h0_matrix(1:K)=h_matrix;
    H0 = (myFFT(h0_matrix));
    for k=0:N-1
       m=mod(k-N/2,N)+1;
       H1(k+1)=-exp(-1i*2*pi*k/N)*conj(H0(m)); % perform wavelet transform as formula
    end
    a=x; 
    y=zeros(N);
    n=length(a);
    while n>64 % level 3 (512/64 = 8 = 2^3)
        temp = zeros(n,n); % init output matrix to zeros
        for k=1:n          % Go for all n columns
            A = myFFT(double(a(:,k)));
            D = real(IFFT1D(A.*H1'));  % Perform as in freqency domain
            A = real(IFFT1D(A.*H0'));  % Perform as in freqency domain
            temp(:,k) = [A(2:2:n); D(2:2:n)];
        end
        for k=1:n                  % Go for all n rows
            A = (myFFT(double(temp(k,:))));
            D = (real(IFFT1D(A.*H1)));  % convolution as in freqency domain
            A = (real(IFFT1D(A.*H0)));  % convolution as in freqency domain
        temp(k,:)=[A(2:2:n) D(2:2:n)];
        end
        y(1:n,1:n)=temp;              % Perform concatenate wavelet coefficients
        H0=H0(1:2:length(H0));     % Perform subsampling H0
        H1=H1(1:2:length(H1));     % Perform subsampling H1
        n=n/2;
        a=temp(1:n,1:n);     
    end
end

