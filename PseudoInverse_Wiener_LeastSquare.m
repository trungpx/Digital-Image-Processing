
%% How to run: Just click "Run" button and hit "Select folder"

function PseudoInverse_Wiener_LeastSquare()
close all;
clc;
imgdir = uigetdir('Test_images');
%% Loading directory, gray image
file = fopen(fullfile(imgdir,'\monarch_gray_512x512.raw'),'rb');
gray_image = fread(file,fliplr([512,512]),'*uint8')';
fclose(file); gray_image = double(gray_image);
%% Make LPF and add noise to image
disp('Image restoration');
disp('Blurring and adding noise...');
[M,N] = size(gray_image);
fft_H = LPF_add(M,N); % LPF window, bluring image
fft_ori_img = FFT2D(gray_image); % INPUT in frequency
LPF = real(IFFT2D(double(fft_H.*fft_ori_img)));
LPF_noise = add_gauss_noise(LPF); % add gausian noise
disp('Computing psnr and filtering...');

% Computes output img
out_inverse = m_pseudo_inv(LPF_noise);
out_wiener = m_wiener(LPF_noise);
lambda = 0.001;
out_least = m_least_sqr(LPF_noise,lambda);

% Computes psnr
psnr_LPF_noise = psnr(uint8(LPF_noise),uint8(LPF));
psnr_p_inverse = psnr(uint8(out_inverse),uint8(LPF));
psnr_wiener = psnr(uint8(out_wiener),uint8(LPF));
psnr_const = psnr(uint8(out_least),uint8(LPF));

% Displays psnr of each image
disp(['psnr_LPF_noise = ', num2str(psnr_LPF_noise), ' dB']);
disp(['psnr_p_inverse = ',num2str(psnr_p_inverse), ' dB']);
disp(['psnr_wiener = ', num2str(psnr_wiener), ' dB']);
disp(['psnr_const = ', num2str(psnr_const), ' dB']);

% Computes spectra images
disp('Drawing the power spectrum...');
power_ori_img = log(abs(fft_ori_img.^2) + eps);
power_noise_img = log(abs(FFT2D(double(LPF_noise)).^2) + eps);
power_inversed_img = log(abs(FFT2D(double(out_inverse)).^2) + eps);
power_wiener_img = log(abs(FFT2D(double(out_wiener)).^2) + eps);

% Displays figures
figure; imshow(gray_image,[]); title('Origional'); % show original image
figure; imshow(real(out_inverse),[]); title('image-pseudo-inverse-filtered');
figure; imshow(real(out_wiener),[]); title('image-wiener-filtered-filterred');
figure; imshow(real(out_least),[]); title('image-out-least-square');

% Plots the power spectrum
figure; clf; imagesc(power_ori_img), axis image; colorbar; title('Power spectrum of original image');
figure; clf; imagesc(power_noise_img), axis image; colorbar; title('Power spectrum of Gaussian noise image');
figure; clf; imagesc(power_inversed_img), axis image; colorbar; title('Power spectrum of inversed image');
figure; clf; imagesc(power_wiener_img), axis image; colorbar; title('Power spectrum of wiener image');

disp('Computing Completed!');
end
%% Pseudo-inverse filtered image function
function OUTPUT = m_pseudo_inv(INPUT)
% INPUT is the noise image
sigma = sqrt(4*10^(-5));
alpha = 0.5; % Indicates inverse filter with specific param
[M,N] = size(INPUT);
H_freq = LPF_add(M,N); % LPF window
X_freq = FFT2D(INPUT); % change INPUT to freq domain
powfreqx = X_freq.^2/(M*N);
G_freq = ((H_freq.')').*abs(powfreqx)./(abs(H_freq.^2).*abs(powfreqx)+alpha*sigma^2);
Y_freq = G_freq.*X_freq;
OUTPUT = uint8(real(IFFT2D(Y_freq))); % taking the inverse image
end
%% Wiener filtered image function
function OUTPUT = m_wiener(INPUT)
[M,N] = size(INPUT);
IMG_FFT = FFT2D(INPUT); % change INPUT to freqency domain
sigma = 100;
H = LPF_add(M,N); % apply the LPF to get matrix in freq
PY = abs(IMG_FFT).^2/M^2;
SH = H.*(abs(H)>0)+(abs(H)==0);
IVH = 1./SH;
IVH = IVH.*(abs(H)>1)+abs(SH).*IVH.*(abs(SH)<=1);
PY = PY.*(PY>sigma^2)+sigma^2*(PY<=sigma^2);
temp = IVH.*(PY-sigma^2)./PY;
pred = temp.*IMG_FFT; % multiply in freq
OUTPUT = real(IFFT2D(pred));% Restorated image without denoising
end
%% Constrained least square image function
function OUTPUT = m_least_sqr(INPUT, lambda)
[M,N] = size(INPUT); % image sizes
tmp = (1:M)'*ones(1,N)+ones(M,1)*(1:N);
tmp = 1-2*rem(tmp,2);
Hh = LPF_add(M,N); % Using LPF filter to blur image
H2 = abs(Hh).^2;
Gg = FFT2D(INPUT); % Perform FFT2D of blurred+noise image
pa = [0 -1 0; -1 4 -1; 0 -1 0]; % Constrained least square method using smoothing operator
pad_p = zeros(M,N);
pad_p(1:size(pa,1),1:size(pa,2)) = pa.*tmp(1:size(pa,1),1:size(pa,2)); % pads matrix X with zeros to size MROWS-by-NCOLS before transforming
P = FFT2D(pad_p);
K = conj(Hh).*Gg./(H2+lambda*abs(P).^2 + eps);
OUTPUT = abs(IFFT2D(K));
end

%% Code for FFT2D
function y = FFT2D(x)
y = myb(myb(x).').'; % perform FFT for the vector that is column of the matrix x
end
function y = IFFT2D(x)
[M,N] = size(x);
y = conj(FFT2D(conj(x)))/(M*N); % using the formulation of inverse FFT in slide
end
function y = myb(x)
y = zeros(size(x,1),1);
for i = 1:size(x,2)
    y = [y myFFT(x(:,i))];
end
y = y(:,2:end);
end
function y = myFFT(x)
%function to get exactly the output is row or column
y = mya(x);
if (isrow(x)~= 0)
    y = y';
end
end
function y = mya(x)
x = x(:);
N = length(x);
even = x(1:2:end); % get the even part
odd = x(2:2:end); % get the odd part
if (N >= 8) %perform the incurve algorithm
    x_even = mya(even);
    x_odd = mya(odd);
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
%% Additional functions
function y = add_gauss_noise (x)
SNR_dB = 6.2;
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
function y = LPF_add(M,N)
% M,N is size of filter
F = 60; %The cut off frequency initialization
a = 0:(M-1);
b = 0:(N-1);
index = find(a>M/2);
a(index) = a(index)-M;
indey = find(b>N/2);
b(indey) = b(indey)-N;
[B,A] = meshgrid(b,a);
T = sqrt(A.^2+B.^2);
y = double(T<=F); % Cut-off frequency comparison, return the freq LPF matrix
end

