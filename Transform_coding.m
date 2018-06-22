
%% How to run: Just click "Run" button and hit "Select folder"

function Transform_coding()
close all;
clc;
imgdir = uigetdir('Test_images');
%% Loading directory, gray image
file = fopen(fullfile(imgdir,'\cameraman_gray_256x256.raw'),'rb');
gray_image = fread(file,fliplr([256,256]),'*uint8')';
fclose(file); gray_image = double(gray_image);
%% Starting
disp('Transform coding (8x8 block DCT + quantization)');
disp('Computing the 8x8 block DCT image...');
disp('Quantizing image...');
disp('Dequantizing image...');
disp('Inversing DCT 8x8 block image...');

% Computes output img
DCT8_img = DCT8(gray_image);
quantization_img = quantization(DCT8_img);
dequantization_img = dequantization(quantization_img);
IDCT8_img = IDCT8(dequantization_img);

% Computes psnr
psnr_inversed_img = psnr(uint8(IDCT8_img),uint8(gray_image));
disp(['PSNR: ', num2str(psnr_inversed_img), ' dB']);

% Displays figures
figure; imshow(gray_image,[]); title('Origional'); % show original image
figure; imshow(DCT8_img,[]); title('DCT block 8x8 image');
figure; imshow(quantization_img,[]); title('DCT block 8x8 quantization image');
figure; imshow(dequantization_img,[]); title('DCT block 8x8 dequantization image'); % show dequantization image
figure; imshow(IDCT8_img,[]); title('DCT block 8x8 inversed image'); % show inversed image

disp('Computing Completed!');
end
%% function to calculate DCT8x8
function y = DCT8(x)
    [M,N] = size(x);
    y = zeros(M,N);
    for i=1:8:M
        for j=1:8:N %jump to each block 8x8
            block = myDCT2D(x(i:i+7,j:j+7));
            y(i:i+7,j:j+7) = block; % give the output the value when apply DCT for each block
        end
    end
end
%% function quantization 8x8 block
function y = quantization(x)
% input must be 8x8 block DCT
[M,N] = size(x);
q_matrix = ... % the quantization parameter from 1 to 31
  [	31	22	23	21	24	5	22	4	;
	19	25	5	20	1	4	22	13	;
	26	10	30	20	23	26	2	11	;
	13	19	30	9	4	20	30	14	;
	31	2	29	14	22	22	8	11	;
	12	10	11	1	14	29	25	25	;
	2	21	24	13	3	25	21	7	;
    25	3	8	12	24	22	19	9	];
% you can change an abitrary matrix for quantization
y = zeros(M,N);
for i=1:8:M
    for j=1:8:N %jump to each block 8x8
        block = x(i:i+7,j:j+7); % block 8x8
        y(i:i+7,j:j+7) = round( block./(2*q_matrix)); % quantizing by dividing image block 8x8 to q_matrix and round them
    end
end
end
%% function dequantization 8x8 block image
function y = dequantization(x)
% input must be 8x8 block DCT
[M,N] = size(x);
q_matrix = ...
  [	31	22	23	21	24	5	22	4	;
	19	25	5	20	1	4	22	13	;
	26	10	30	20	23	26	2	11	;
	13	19	30	9	4	20	30	14	;
	31	2	29	14	22	22	8	11	;
	12	10	11	1	14	29	25	25	;
	2	21	24	13	3	25	21	7	;
    25	3	8	12	24	22	19	9	];
% this q-matrix must be same as the forward transform
y = zeros(M,N);
for i=1:8:M
    for j=1:8:N %jump to each block 8x8
        block = x(i:i+7,j:j+7); % block 8x8
        y(i:i+7,j:j+7) = block.*(2*q_matrix); % dequantizing by multiplying image block 8x8 to q_matrix
    end
end
end
%% function to calculate IDCT8x8
function y = IDCT8(x)
    [M,N] = size(x);
    y = zeros(M,N);
    for i=1:8:M
        for j=1:8:N %jump to each block 8x8
            block = myIDCT2D(x(i:i+7,j:j+7));
            y(i:i+7,j:j+7) = block; % give the output the value when apply IDCT for each block
        end
    end
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