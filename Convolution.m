
function Convolution()
clear variables;
 %% Loading directory for image files
 imgdir = uigetdir('Test_images');
 file = fopen(fullfile(imgdir,'Gray_monarch_noisy_512x512.raw' ), 'rb' );
 gray_image = fread(file,fliplr([512,512]), '*uint8' )'; % input matrix for image Gray_monarch_noisy_512x512.raw
 fclose(file);
 %%
H1 = [1 2 1; 2 4 2; 1 2 1]/16; % kernel H1
H2 = [1 4 1; 4 -20 4; 1 4 1]/6; % kernel H2
output_image_kelnel1 = InnerFunction(gray_image, H1); % Output after convolution with kernel 1
output_image_kelnel2 = InnerFunction(gray_image, H2); % Output after convolution with kernel 2
 %% Displaying figures (Edit this part as needed)
 figure('Name','DIP', 'NumberTitle','off'); imshow(output_image_kelnel1,[]); title('Conv of Gray\_monarch\_noisy (kernel H1)' ); % display the image after convolution H1
 figure('Name','DIP', 'NumberTitle','off'); imshow(output_image_kelnel2); title('Conv of Gray\_monarch\_noisy (kernel H2)' ); % display the image after convolution H2
 %%---------------------------------------------------------------
end
function OUTPUT = InnerFunction(INPUT, kernel)
%%------------------------------------------------------------------
%PAD THE MATRIX WITH ZEROS
[a,b] = size(INPUT);
 B = INPUT;
 B = [zeros(1,b); B];
 B = [B; zeros(1,b)];
 B = [B zeros(a+2,1)];
 B = [zeros(a+2,1) B];
 
%B = padarray(INPUT, [1 1]); Don't use
% PRE-ALLOCATE THE MATRIX
OUTPUT = zeros([size(INPUT,1) size(INPUT,2)]);
%PERFORM COONVOLUTION
for i = 1:size(B,1)-2
    for j = 1:size(B,2)-2
        current = double(B(i:i+2,j:j+2)).*kernel;
        OUTPUT(i,j) = sum(current(:));
    end
end
%%------------------------------------------------------------------
end

