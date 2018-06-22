
    %% How to run: Just click "Run" button and hit "Select folder"
    
function Histogram_Filter_Salt_Pepper_Gaussian_Noise()
    clc;
    clear;
    close all;
    %% read color image
    imgdir = uigetdir('Test_images');
    file = fopen(fullfile(imgdir,'\Color_monarch_512x512.raw'),'rb');
    color_image = fread(file,fliplr([512,512*3]),'*uint8')';
    fclose(file);
    %% get the r,g,b channels
    r = color_image(:,1:3:end); %R matrix
    g = color_image(:,2:3:end); %G matrix
    b = color_image(:,3:3:end); %B matrix
    color_image = cat(3, r, g, b); %RGB image
    %% read gray noise image
    file = fopen(fullfile(imgdir,'\Gray_monarch_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);
    %% gray histogram
    disp('Calculating for Point operations...');
    figure; imshow(gray_image,[]); title('Gray Origional');
    gray_image_1 = double(gray_image);
    gray_histo = histEqualization(gray_image_1);
    figure; imshow(gray_histo,[]); title('Gray histogram');
    %% plot histogram gray origional
    Hist1 = histogram_gray(gray_image_1);
    figure; bar(Hist1(:,1),Hist1(:,2));
    grid on;
    ylabel('No. of pixels with intensity levels ---->');
    xlabel('Intensity levels (gray) ---->'); title('HISTOGRAM OF THE gray Origional IMAGE')
    %% plot histogram gray equalization
    Hist2 = histogram_gray(gray_histo);
    figure; bar(Hist2(:,1),Hist2(:,2));
    grid on;
    ylabel('no. of pixels with intensity levels---->');
    xlabel('intensity levels---->'); title('HISTOGRAM OF THE gray equalization IMAGE')
    %% RGB histogram
    figure; imshow(color_image,[]); title('RGB Origional');
    r_histo = histEqualization(r);
    g_histo = histEqualization(g);
    b_histo = histEqualization(b);
    rgb_histo = cat(3, r_histo, g_histo, b_histo);
    figure; imshow(rgb_histo,[]); title('RGB histogram');
    %% get HSI
    R = im2double(r);
    G = im2double(g);
    B = im2double(b);
    theta = acos((0.5*((R-G)+(R-B)))./((sqrt((R-G).^2+(R-B).*(G-B)))+eps));
    H = theta;
    H(B>G) = 2*pi-H(B>G);
    H = H/(2*pi);
    S = 1-3.*(min(min(R,G),B))./(R+G+B+eps);
    I = (R+G+B)/3;
    HSI = cat(3,H,S,I); %HSI image
    H = uint8(255*H); S = uint8(255*S); I = uint8(255*I);
    % Histogram of HSI
        H_histo = histEqualization(H);
        S_histo = histEqualization(S);
        I_histo = histEqualization(I);
    HSI_histo = im2double(cat(3, H_histo, S_histo, I_histo));
    % display histogram HSI
    figure; imshow(HSI_histo,[]); title('HSI histogram');
    % display image HSI
    figure; imshow(HSI, []); title('HSI origional');
    %% get CMY
    C = 1 - R;
    M = 1 - G;
    Y = 1 - B;
    CMY = (cat(3, C, M, Y));
    C=uint8(255*C);M=uint8(255*M);Y=uint8(255*Y);
    C_histo = histEqualization(C);
    M_histo = histEqualization(M);
    Y_histo = histEqualization(Y);
    %combine C,M,Y to CMY image
    CMY_histo = cat(3, C_histo, M_histo, Y_histo);
    CMY_histo = im2double(CMY_histo);
    % DISPLAY CMY
    figure; imshow(CMY,[]); title('CMY ORIGIONAL');
    figure; imshow(CMY_histo,[]); title('CMY histogram');
    %% 2.1b Histogram modification
    alpha = 0.005;
    modification_image = histogramModification(gray_image_1,alpha);
    figure; imshow(modification_image,[]); title(['Histogram modification with alpha = ',num2str(alpha)]);
    %% 2.2a) Gamma correction
    gamma_init = 5.0;
    gamma_image = gammacorrection(gray_image_1,gamma_init);
    figure; imshow(gamma_image,[]); title(['Gamma correction = ',num2str(gamma_init)]);

    %% 2.2b) Spatial filtering
    tic; disp('Calculating for Spatial filtering...');
    %% 2.2b)- Salt & pepper noise, kernel 3, mean filtering
    SaltPepper_image_3 = AddSaltPepper(gray_image_1);
    meanFiltering_image_3 = meanFiltering(SaltPepper_image_3,3);
    figure; imshow(SaltPepper_image_3,[]); title('Gray SaltPepper noise');
    figure; imshow(meanFiltering_image_3,[]); title('Gray SaltPepper meanFiltering (kernel 3)');

    SP_3_mean_filter = psnr(uint8(meanFiltering_image_3),uint8(gray_image_1))
    
    %% 2.2b)- Salt & pepper noise, kernel 7, mean filtering
    SaltPepper_image_7 = AddSaltPepper(gray_image_1);
    meanFiltering_image_7 = meanFiltering(SaltPepper_image_7,7);
    figure; imshow(meanFiltering_image_7,[]); title('Gray SaltPepper meanFiltering (kernel 7)');

    SP_7_mean_filter = psnr(uint8(meanFiltering_image_7),uint8(gray_image_1))
    
    %% 2.2b)- Gaussian noise, kernel 3, mean filtering
    Gaussian_image_3 = AddSaltPepper(gray_image_1);
    meanFiltering_image_3 = meanFiltering(Gaussian_image_3,3);
    figure; imshow(Gaussian_image_3,[]); title('Gray Gaussian noise');
    figure; imshow(meanFiltering_image_3,[]); title('Gray Gaussian meanFiltering (kernel 3)');
 
    Gau_3_mean_filter = psnr(uint8(meanFiltering_image_3),uint8(gray_image_1))

    %% 2.2b)- Gaussian noise, kernel 7, mean filtering
    Gaussian_image_7 = AddSaltPepper(gray_image_1);
    meanFiltering_image_7 = meanFiltering(Gaussian_image_7,7);
    figure; imshow(meanFiltering_image_7,[]); title('Problem 2.2b) Gray Gaussian meanFiltering (kernel 7)');

    Gau_7_mean_filter = psnr(uint8(meanFiltering_image_7),uint8(gray_image_1))
    
    %% 2.2b)- Salt & pepper noise, kernel 3, median filtering
    SaltPepper_image_3 = AddSaltPepper(gray_image_1);
    medianFiltering_image_3 = medianFiltering(SaltPepper_image_3,3);
    figure; imshow(medianFiltering_image_3,[]); title('Gray SaltPepper medianFiltering (kernel 3)');

    SP_3_med_filter = psnr(uint8(medianFiltering_image_3),uint8(gray_image_1))
    
    %% 2.2b)- Salt & pepper noise, kernel 7, median filtering
    SaltPepper_image_7 = AddSaltPepper(gray_image_1);
    medianFiltering_image_7 = medianFiltering(SaltPepper_image_7,7);
    figure; imshow(medianFiltering_image_7,[]); title('Gray SaltPepper medianFiltering (kernel 7)');

    SP_7_med_filter = psnr(uint8(medianFiltering_image_7),uint8(gray_image_1))
    
    %% 2.2b)- Gaussian noise, kernel 3, median filtering
    Gaussian_image_3 = AddGaussian(gray_image_1);
    medianFiltering_image_3 = medianFiltering(Gaussian_image_3,3);
    figure; imshow(Gaussian_image_3,[]); title('Gray Gaussian noise');
    figure; imshow(medianFiltering_image_3,[]); title('Gray Gaussian medianFiltering (kernel 3)');
    
    Gau_3_med_filter = psnr(uint8(medianFiltering_image_3),uint8(gray_image_1))
    
    %% 2.2b)- Gaussian noise, kernel 7, median filtering
    Gaussian_image_7 = AddGaussian(gray_image_1);
    medianFiltering_image_7 = medianFiltering(Gaussian_image_7,7);
    figure; imshow(medianFiltering_image_7,[]); title('Gray Gaussian medianFiltering (kernel 7)');
    
    Gau_7_med_filter = psnr(uint8(medianFiltering_image_7),uint8(gray_image_1))
    
    %% 2.2b)- Salt & pepper noise, kernel 3, smooth filtering
    SaltPepper_image_3 = AddSaltPepper(gray_image_1);
    smoothFiltering_image_3 = smoothFiltering(SaltPepper_image_3,3);
    figure; imshow(smoothFiltering_image_3,[]); title('Gray SaltPepper smoothFiltering (kernel 3)');
    
    SP_3_dir_filter = psnr(uint8(smoothFiltering_image_3),uint8(gray_image_1))
    
    %% 2.2b)- Salt & pepper noise, kernel 7, smooth filtering
    SaltPepper_image_7 = AddSaltPepper(gray_image_1);
    smoothFiltering_image_7 = smoothFiltering(SaltPepper_image_7,7);
    figure; imshow(smoothFiltering_image_7,[]); title('Gray SaltPepper smoothFiltering (kernel 7)');
    
    SP_7_dir_filter = psnr(uint8(smoothFiltering_image_7),uint8(gray_image_1))
    
    %% 2.2b)- Gaussian noise, kernel 3, smooth filtering
    Gaussian_image_3 = AddGaussian(gray_image_1);
    smoothFiltering_image_3 = smoothFiltering(Gaussian_image_3,3);
    figure; imshow(smoothFiltering_image_3,[]); title('Gray Gaussian smoothFiltering (kernel 3)');    
    
    Gau_3_dir_filter = psnr(uint8(smoothFiltering_image_3),uint8(gray_image_1))
    
    %% 2.2b)- Gaussian noise, kernel 7, smooth filtering
    Gaussian_image_7 = AddGaussian(gray_image_1);
    smoothFiltering_image_7 = smoothFiltering(Gaussian_image_7,7);
    figure; imshow(smoothFiltering_image_7,[]); title('Gray Gaussian smoothFiltering (kernel 7)');  

    Gau_7_dir_filter = psnr(uint8(smoothFiltering_image_7),uint8(gray_image_1))
    
    toc;
end

%% Getting the Histogram of a particular image
function  output  = histEqualization(input)
    % input: original grayscale image (uint8)
    % output: Adjusted grayscale image after applying histogram equalization(uint8)
    Level = 256;     %number of levels
    counts = zeros(256,1);
    for i=0:255
        counts(i+1)=sum(sum(input==i));
    end
    %compute PDF of the image
    proba = counts/(size(input,1)*size(input,2));
    %compute s=(L-1)*CDF(p)
    s = (Level-1)*cumsum(proba);
    %Perform the round the float into the nearest integare
    s = round(s);
    %Perform Map the value of s to corresponding pixels in the image
    output = uint8(zeros(size(input)));
    for k = 1:size(s,1)
        output(input==k-1) = s(k);
    end
end

%% Get the histogram of a particular image without using inbuilt function
function output = histogram_gray(x)
    [M,N]=size(x);
    L = 1:256;
    n = 0:255;
    count=0;
    for z=1:256
        for i=1:M
            for j=1:N
                if x(i,j)==z-1
                    count = count+1; %calculate the number of pixels
                end
            end
        end
                L(z)=count;
                count=0;
    end
    output = [n' L']; 
end
%% 2.2a) Gamma correction
function Correction = gammacorrection(x,Gamma)
% Function to calculate the Gamma Correction for a particular image
if (Gamma<0)
    Gamma = 1;
end
x = double(x);
Correction = 255 * (x/255).^ Gamma;
end
%% Add salt & pepper noise
function y = AddSaltPepper (x)
    %customize black and white to adjust noise
    black = 10;
    white = 250;
    noiseImg = x;
    randMatrix = randi(255,size(x,1),size(x,2));
    noiseImg(randMatrix <= black) = 0;
    noiseImg(randMatrix >= white) = 255;
    y = noiseImg;
end
%% Add Gaussian noise
function y = AddGaussian (x)
    %y=awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
    %resulting signal vector y of specified SNR in dB
    SNR_dB = 10;
    SNR = 10^(SNR_dB/10); %SNR to linear scale
    [M,L] = size(x);
    rng('default');%set the random generator seed to default (for comparison only)
    y = zeros(M,L);
    for i=1:M
        Esym = sum(abs(x(i,:)).^2)/(L); %Calculate actual symbol energy
        N0 = Esym/SNR; %Find the noise spectral density
        noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
        n = noiseSigma*randn(1,L);%computed noise
        y(i,:) = x(i,:) + n; %received signal
    end
end

%% Median filtering
function y = medianFiltering(x,windowSize)
    %windowSize = 3 or 7 (3x3 or 7x7)
    [m,n] = size(x);
    M = windowSize; N = windowSize;
    %PAD THE MATRIX WITH ZEROS ON ALL SIDES
    %x_pad = padarray(x,[floor(M/2),floor(N/2)]);
    x_pad = zeros(m+M-1,n+N-1);
    x_pad(floor(M/2)+1:1:floor(M/2)+m,floor(N/2)+1:1:floor(N/2)+n) = x;
    y = zeros([size(x,1) size(x,2)]);
    med_indx = round((M*N)/2); %MEDIAN INDEX
    for i=1:size(x_pad,1)-(M-1)
        for j=1:size(x_pad,2)-(N-1)
            temp = x_pad(i:i+(M-1),j:j+(N-1),:);
            tmp_sort = sort(temp(:)); %tmp(:) converts 2D matrix to 1D matrix
            y(i,j) = tmp_sort(med_indx); 
        end
    end
end

%% Mean filtering 
function y = meanFiltering(x,windowSize)
    [m,n] = size(x);
    M = windowSize; N = windowSize;
    x = double(x);
    sz = m*n;
    %PAD THE MATRIX WITH ZEROS ON ALL SIDES
    x_pad = zeros(m+M-1,n+N-1);
    x_pad(floor(M/2)+1:1:floor(M/2)+m,floor(N/2)+1:1:floor(N/2)+n) = x;
    C = x_pad;
    %Define the mean matrix initiation
    lvar = zeros([size(x,1) size(x,2)]);
    lmean = zeros([size(x,1) size(x,2)]);
    %temp = zeros([size(x,1) size(x,2)]);
    %NewImg = zeros([size(x,1) size(x,2)]);
    
    for i = 1:size(C,1)-(M-1)
        for j = 1:size(C,2)-(N-1)
            temp = C(i:i+(M-1),j:j+(N-1));
            tmp =  temp(:);
            %Find the local mean and local variance for the local region        
            lmean(i,j) = mean(tmp);
            lvar(i,j) = mean(tmp.^2)-mean(tmp).^2;
        end
    end
    %Noise variance and average of the local variance
    nvar = sum(lvar(:))/sz;
    %If noise_variance > local_variance then local_variance=noise_variance
    lvar = max(lvar,nvar);     
    %Final_Image = x - (noise variance/local variance)*(x-local_mean);
    NewImg = nvar./lvar;
    NewImg = NewImg.*(x-lmean);
    NewImg = x-NewImg;
    %Convert the image to uint8 format.
    y = uint8(NewImg);
end

%% Smooth directional

function output = smoothFiltering(x,windowSize)	
    %[m,n] = size(x);
    M = windowSize; N = windowSize;
    A = x;
    wind = ones(M,N)/(M*N);
    %padding zero for input
    p1 = size(wind,1)-1;
    p2 = size(wind,2)-1;
    y = uint8(zeros(size(A)));
    if(size(wind,1)==1) % check window rows to 1
     B = zeros(size(A,1),size(A,2)+p2);
     m = 0;
     n = floor(size(wind,2)/2);
     sz1 = size(B,1);
     sz2 = size(B,2)-p2;
    elseif(size(wind,2)==1) % check window columns to 1
        B = zeros(size(A,1)+p1,size(A,2));
        m = floor(size(wind,1)/2);
        n = 0;
        sz1 = size(B,1)-p1; % minus padding size
       sz2 = size(B,2);
    else
        B = zeros(size(A,1)+p1,size(A,2)+p2);
        m = floor(size(wind,1)/2);
        n = floor(size(wind,2)/2);
        sz1 = size(B,1)-p1; % minus padding size
     sz2 = size(B,2)-p2;
    end
     for i=1:size(A,1)
         for j=1:size(A,2)
             B(i+m,j+n) = A(i,j);
         end
     end
    szcorr1 = size(wind,1);
    szcorr2 = size(wind,2);
    for i=1:sz1
        for j=1:sz2
            sum=0;
            m=i;
            n=j;
            for x=1:szcorr1
              for y=1:szcorr2
           %Performed the weighted sum of the neighborhood pixels
                   sum = sum+(B(m,n)*wind(x,y));
                   n = n+1;                    
              end
                 n=j;
                m=m+1;
            end
            output(i,j)= sum;
        end
    end
end
%% Histogram modification
%% histogram modification(with exponential distribution)
function out = histogramModification(x,alpha)
    %alpha = 0.01;
    factor = histogram_factor(x);
    P = distributionEXP(alpha);
    for i = 1:256
        val(i) = sum(P(1:i));
    end
    val = val./val(256);

    for i = 1:256
        temp = factor(i) <= val;
        index = find(temp == 1);
        n(i) = index(1) - 1;
    end
    for i = 1:512
        for j = 1:512
            out(i,j) = n(x(i,j)+1);
        end
    end
end
% function for return the distribution of exponential
function P = distributionEXP(alpha)
    vmin = 0;
    for i = 1:256
        v = i;
        P(i) = alpha*exp(-alpha*(v-vmin));
    end
end
% function for return the coefficients of histogram
function out = histogram_factor(x)
    for i = 1:256
        number(i) = sum(sum(x==i-1));
    end
    count = number/(512*512);
    for i = 1:256
        out(i) = sum(count(1:i));
    end
end
