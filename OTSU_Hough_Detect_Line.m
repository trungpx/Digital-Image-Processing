
%% How to run: Just click "Run" button and hit "Select folder"

function OTSU_Hough_Detect_Line()
clc;
clear;
close all;
%% Loading directory and read gray images
imgdir = uigetdir('Test_images');
% read gray image1
file = fopen(fullfile(imgdir,'texture1_gray_256x256.raw'),'rb');
gray_image1 = fread(file,fliplr([256,256]),'*uint8')';
fclose(file);
% read gray image2
file = fopen(fullfile(imgdir,'texture2_gray_256x256.raw'),'rb');
gray_image2 = fread(file,fliplr([256,256]),'*uint8')';
fclose(file);
% read gray image3
file = fopen(fullfile(imgdir,'hough_gray_256x256.raw'),'rb');
gray_image4 = fread(file,fliplr([256,256]),'*uint8')';
fclose(file);
% read gray image4
file = fopen(fullfile(imgdir,'object_gray_360x285.raw'),'rb');
gray_image3 = fread(file,fliplr([285,360]),'*uint8')';
fclose(file);
%% Starting
gray_image1 = double(gray_image1);
gray_image2 = double(gray_image2);
gray_image3 = double(gray_image3);
gray_image4 = double(gray_image4);

disp('Image analysis');
disp('Computing...');

% Computes the output images
[slit_45, slit_135] = slit_image(256);
% Compute FFT 2.1
fft_image1 = FFT2D(gray_image1);
fft_image2 = FFT2D(gray_image2);
% Multiply slit 2.1
fft_image1_45 = fft_image1.*slit_45;
fft_image1_135 = fft_image1.*slit_135;
fft_image2_45 = fft_image2.*slit_45;
fft_image2_135 = fft_image2.*slit_135;

% The ratios of the energy in the slit 2.1
t1_45 = sum(sum(abs(fft_image1_45)))/sum(sum(abs(fft_image1)));
t1_135 = sum(sum(abs(fft_image1_135)))/sum(sum(abs(fft_image1)));
t2_45 = sum(sum(abs(fft_image2_45)))/sum(sum(abs(fft_image2)));
t2_135 = sum(sum(abs(fft_image2_135)))/sum(sum(abs(fft_image2)));

% Computes IFFT 2.1
fft_image1_45_res = IFFT2D(fft_image1_45);
fft_image1_135_res = IFFT2D(fft_image1_135);
fft_image2_45_res = IFFT2D(fft_image2_45);
fft_image2_135_res = IFFT2D(fft_image2_135);

% Computes 2.2 and 2.3 parts
out_threshold_img = extract_otsu(gray_image3);
[out_hough,~,~] = hough_trans(gray_image4);

% Displays the ratios of the energy in the slit
disp(['t1_45 = ', num2str(t1_45)]);
disp(['t1_135 = ',num2str(t1_135)]);
disp(['t2_45 = ', num2str(t2_45)]);
disp(['t2_135 = ', num2str(t2_135)]);

% Displaying figures
figure; imshow(slit_45,[]); title('Image-slit-45');
figure; imshow(slit_135,[]); title('Image-slit-135');
figure; imshow(gray_image1,[]); title('Original image1');
figure; imshow(log(1+abs(fft_image1)),[]); title('fft-image1');
figure; imshow(gray_image2,[]); title('Original image2');
figure; imshow(log(1+abs(fft_image2)),[]); title('fft-image2');
figure; imshow(log(1+abs(fft_image1_45)),[]); title('fft-image1-45');
figure; imshow(log(1+abs(fft_image1_135)),[]); title('fft-image1-135');
figure; imshow(log(1+abs(fft_image2_45)),[]); title('fft-image2-45');
figure; imshow(log(1+abs(fft_image2_135)),[]); title('fft-image2-135');
figure; imshow(real(fft_image1_45_res),[]); title('image1 recovered from slit45');
figure; imshow(real(fft_image1_135_res),[]); title('image1 recovered from slit135');
figure; imshow(real(fft_image2_45_res),[]); title('image2 recovered from slit45');
figure; imshow(real(fft_image2_135_res),[]); title('image2 recovered from slit135');
figure; imshow(gray_image3,[]); title('object\_gray\_360x285.raw'); disp('Processing...');
figure; imshow(out_threshold_img,[]); title('Output OTSU object '); disp('Processing...');
figure; imshow(gray_image4,[]); title('Input hough\_gray\_256x256.raw');
figure; imshow(rescale_img(out_hough,[256 256]),[]); title('Parameter space image(Hough)');
draw_detected_lines(gray_image4);

disp('Completed...');
%% Finish
end

%% Define required functions
function [OUTPUT, theta_out, rho_out] = hough_trans(INPUT)
% Peforms the Hough transform. Compute Hough accumulator array for finding lines.
[rows, cols] = size(INPUT);
rho_maximum = floor(sqrt(rows^2 + cols^2)) - 1;
theta_max_degree = 90;
theta_out = -theta_max_degree:theta_max_degree - 1;
rho_out = -rho_maximum:rho_maximum;
OUTPUT = zeros(length(rho_out), length(theta_out));
% Starting Hough algorithm
for row = 1:rows
    for col = 1:cols
        if INPUT(row, col) > 0
            x = col - 1;
            y = row - 1;
            for theta = theta_out
                rho_out = round((x * cosd(theta)) + (y * sind(theta)));
                rho_index = rho_out + rho_maximum + 1;
                theta_index = theta + theta_max_degree + 1;
                OUTPUT(rho_index, theta_index) = OUTPUT(rho_index, theta_index) + 1;
            end
        end
    end
end
end
function OUTPUT = hough_peaks(H, varargin)
% Find peaks in a Hough accumulator array
temp = inputParser;
addOptional(temp, 'numpeaks', 1, @isnumeric);
addParameter(temp, 'Threshold', 0.5 * max(H(:))); % add some options
addParameter(temp, 'NHoodSize', floor(size(H) / 100.0) * 2 + 1);  % Make sure odd values above size(H)/50
parse(temp, varargin{:});
countpeaks = temp.Results.numpeaks;
thres = temp.Results.Threshold;
nHoodSize = temp.Results.NHoodSize;
OUTPUT = zeros(countpeaks, 2);
count = 0;
while(countpeaks>count)
    maxH = max(H(:));
    if (maxH >= thres)
        count = count + 1; [a,b] = find(H == maxH);
        OUTPUT(count,:) = [a(1),b(1)];
        aStart = max(1, a - (nHoodSize(1) - 1) / 2);
        aEnd = min(size(H,1), a + (nHoodSize(1) - 1) / 2);
        bStart = max(1, b - (nHoodSize(2) - 1) / 2);
        bEnd = min(size(H,2), b + (nHoodSize(2) - 1) / 2);
        for i = aStart : aEnd
            for j = bStart : bEnd
                H(i,j) = 0;
            end
        end
    else
        break;
    end
end
OUTPUT = OUTPUT(1:count, :);
end
function OUTPUT = extract_otsu(INPUT)
% Change histogram
histnum = zeros(256,1);
for i=0:255
    histnum(i+1) = sum(INPUT(:)==i);%length(find(gray_image=i))
end
all = sum(histnum); % Number of pixels in the given image
%% OTSU starting algorithm
temp = 0; w_B = 0; M = 0.0;
dot_v = dot((0:255), histnum);
for i=1:256
    w_B = w_B + histnum(i);
    wF = all - w_B;
    if (w_B == 0 || wF == 0)
        continue;
    end
    temp = temp +  (i-1) * histnum(i); % calculate temporary value
    mF = (dot_v - temp)/wF;
    middle = w_B*wF*((temp/w_B)-mF)*((temp/w_B)-mF);
    if ( middle >= M-1)
        OUTPUT = i;
        M = middle;
    end
end
OUTPUT = INPUT > OUTPUT-2; % determine level 0 and 1
end
function draw_detected_lines(INPUT)
[gx,gy] = gradient(INPUT);
% Angles range from -pi/2 to +pi/2
Gangle = atan(gy./gx);
bound = find_edge(INPUT); % define new function to find edge
%bound = edge(INPUT, 'canny'); Use ML func easily
[rows,cols] = size(bound);
distMax = sqrt(rows^2 + cols^2); % Max possible distance from origin
theta = -90:1:89; % range of theta values
rho = -distMax:1:distMax; % range of rho values
H_img = zeros(length(rho),length(theta)); % Allocate accumulator array
% Scan through edge image
for i1=1:cols
    for i2=1:rows
        if bound(i2,i1) ~= 0
            t_i = Gangle(i2,i1);
            % Calculate distance from origin, given this angle
            dist = i1*cos(t_i) + i2*sin(t_i);
            % Find theta value that is closest to this
            [~, iTheta] = min(abs(theta-t_i*180/pi));
            % Find rho value that is closest to this
            [d, rho_i] = min(abs(rho-dist));
            if d <= 1
                H_img(rho_i,iTheta) = H_img(rho_i,iTheta) + 1; % Inc accumulator
            end
        end
    end
end

peaks_out = hough_peaks(H_img,50,"threshold",ceil(0.08*max(H_img(:)))); % Find all the peaks
% Draw lines. Equation of a line: rho = x*cos(theta) + y*sin(theta)
figure; imshow(INPUT, []); title('Detected straight lines');
for i=1:size(peaks_out,1)
    r_i = rho(peaks_out(i,1));
    t_i = theta(peaks_out(i,2));
    t_i = t_i * pi/180; % convert to radians
    if t_i == 0
        m0 = r_i; m1 = r_i;
        n0 = 1; n1 = size(bound,2);
    else
        m0 = 1;
        n0 = (r_i - m0*cos(t_i)) / sin(t_i);
        m1 = size(bound,2);
        n1 = (r_i - m1*cos(t_i)) / sin(t_i);
    end
    line([m0 m1], [n0 n1], 'Color', 'c'); % drawing the dettected lines
end
end
function OUTPUT = rescale_img(INPUT, out_dims)
% %resize image as bilinear Interpolation. Get some necessary variables first
out_rows = out_dims(1); out_cols = out_dims(2);
[in_rows,in_cols] = size(INPUT);
R_scale = in_rows/out_rows; C_scale = in_cols/out_cols;
% Define grid of co-ordinates in our image. Generate (x,y) pairs for each point in our image
[CF, RF] = meshgrid(1:out_cols, 1:out_rows);
RF = RF * R_scale; CF = CF * C_scale; R = floor(RF); C = floor(CF);
R(R < 1) = 1; C(C < 1) = 1; R(R > in_rows - 1) = in_rows - 1;
C(C > in_cols - 1) = in_cols - 1;
delta_R = RF - R; delta_C = CF - C;
% Final line of algorithm. Get column major indices for each point we wish to access
in1 = sub2ind([in_rows, in_cols], R, C);
in2 = sub2ind([in_rows, in_cols], R+1,C);
in3 = sub2ind([in_rows, in_cols], R, C+1);
in4 = sub2ind([in_rows, in_cols], R+1, C+1);
% Now interpolate; Go through each channel for the case of colour; Create output image that is the same class as input
OUTPUT = zeros(out_rows, out_cols, size(INPUT, 3));
OUTPUT = cast(OUTPUT, class(INPUT));
for dx = 1 : size(INPUT, 3)
    channel = double(INPUT(:,:,dx)); % Get i'th channel
    % Interpolate the channel
    temp = channel(in1).*(1 - delta_R).*(1 - delta_C) + ...
        channel(in2).*(delta_R).*(1 - delta_C) + ...
        channel(in3).*(1 - delta_R).*(delta_C) + ...
        channel(in4).*(delta_R).*(delta_C);
    OUTPUT(:,:,dx) = cast(temp, class(INPUT));
end
end
function [OUT45,OUT135] = slit_image(n)
% Make filter
OUT45 = zeros(n);
OUT135 = zeros(n);
for i = 1:n
    OUT135(i,i) = 1;
    OUT45(i,n+1-i) = 1;
end
end
function OUTPUT = find_edge(INPUT)
C = INPUT;
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Sobel mask for x-direction:
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Sobel mask for y-direction:
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        OUTPUT(i,j)=sqrt(Gx.^2+Gy.^2);
    end
end

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
