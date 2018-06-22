
function RGB_HSI_YCbCr()
clc;
clear;
close all;
 %% Loading directory for image files
 %RGB image components
 imgdir = uigetdir('Test_images' );
 file = fopen(fullfile(imgdir,'\Color_monarch_512x512.raw' ), 'r' );
 color_image = fread(file,fliplr([512,512*3]), '*uint8' )'; % input matrix for image Color_monarch_512x512.raw
 fclose(file);
 r = color_image(:,1:3:end); %R matrix
 g = color_image(:,2:3:end); %G matrix
 b = color_image(:,3:3:end); %B matrix
 color_image = cat(3, r, g, b); %RGB image
 %HSI transformation components
 R = im2double(r);
 G = im2double(g);
 B = im2double(b);
 theta=acos((0.5*((R-G)+(R-B)))./((sqrt((R-G).^2+(R-B).*(G-B)))+eps));
 H=theta;
 H(B>G)=2*pi-H(B>G);
 H=H/(2*pi);
 S=1-3.*(min(min(R,G),B))./(R+G+B+eps);
 I=(R+G+B)/3;
 hsi = cat(3,H,S,I); %HSI image
 %Y,Cr,Cb
 %YUV = [0.299 0.587 0.114; -0.169 -0.331 0.5; 0.5 -0.419 -0.081];
 Y = 0.299.*R + 0.587.*G + 0.114.*B;
 Cb = -0.169.*R - 0.331.*G + 0.5.*B;
 Cr = 0.5.*R - 0.419.*G - 0.081.*B;
 YCbCr = cat(3,Y,Cb,Cr); %YCbCr image
 
 %% Displaying figures (Edit this part as needed)

 figure('Name','DIP', 'NumberTitle','off'); imshow(color_image,[]); title('Original Color\_monarch\_512x512 (RGB)' ); % display the origional image
 figure('Name','DIP', 'NumberTitle','off');
 subplot(3,3,1); imshow(r,[]); title('(R)' );
 subplot(3,3,2); imshow(g,[]); title('(G)' );
 subplot(3,3,3); imshow(b,[]); title('(B)' );
 subplot(3,3,4); imshow(H,[]); title('(H)' );
 subplot(3,3,5); imshow(S,[]); title('(S)' );
 subplot(3,3,6); imshow(I,[]); title('(I)' );
 subplot(3,3,7); imshow(Y,[]); title('(Y)' );
 subplot(3,3,8); imshow(Cb,[]); title('(Cb)' );
 subplot(3,3,9); imshow(Cr,[]); title('(Cr)' );

 %%---------------------------------------------------------------
end
