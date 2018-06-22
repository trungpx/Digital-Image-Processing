
function AffineTransformImage2D()
clc;
clear;
close all;
 %% Loading directory for image files
 %RGB image components
 imgdir = uigetdir('Test_images' );
 file = fopen(fullfile(imgdir,'\Color_baboon_256x256.raw' ), 'r' );
 color_image = fread(file,fliplr([256,256*3]), '*uint8' )'; % input matrix for image Color_monarch_512x512.raw
 fclose(file);
 r = color_image(:,1:3:end); %R matrix
 g = color_image(:,2:3:end); %G matrix
 b = color_image(:,3:3:end); %B matrix
 color_image = cat(3, r, g, b);
 inputImage = color_image;
 Result = AffineImage2D(inputImage,0,0,1,1,40); %Input parameter here to translate, rotate, scaling...
 %% Displaying figures (Edit this part as needed)
 figure('Name','DIP', 'NumberTitle','off'); imshow(color_image,[]); title('Color\_baboon\_256x256 (origional)' ); % display the origional image
 figure('Name','DIP', 'NumberTitle','off'); imshow(Result,[]); title('Color\_baboon\_256x256 (Translation)' ); % display the translation image
 %%---------------------------------------------------------------
end

function OUTPUT = AffineImage2D(InputImage, Tr_x, Tr_y, Sc_x, Sc_y, Ang)
%%------------------------------------------------------------------
    T_matrix = [Sc_x*cosd(Ang) Sc_y*sind(Ang) Tr_x;  -Sc_x*sind(Ang) Sc_y*cosd(Ang) Tr_y; 0 0 1];
    X_axis = round(size(InputImage,1)/2);
    Y_axis = round(size(InputImage,2)/2);
    OUTPUT = uint8(ones(size(InputImage)));
    for i = 1:size(InputImage,1)
        for j = 1:size(InputImage,2)
            point =[i-X_axis;j-Y_axis;1] ;
            temp = T_matrix * point;
            temp = int16(temp);
            newX_axis =temp(1,1);
            newY_axis=temp(2,1);
            if newX_axis+X_axis<=size(InputImage,1)&& newX_axis+X_axis >= 1 && newY_axis+Y_axis<=size(InputImage,2)&& newY_axis+Y_axis >= 1 
            OUTPUT(newX_axis+X_axis,newY_axis+Y_axis,:) = InputImage(i,j,:);
            end
        end
    end
%%------------------------------------------------------------------
end

