%CompE 565: Multimedia Communication Systems
%Project 1: Basic Digital Image Processing Operations
%Vyshak Athreya BK (51523409) , Vishakha Vijaykumar (820535151)

clear all;
clc;

% Question 1
% Reading the image
%Read and display the image using Matlab
img = imread('D:\Course work\MCS\Assignments\Landscape.jpg');
% Displaying the image
figure(1);
imshow(img);
title('Original RGB image');

% Question 2 
% Display each band (Red, Green and Blue) of the image file
%Display each band (Red, Green and Blue) of the image file
z = zeros(size(img, 1), size(img, 2));
red_comp=img(:,:,1);
red_img = cat(3, red_comp, z, z);
green_comp=img(:,:,2);
green_img = cat(3,z,green_comp,z);
blue_comp=img(:,:,3);
blue_img = cat(3,z,z,blue_comp);
figure(2);
subplot(2,2,1);
imshow(img);
title('Original RGB image');
subplot(2,2,2);
imshow(red_img);
title('Red Component');
subplot(2,2,3);
imshow(green_img);
title('Green Component');
subplot(2,2,4);
imshow(blue_img);
title('Blue Component');

% Question 3
%Convert the image into YCbCr color space
ycbcr=rgb2ycbcr(img);
img_rgb=ycbcr2rgb(ycbcr);

% Question 4
%Display Y, Cb and Cr bands separately 
figure(3);
subplot(2,2,1);
imshow(ycbcr);
title('YCbCr image');
Y_component=ycbcr(:,:,1); % to get Y component
Cb_component=ycbcr(:,:,2); % to get Cb component
Cr_component=ycbcr(:,:,3); % to get Cr component
subplot(2,2,2);
imshow(Y_component);
title('Y Component');
subplot(2,2,3);
imshow(Cb_component);
title('Cb component');
subplot(2,2,4);
imshow(Cr_component);
title('Cr component');

% Question 5
% Subsample Cb and Cr bands using 4:2:0 and display both bands.
Cr_copy = Cr_component;
Cb_copy = Cb_component;

original = size(Cr_copy)
r = size(Cr_copy,1);
c = size(Cr_copy,2);

% Making alternate rows and columns zero
for row = 2:2:r
    for col = 2:2:c
    Cb_copy(row,col) = 0;
    Cr_copy(row,col) = 0;
    end
end

Cb_ip = Cb_copy;
Cr_ip = Cr_copy;
Cb_replica = Cb_copy;
Cr_replica = Cr_copy;

% Removing alternate rows and columns since we dont send zeros
Cb_copy(2:2:r,:) = [];
Cb_copy(:,2:2:c) = [];
Cr_copy(2:2:r,:) = [];
Cr_copy(:,2:2:c) = []; 

figure(4);
subplot(2,2,1);
imshow(Cb_component);
title('Original Cb component');
subplot(2,2,2);
imshow(Cb_copy);
title('4:2:0 subsampled Cb component');
subplot(2,2,3);
imshow(Cr_component);
title('Original Cr component');
subplot(2,2,4);
imshow(Cr_copy);
title('4:2:0 subsampled Cr component');
compr=size(Cr_copy)

% Question 6
% 6.1 Upsampling using Liner Interpolation
    % Averaging the rows and coloumns
    if mod(r,2) == 0
        for row = 2:2:r-2
        for col = 2:2:c-2
                Cb_ip(row,col) = round(Cb_ip(row-1,col-1)/2 + Cb_ip(row+1,col+1)/2);
                Cr_ip(row,col) = round(Cr_ip(row-1,col-1)/2 + Cr_ip(row+1,col+1)/2);
        end    
        end
        for row = r
            for col = c
            Cb_ip(row,col) = Cb_ip(r-1,c-1);
            Cr_ip(row,col) = Cr_ip(r-1,c-1);
            end       
        end     
    else
        for row = 2:2:r-1
            for col = 2:2:c-1    
            Cb_ip(row,col) = round((Cb_ip(row-1,col-1) + Cb_ip(row+1,col+1))/2);
            Cb_ip(row,col) = round((Cb_ip(row-1,col-1) + Cb_ip(row+1,col+1))/2);
            end
        end
    end
    
    
  
%6.2 Upsample and display the Cb and Cr bands using Simple row or column replication
%First replicate the missing entries in the odd rows. Then replicate the
%columns

 for row = 2:2:r
    Cb_replica(row,:) = Cb_replica(row-1,:);
    Cr_replica(row,:) = Cr_replica(row-1,:);
    for col = 2:2:c
    Cb_replica(:,col) = Cb_replica(:,col-1);
    Cr_replica(:,col) = Cr_replica(:,col-1);
    end
 end
     
figure(5);
subplot(2,2,2);
imshow(Cb_ip);
title('Up sampled Cb by Linear Interpolation');
subplot(2,2,4);
imshow(Cr_ip);
title('Up sampled Cr by Linear Interpolation');
subplot(2,2,1);
imshow(Cb_copy);
title('Compressed Cb');
subplot(2,2,3);
imshow(Cr_copy);
title('Compressed Cr');
Us_Li=size(Cr_ip)

figure(6);
subplot(2,2,2);
imshow(Cb_replica);
title('Up sampled Cb by Column Replication');
subplot(2,2,4);
imshow(Cr_replica);
title('Up sampled Cr by Column Replication');
subplot(2,2,1);
imshow(Cb_ip);
title('Compressed Cb');
subplot(2,2,3);
imshow(Cr_copy);
title('Compressed Cr'); 
 
 % Question 7
 % Converting the image into RGB format
 ycbcr_interpol(:,:,1) = Y_component;
 ycbcr_interpol(:,:,2) = Cb_ip;
 ycbcr_interpol(:,:,3) = Cr_ip;
 rgb_interpol = ycbcr2rgb(ycbcr_interpol);
  
 ycbcr_replica(:,:,1) = Y_component;
 ycbcr_replica(:,:,2) = Cb_replica;
 ycbcr_replica(:,:,3) = Cr_replica;
 rgb_replica = ycbcr2rgb(ycbcr_replica);
 
 % Question 8
 % To display the original and reconstructed images restored from the YCbCr coordinate
 figure(7);
 subplot(2,1,1);
 imshow(img);
 title('Original image');
 subplot(2,1,2);
 imshow(rgb_interpol);
 title('Reconstructed image after Linear Interpolation');
 
 figure(8);
 subplot(2,1,1);
 imshow(img);
 title('Original image');
 subplot(2,1,2);
 imshow(rgb_replica);
 title('Reconstructed image after Row Column Replication');
 
% Question 9
%From the above images it can be observed that we get a better
%reconstructed image by using the method of linear interpolation than row
%column replication. Since linear interpolation method uses the concept of
%taking the average of neighbouring cells, it gives a better pixel value
%thereby improving the image quality. In the above image it is evident that
%the freqency is very high. There are many pixels of varied intensity
%than that of similar adjacent pixels. In row column replication we
%reproduce the adjacent values and hence it is not effecient here.

% Question 10
% MSE between the original and reconstructed images
Difference = (img - rgb_interpol).^2;
S = sum(sum(Difference));
MSE = S/(r*c)
%Mean Square error is the average of the squared difference between the
%pixels and is an objecive measure of the signal quality. Here it ranges
%from 0.1505 to 0.8086 . We can conclude that the up sampled signals are
%good and image reconstruction is effecient.
 
% Question  11
%Post compression of the Cb and Cr components, the image size was reduced
%to 240*320 from 480*640. This shows the compression ratio is 1:2. This
%ratio indicates that we have reduced the transmission data by half of the
%actual data and henceforth the gain is appreciable.
 