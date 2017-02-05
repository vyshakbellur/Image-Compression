%CompE 565: Multimedia Communication Systems
%Project 2: JPEG based Image Compression
%Vyshak Athreya BK (51523409) , Vishakha Vijaykumar (820535151)

clear all;
clc;

%At the encoder using 4:2:0 YCbCr component image,
%(a) Compute the 8x8 block DCT transform coefficients of the luminance and chrominance components of the image.
%Please display the DCT coefficient matrix as well as image of the DCT transformed image blocks of the first 
%2 blocks in the 5th row (of blocks) from top for the luminance component.

img =imread('Landscape.jpg');
r=size(img,1);
c=size(img,2);
% Displaying the image
figure(1);
imshow(img);
title('Original RGB image');

%Convert the image into YCbCr color space
ycbcr=rgb2ycbcr(img);

Y_component=ycbcr(:,:,1); % to get Y component
Cb_component=ycbcr(:,:,2); % to get Cb component
Cr_component=ycbcr(:,:,3); % to get Cr component

%4:2:0 subsample Cb and Cr component
Cb_subsample = Cb_component(1:2:r,1:2:c); 
Cr_subsample = Cr_component(1:2:r,1:2:c); 
figure(2);
subplot(2,1,1), subimage(Cb_subsample), title('subsampled Cb element'); % subsampled Cb element is displayed
subplot(2,1,2), subimage(Cr_subsample),  title('subsampled Cr element');% subsampled Cr element is displayed
 
DCT_y = @dct2;  % to obtain DCT coefficients for Y component
DCT_coef_Y = blkproc(Y_component, [8 8], DCT_y); % block processing  for 8*8blocks where  the DCT function is performed
DCT_coef_Y = fix(DCT_coef_Y); % each element is rounded to the nearest integer towards zero
 
DCT_Cb = @dct2; %dct2 command is used for doing DCT for Cb element
DCT_coef_Cb = blkproc(Cb_subsample, [8 8], DCT_Cb); % block processing  for 8*8blocks the DCT function is performed
DCT_coef_Cb = fix(DCT_coef_Cb);%each element is rounded to the nearest integer towards zero
 
DCT_Cr = @dct2; %performing the DCT using the command dct2 for Cr element
DCT_coef_Cr = blkproc (Cr_subsample, [8 8], DCT_Cr); %performing block processing operation for 8*8blocks where in interior the DCT function is performed
DCT_coef_Cr = fix(DCT_coef_Cr);%rounds each element to nearest intezer toward zero
 
figure(3); % Y, Cb , Cr DCT elements are displayed
subplot(3,1,1), subimage(DCT_coef_Y), title('DCT of  Y element');
subplot(3,1,2), subimage(DCT_coef_Cb), title('DCT  of Cb element');
subplot(3,1,3), subimage(DCT_coef_Cr), title('DCT of Cr element');
 
%   DCT coefficient matrix as well as image of the DCT transformed image blocks of the 1st and 2nd block in the 5th row from top for the luminance element. 
Y_DCT_block1_row5 = zeros(8,8); %to extract the 5th row first block 8*8matrix
Y_DCT_block1_row5(1:8,1:8) = DCT_coef_Y(33:40,9:16);

Y_DCT_block2_row5 = zeros(8,8); %to extract the second block 8*8 matrix
Y_DCT_block2_row5=DCT_coef_Y(33:40,1:8);
figure(5); % forth and fifth 8*8 blocks of 4th row are displayed
subplot(2,1,1), subimage(Y_DCT_block1_row5), title('first block in the 5th row of DCT transformed image');
subplot(2,1,2), subimage(Y_DCT_block2_row5), title('second block in the 5th row of DCT transformed image');
 
% (b) Quantize the DCT image by using the JPEG luminance and chrominance quantizer matrix from the lecture notes.
% Report the following output only for the first 2 blocks in the 5th row from top of the luminance component: 
% (a) DC DCT coefficient; (b) Zigzag scanned AC DCT coefficients.
 
Y_Quant_matrix = [16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 14 13 16 24 40 57 69 56; 14 17 22 29 51 87 89 62; 18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 49 64 78 87 108 121 120 101; 72 92 95 98 112 100 103 99];
 
quant_y = @(DCT_Yelement) round(DCT_Yelement ./ Y_Quant_matrix);   % Quantization is done for the luminance element (Y) by using the quantization matrix     
quantized_Y = blkproc(DCT_coef_Y,[8 8],quant_y); %performing block processing operation for 8*8blocks
 
Cr_Quant_matrix = [17 18 24 47 99 99 99 99; 18 21 26 66 99 99 99 99; 24 26 56 99 99 99 99 99; 47 66 99 99 99 99 99 99; 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99];
 
quant_Cb = @(DCT_Cbelement) round(DCT_Cbelement ./ Cr_Quant_matrix); %Quantization is done for the chrominance element (Cb) by using the quantization matrix   
quantized_Cb = blkproc(DCT_coef_Cb,[8 8],quant_Cb); %performing block processing operation for 8*8blocks
 
quant_Cr = @(DCT_Crelement) round(DCT_Crelement ./ Cr_Quant_matrix); %Quantization is done for the chrominance element (Cr) by using the quantization matrix 
quantized_Cr = blkproc(DCT_coef_Cr,[8 8],quant_Cr); %performing block processing operation for 8*8blocks
 
% Report the following output only for the first 2 blocks in the 5th row from top of the luminance component: 
%(a) DC DCT coefficient; (b) Zigzag scanned AC DCT coefficients.
 
Y_DCT_block1_row5 = zeros(8,8); %extracting the fifth row 8*8 block of quantized matrix
Y_DCT_block1_row5(1:8,1:8) = quantized_Y(33:40,1:8)
figure(51);
subplot(2,1,1); 
imshow(Y_DCT_block1_row5);
title('block 1 row 5');
 
Y_DCT_block2_row5 = zeros(8,8); %extracting the fifth 8*8 block of quantized matrix
Y_DCT_block2_row5(1:8,1:8) = quantized_Y(25:32,33:40)
subplot(2,1,2); 
imshow(Y_DCT_block2_row5);
title('block 2 row 5');

DCcoeff_block1_row5 =  Y_DCT_block1_row5(1,1); % (1,1) element of  forth block quantized matrix
DC_coeff_block2_row5 = Y_DCT_block2_row5(1,1);% (1,1) element of fifth block quantized matrix
 
message1 = sprintf('DCT DC coefficient of the first block in the 5th Row of the quantized image is %.2f\n',DCcoeff_block1_row5) %displaying DC coefficient 
       
message2 = sprintf('DCT DC coefficient for the second block in the 5th Row of the quantized image is %.2f\n',DC_coeff_block2_row5) %displaying  DC coefficient 
 
% zigzag scan for first block in the fifth row
% use four variables, one to switch between rows and columns , other to
% move across the diagonals and the rest two to get the x and y positions
 
temp1=0;
size_block1_row5=size(Y_DCT_block1_row5);
sum_block1=size_block1_row5(2)*size_block1_row5(1);  %Size of the first block computed by multiplying columns and rows
for i1=2:sum_block1
 j1=rem(i1,2);  % to check whether the value is even or odd 
    for counti1=1:size_block1_row5(1)
        for countj1=1:size_block1_row5(2)
            if((counti1+countj1)==i1)
                temp1=temp1+1;
                if(j1==0)
                zigzagblock5_1(temp1)= Y_DCT_block1_row5(countj1,i1-countj1);
                else          
                zigzagblock5_1(temp1)= Y_DCT_block1_row5(i1-countj1,countj1);
                end
             end    
         end
     end
end
Y_DCT_block1_row5; 
disp(zigzagblock5_1);
 
% zigzagscan for second block in the fifth row

temp2=0; 
size_block2_row5=size(Y_DCT_block2_row5);
sum_block2=size_block2_row5(2)*size_block2_row5(1);  % M*N value is calculated
for i2=2:sum_block2
 j2=rem(i2,2);  % to check whether the value is even or odd 
    for count_2i=1:size_block2_row5(1)
        for count_2j=1:size_block2_row5(2)
            if((count_2i+count_2j)==i2)
                temp2=temp2+1;
                if(j2==0)
                zigzagblock5_2(temp2)=Y_DCT_block2_row5(count_2j,i2-count_2j);
                else          
                zigzagblock5_2(temp2)=Y_DCT_block2_row5(i2-count_2j,count_2j);
                end
             end    
         end
     end
end
Y_DCT_block2_row5; 
disp(zigzagblock5_2);
 
%DECODER
% (c) Compute the inverse Quantized images obtained in Step (b).
 
Inv_Quant_y = @(quantized_Y) round(quantized_Y .* Y_Quant_matrix);   % Inverse quantization for luminance element 
Inverse_Quantized_Y = blkproc(quantized_Y,[8 8],Inv_Quant_y); % block processing operation
 
InverseQuant_Cbelement = @(quantized_Cb) round(quantized_Cb .* Cr_Quant_matrix);  % inverse quantization for Cb element
Inverse_Quantized_Cb = blkproc(quantized_Cb,[8 8],InverseQuant_Cbelement); % block processing operation 
 
InverseQuant_Crelement = @(quantized_Crelement) round(quantized_Crelement .* Cr_Quant_matrix);  % inverse quantization for Cr element
Inverse_Quantized_Cr = blkproc(quantized_Cr,[8 8],InverseQuant_Crelement); %block processing operation
 
figure(6); 
subplot(3,1,1), subimage(Inverse_Quantized_Y), title('Y element of inverse quantized image') 
subplot(3,1,2), subimage(Inverse_Quantized_Cb), title('Cb element of inverse quantized image') 
subplot(3,1,3), subimage(Inverse_Quantized_Cr), title('Cr element of inverse quantized image') 
 
% (d) Reconstruct the image by computing Inverse DCT coefficients. 
 
a_DCT = @idct2;  % inverse DCT operation 
IDCT_Y = blkproc(Inverse_Quantized_Y, [8 8], a_DCT); % block processing operation on 8*8 Y block
IDCT_Y = uint8(fix(IDCT_Y));
 
b_DCT = @idct2; % inverse DCT operation 
IDCT_Cb = blkproc(Inverse_Quantized_Cb, [8 8], b_DCT); % block processing operation on 8*8 Cb block
IDCT_Cb = uint8(fix(IDCT_Cb));
 
c_DCT = @idct2; % inverse DCT operation 
IDCT_Cr = blkproc(Inverse_Quantized_Cr, [8 8], c_DCT); % block processing operation on 8*8 Cr block
IDCT_Cr = uint8(fix(IDCT_Cr));
 
figure(7); 
subplot(3,1,1), subimage(IDCT_Y), title('Y element of inverse DCT image'); 
subplot(3,1,2), subimage(IDCT_Cb), title('Cb element of inverse DCT image'); 
subplot(3,1,3), subimage(IDCT_Cr); title('Cr element of inverse DCT image'); 
 
% ROW COLUMN REPLICATION
 
temp2=zeros(r,c,3); % Temporary 3 dimensional matrix with element zeros
temp2(1:2:r,1:2:c,2)=IDCT_Cb(:,:); % Cb subsampled matrix elements to the temporary matrix
temp2(1:2:r,1:2:c,3)=IDCT_Cr(:,:); % Cr subsampled matrix elements to the temporary matrix
Cb_upsample=temp2(:,:,2);
Cr_upsample=temp2(:,:,3); 
Cb_upsample=uint16(Cb_upsample); 
Cr_upsample=uint16(Cr_upsample);
for i=1:2:r 
    for j=2:2:c
        Cb_upsample(i,j)=Cb_upsample(i,j-1);
        Cr_upsample(i,j)=Cr_upsample(i,j-1);
    end
end
for i=2:2:r % performing upsampling using rowcolumn replication for even number of i for Cb,Cr ie., assign previous row element to respective row element
    for j=1:c
        Cb_upsample(i,j)=Cb_upsample(i-1,j);
        Cr_upsample(i,j)=Cr_upsample(i-1,j);
    end
end
Cb_upsample=uint8(Cb_upsample); 
Cb_upsample=uint8(Cb_upsample); 
 
ycbcr__upsample(:,:,1)=IDCT_Y; % assigning y element
ycbcr__upsample(:,:,2)=Cb_upsample; % assigning the Cb element
ycbcr__upsample(:,:,3)=Cr_upsample; %assigning the Cr element
 
% RGB IMAGE
 
rgb__upsample=ycbcr2rgb(ycbcr__upsample); %converting ycbcrimage obtained by row column upsampling to RGB image
figure(8);
imshow(rgb__upsample);
title('Reconstructed RGB image');
 
figure(9);
subplot(2,1,1), subimage(img), title('Original image'); %Subplotting the original rgb image and rgbimg obtained from rowcolumn replication upsampling
subplot(2,1,2), subimage(rgb__upsample),title('Reconstrructed image');
 
% Display the Error Image (by subtracting the reconstructed image form the original) for the luminance image. 
 
Error_luminance = Y_component - IDCT_Y; % difference between initial y elemennt and inverse DCT  of Y element
figure(10);
imshow(Error_luminance);
title('error image');% error image of luminance 
 
% MEAN SQUARE ERROR OF Y element
Difference = (Y_component - IDCT_Y).^2;
S = sum(sum(Difference));
MSE = S/(r*c) 
Mse=sprintf ('The Mean Square Error for luminance component is %.2f\n',MSE)
% PSNR OF luminance(y)
PSNR_Y1=10*(log10(((255)^2)/MSE)); % Calculation of PSNR
psnr_Y = sprintf('PSNR for luminance component = %.2f dB\n',PSNR_Y1)
