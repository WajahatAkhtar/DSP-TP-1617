
%Lab 5 Digital Signal Processing : 
%Wajahat Akhtar.


%*****************************************************************%   
% Exercice 1 : 2D - DFT
img = zeros(301,301);  %Creating an image with zeros 
img(100:200, 140:160) = 255;  %Adding a white box in the image above by setting so,e rows and coulmns to one
figure(1); %Displaying the figure
subplot(131) 
imshow(img);

%*****************************************************************%  

%*****************************************************************%   
%we can compute the FFT as follows:
imgFreq = fftshift(fft2(img)); %Taking fast fourier transform
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude') %Plotting magnitude 
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase') %Plotting Phase 

%*****************************************************************% 

%*****************************************************************%   
%Compute the FFT of translated, and rotated image
img = zeros(301,301);    %Creating an image with zeros 
img(100:200, 140:160) = 255;  %Adding a white box in the image above by setting so,e rows and coulmns to one

imgTrans = zeros(301,301);  %translating the image
imgTrans(150:250, 160:180) = 255;

imgRot = imrotate(img, 45); %Rotating the image with an angle of 45 degree

% Creating different images and checking the difference in phase and magnitude after taking fft. 
img2 = zeros(301,301); 
img2(20:120, 140:160) = 255;
img2(180:280, 140:160) = 255;

img3 = zeros(301,301); %Third image with same size
img3(100:200, 145:155) = 255;

%Taking fourier transform all the images 
imgFreq1 = fftshift(fft2(imgTrans)); %FFT of translated image
imgFreq2 = fftshift(fft2(imgRot));  %FFT of rotated image
imgFreq3 = fftshift(fft2(img2));  %FFT of image2
imgFreq4 = fftshift(fft2(img3));  %FFT of iamge 3


figure(2);
subplot(131); imshow(imgTrans); %Ploting Translated image
subplot(132); imagesc(abs(imgFreq1)); colormap('gray'); title('Magnitude') %Plotting magnitude 
subplot(133); imagesc(angle(imgFreq1)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

figure(3);
subplot(131); imshow(imgRot);%Ploting rotated image
subplot(132); imagesc(abs(imgFreq2)); colormap('gray'); title('Magnitude') %Plotting magnitude 
subplot(133); imagesc(angle(imgFreq2)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

figure(4)
subplot(131);   imshow(img2); %Ploting image2
subplot(132); imagesc(abs(imgFreq3)); colormap('gray'); title('Magnitude')%Plotting magnitude 
subplot(133); imagesc(angle(imgFreq3)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

figure(5)
subplot(131);   imshow(img3);%Ploting image3
subplot(132); imagesc(abs(imgFreq4)); colormap('gray'); title('Magnitude')%Plotting magnitude 
subplot(133); imagesc(angle(imgFreq4)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

%*****************************************************************%   

%*****************************************************************%   
%Considering Synthetic image and compute the phase and the magnitude of the
%following image and displaying the normalized center frequency.
Im=0; 
N=64;
T=1;
Ts=T/N;
Fs=1/Ts;
df=Fs/N;
Im(N/8:N/4,N/4+1:N/2)=1; 
Im(1:N/4,N/2+1:N)=Im;
Im(N/4+1:N/2,:) = Im;
Im(N/2+1:3*N/4,:) = Im(1:N/4,:);
Im(3*N/4+1:N,:) = Im(1:N/4,:);

imgFreq1 = fftshift(fft2(Im)); %Taking fast fourier transform of the image created above

%Normalized central Frequency by finding mean
Centeral_frequency = mean (Im(:))  ;
disp('Centeral_frequency = ') % Displaying central Frequency by finding mean
disp(Centeral_frequency);

figure(6);
subplot(1,3,1); imshow(Im)
subplot(1,3,2); imagesc(abs(imgFreq1)); colormap('gray'); title('Magnitude')%Plotting magnitude 
subplot(1,3,3); imagesc(angle(imgFreq1)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

%*****************************************************************%   

%*****************************************************************%   
% 1.5 Plot If(u,0) and If(0,v) with the correct frequency range. 

central_row = imgFreq1(N/2+1,:); %Taking FFT of central_row 
central_column = imgFreq1(:,N/2+1); % Taking FFT of central_column
center = (-N/2 : N/2-1);

figure (7); 
subplot(1,2,1); plot(center,abs(central_row)); title('If(u,0)')%Plotting magnitude 
subplot(1,2,2); plot(center,abs(central_column)); title('If(0,v)')%Plotting magnitude 

% Comments : We can see after plotting we figure 7 above we see two
% profiles one for central row of 2D image and one for central coulmn of 2D
% image with If(u,0) both having higher intensity at the center with small
% difference in frequency.

%*****************************************************************%   

%*****************************************************************%   
%1.6 Load the \lena" image, and show the phase and magnitude, then reconstruct the image using
%either frequency or phase.

figure(8);
I = imread ('lena-grey.bmp');
%I= rgb2gray(I);
imgFreq_lena = fftshift(fft2(I));
subplot(1,3,1); imshow(I)
subplot(1,3,2); imagesc(log(abs(imgFreq_lena))); colormap('gray'); title('Magnitude')%Plotting magnitude 
subplot(1,3,3); imagesc(angle(imgFreq_lena)/pi*180); colormap('gray'); title('Phase')%Plotting Phase

%In order to reconstrcut lena image we can just take inverse fast fourier
%transform and we dont need to shift to the center as we are already
%shifting to the center by taking fftshift above so in order to reconstruct
%we simply take ifft as below but in case we take fftshift so then we have
%to move to the center of the image in order to reconstruct the image.

figure(9);
imgFreq_lena2 = ifft2(imgFreq_lena);
subplot(1,2,1); imshow(I);
subplot(1,2,2); imagesc(log(abs(imgFreq_lena2))); colormap('gray'); title('Reconstructed image from Frequency')%Plotting image 

%*****************************************************************%   

%*****************************************************************%   
%1.7 Apply the sobel filter only in vertical direction to `lena` image in the frequency domain.
sobel_filter_matrix=[-1 0 1; -2 0 2; -1 0 1]; %Sobel Filter we already know that 
padding_zeros = size(I) + size(sobel_filter_matrix) - 1; %Finding the size for padding zeros
image_I = double(I); %Changing data type of image ;

%Compute the DFT to obtain H(u; v)
DFT_I = fft2(image_I,padding_zeros(1), padding_zeros(2)); %padding zeros and taking DFT
DFT_S= fft2(double(sobel_filter_matrix), padding_zeros(1), padding_zeros(2)); %padding zeros and taking DFT in order to convert spatial filter to frequency filter

%Apply the multiplication in the Fourier space
DFT_Conv = DFT_S.*DFT_I;  %In Frequency domain we multiply for convolution
%Compute the inverse Fourier transform
FFT_inv = ifft2(DFT_Conv); 
%Crop the image at its original size
FFT_inv = FFT_inv(2:size(I,1)+1, 2:size(I,2)+1);

figure(10);
imshow(FFT_inv,[]);%Displaying Fianl Image after all operations

%*****************************************************************%   
                       % END
%*****************************************************************%   











