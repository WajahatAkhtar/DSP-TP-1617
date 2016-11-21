
% Wajahat Akhtar : Lab 4 Digital Signal Processing
%Discrete Fourier Transform

function lab4() %Function Defination 

clc ;
clear all ;

%*****************************************************************%  
%Exercice 1 : DFT

f = 5; fs = 50;  %Defining Sampling frequency of the function
t = 0: 1/fs : 1; 
xn1 = sin(2 * pi *f * t);  %Simple sin function
N = length(xn1);
fr = (-N/2 : N/2-1)* fs/N;
xf1 = fftshift(fft(xn1));  %taking fourier transform of the signal

figure(1);
subplot(221); plot(t, xn1); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf1)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf1)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf1)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');


%1.2
xn2 = cos(2 * pi *f * t);  %Simple cos function
N = length(xn2);
fr = (-N/2 : N/2-1)* fs/N;
xf2 = fftshift(fft(xn2));  %taking fourier transform of the signal

figure(2);
subplot(221); plot(t, xn2); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf2)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf2)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf2)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

% Comments
%In cos and sin signals there is no actual difference in amplitude , there
%lies only a phase shift of 90 degree between both the signals.

%1.3
t = 0: 1/fs : 1;
xn3 = square(2 * pi *f * t);  %Simple square function
N = length(xn3);
fr = (-N/2 : N/2-1)* fs/N;
xf3 = fftshift(fft(xn3));  %taking fourier transform of the signal
figure(3);
subplot(221); plot(t, xn3); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf3)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf3)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf3)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');


%the sqauare wave signal is non-sinusoidal periodic waveform (which can be 
%represented as an infinite summation of sinusoidal waves as een in figure),
%in which the amplitude alternates at a steady frequency between fixed 
%minimum and maximum values, with the same duration at minimum and maximum.


%1.4
t = 0: 1/fs : 1;
xn4 = randn(1,10000);
N = length(xn4);
fr = (-N/2 : N/2-1)* fs/N;
xf4 = fftshift(fft(xn4));
figure(4);
subplot(221); plot( xn4); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf4)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf4)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf4)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

%Gaussian noise is statistical noise distributed over the whole frequency 
%spectrum , which is also known as the Gaussian distribution. it can be
%clearly seen in the figure that gausian noise is distributed on whole
%spectrum

%*****************************************************************%

%*****************************************************************%   
%Exercice 2 : Sampling
%2.1
f1 = 5 ; 
f2 = 20;

fs = [10; 20; 25; 40; 50; 100; 150];

for i=1:length(fs)
    rows = 0: 1/fs(i) : 1;
xn5 = 3*cos(2*pi*f1*rows) + 4*sin(2*pi*f2*rows);

%2.4
N = 1000;
fr5 = (-N/2 : N/2-1)* fs(i)/N;
xf5 = fftshift(fft(xn5,1000));

%2.2 Plotting x[n]
figure;
subplot(221); plot( rows , xn5); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr5, abs(xf5)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr5, real(xf5)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr5, imag(xf5)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

end

% 2.3 and 2.4 Discussing Alliasing effect in time Domain and the FFT
% computation.

% Comments : 
%In signal processing, aliasing is an effect that causes different signals 
%to become indistinguishable (or aliases of one another) when sampled. 
%It also refers to the distortion or artifact that results when the signal 
%reconstructed from samples is different from the original continuous signal.
%if sampling frequency is less then Nyquist frequency which is twice the 
%highest frequency present in the signal then there will be aliasing effect
% if it is twice the highest frequency present in the signal then there
% will be no Aliasing and it  can be clearly seen in the results we got
% there is aliasing effect when sampling frequency is less then 40 hz
% because our signal has a frequency of 20hz so sampling frequency should
% be twice of it that is 40 hz or more in orrder to avoid aliasing.

%*****************************************************************%

%*****************************************************************%   

% Exercise 3: 1D DFT for image classification

% We have already been given some information about barcode tif images and
% non-barcode tif images.
barcode_images_class = 0;  %Creatig class for seperating images with and without barcode
barcode_images = [1,2,6, 44:54] ;  %Given image with barcode
nonbarcode_images = [ 3:5, 7:43] ;  %Given image with no barcode


%Loading the images from 1D-DFT folder
tiff_images_folder = 'E:\AVibot 2016-2018\Digital Signal processing\lab4\images\1D-DFT';  % moving to images folder
cd (tiff_images_folder) % comment this if inside tiff images folder 
tif_images = dir;  %Whole directory data is stored in tif images vector
number_of_tiff_images = length(tif_images)-2; 


%Reading all the images
image_size = zeros(1,number_of_tiff_images);
for  i = 1:number_of_tiff_images
    image_read = imread(tif_images(i+2).name);
     [rows,coulmns]= size(image_read); % Finding the size of the images
    image_size(i)= rows*coulmns; 
      
end

%Reading ad Finding the Image in folder with minimum size
Least_size_image = find(image_size == min(image_size(:)));
image_number = Least_size_image ;
image_read = imread(tif_images(image_number).name);
[coulmns,rows,l] = size(image_read); 

%Normalizing and Resizing Images
barcode_images_var = []; 
non_barcod_images =[];

for i = 1:number_of_tiff_images
    image_read = imread(tif_images(i+2).name);
    [ x, y, z]= size(image_read);
        if z == 4 
            image_read = image_read(:,:,1);
        end
        
% Normalizing the images   
image_read = image_read/(max(image_read(:))); 
    
%Resizing Images
image_resized = imresize(image_read, [coulmns rows]); 


%Taking Fourier transform
N = rows; 
fr = (-N/2 : N/2-1);

first_image = image_resized(coulmns/2,:);
DFT_of_image1 = ifftshift(fft(first_image));

second_image= image_resized( coulmns/2+10,:);
DFT_of_image2 = ifftshift(fft(second_image));

third_image= image_resized( coulmns/2-10,:);
DFT_of_image3 = ifftshift(fft(third_image));


 % Cooments :  
 % The main task of this part of the question is to separate
 % high frequency spectrum from low frequency spectrum that can be done
 % through different techniques but the technique we use is to first filter
 % the images by low pas filter or high pass filter and separate high
 % frequencies from image spectrum if high pas filter is used or low frequencies
 % if low pass filter is used , once we apply filter we take fft which will move 
 % high frequenies on right and low frequencies on left , fft will also 
 % result in increase of magnitude of the spectrums but we want to expand
 % the spectrum in order to separate high frequency class from low
 % frequency class and that we do by simply squaring they spectrum which
 % will furthere increase the spectrum amplitude and separate the classes
 % from each other,and once they are separated we apply that as threshold
 % and the places where there is barcode there will be low frequency
 % spectrum so we add it to class with barcode images and the images
 % without barcode as they will be bright so they are in high frequency
 % spectrum so we will call it as images without tiff images and check
 % the accuracy of our computation in order to get the desired results. 
 


% Setting threshold in order to separate high frequency spectrum images from
% low frequency images.

Im1 = abs( DFT_of_image1);
Im2 = abs ( DFT_of_image2);
Im3 = abs ( DFT_of_image3);


if (  Im1(1,190)  > 4 )
    
    barcode_images_var = [ barcode_images_var, i];
    imge_at_position = find(barcode_images == i);
    
    if ~isempty(imge_at_position)
        barcode_images_class = barcode_images_class + 1;
    end
else
    
    if  (  Im1 - Im2 ) < 0.9 % Changing the value allows you to get better results or near to best results.
        barcode_images_var = [ barcode_images_var, i];
        imge_at_position = find( barcode_images == i);
        
        if ~isempty(imge_at_position)
            barcode_images_class = barcode_images_class + 1; %Adding result to barcode image class
        end
        
    else
        
        non_barcod_images = [ non_barcod_images, i ];
        imge_at_position = find( nonbarcode_images == i);
        if ~isempty( imge_at_position )
            barcode_images_class=barcode_images_class + 1;  %Adding result to non-barcode image class
        end
    end
end
end

%Displaying Results
disp('Barcode Images in the folder');
 disp(barcode_images_var);  
disp('Non-barcode Images in the folder');
disp(non_barcod_images);



%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%



end




