%Lab 3 Digital Signal Processing

function Lab3() %Function Defination 



 %*****************************************************************%   
%Exercise 1 : 1D Convolution
      %1.2 a  
         x = [ 1 2 3 4 ] ;
         y_s1 = Dirac(5,10) ;
         y1 = convolution (x,y_s1);
         figure(1) 
         stem(y1);
         xlabel('n')
         ylabel(' y = x * h  ')
         title('Convolution Function') 
      
      %1.2 b  
         h_s2 = [ 0 0 0 0 1 ] ;
         y2 = convolution (x,h_s2); %Convolution Function
         figure(2) 
         stem(y2);
         xlabel('n')
         ylabel(' y = x * h  ')
         title('Convolution Function')           
         
      %1.3 
         y_s3 = [exp(1-5) exp(2-5) exp(3-5) exp(4-5)];
         y3 = convolution (x,y_s3);
         figure(3) 
         stem(y3);
         xlabel('n')
         ylabel(' y = x * h  ')
         title('Convolution Function') 
         
      %1.4
         h_s4 = [ -1 1 ] ;
         y4 = convolution (x,h_s4); %Convolution Function
         figure(4) 
         stem(y4);
         xlabel('n')
         ylabel(' y = x * h  ')
         title('Convolution Function')   
         
%*****************************************************************%

%*****************************************************************%   
%Exercise 2 : 2D Convolution
      %2.1 a = Creating 2D convolution Function  
      %2.2 Using your function, smooth the \ lena" image with the following Gaussian kernel.

      lena_image=imread('E:\AVibot 2016-2018\Digital Signal processing\Lab3\lena-grey.bmp');
      lena_image=mat2gray(lena_image); %Converting matrix to intensity image The returned matrix I contains values in the range 0.0 (black) to 1.0 (full
                                       % intensity or white)
      %Applying Mask or Kernel Given
      kernel=1/256.*[1 4 6 4 1;4 16 24 16 4;6 24 36 24 6;4 16 24 16 4;1 4 6 4 1]; 
      [smooth_lena_image]=convolution2D(lena_image,kernel);
      
      figure(5);
      imshow(lena_image) %Showing lena image
      title('lena Image')
      figure(6) 
      imshow(smooth_lena_image) %Showing smooth lena image after applying kernel
      title('Lena Image after smoothing through Kernel')
      
      %2.3 Apply sobel Filtering using convolution.
      sobelfilter_x = [-1 0 1; -2 0 2; -1 0 1]; % we already know the kernel of sobel filter
      sobelfilter_y = [-1 -2 -1; 0 0 0; 1 2 1];
      [sobelfilter_lena_x]=convFn2(lena_image,sobelfilter_x);
      [sobelfilter_lena_y]=convFn2(sobelfilter_lena_x,sobelfilter_y);
     
      figure(7)
      imshow(sobelfilter_lena_y) %Showing image after convolving through sobel filter
      title('Lena image after applying Sobel Filter')

%*****************************************************************%   
%Exercise 3 : Character recognition using cross-correlation     
      % 3.1
      First_image_A=imread('E:\AVibot 2016-2018\Digital Signal processing\Lab3\a.png'); %Reading Image
      First_image_A=mat2gray(First_image_A); %Converting matrix to intensity image The returned matrix I contains values in the range 0.0 (black) to 1.0 (full
                                       % intensity or white)
      
      second_text_image=imread('E:\AVibot 2016-2018\Digital Signal processing\Lab3\text.png');
      second_text_image=mat2gray(second_text_image);%Converting matrix to intensity image The returned matrix I contains values in the range 0.0 (black) to 1.0 (full
                                       % intensity or white)
                      
      %Applying otsu_threshold function to seperate background and foreground  
      %Image A
      otsu_threshold_First_image_A=graythresh(First_image_A);
      binary_image_A = im2bw(First_image_A, otsu_threshold_First_image_A); %Convert image to binary image by thresholding.
      binary_image_A =+ binary_image_A;
      binary_image_A= ones(16,16)- binary_image_A;
     
     %Image text
      otsu_threshold_second_text_image=graythresh(second_text_image); %Uisng in built Matlab graythresh function 
      binary_image_text = im2bw(second_text_image, otsu_threshold_second_text_image);
      binary_image_text =+ binary_image_text;
      
      %Showing image A
      figure(8)
      imshow(binary_image_A)
      title('Image a after applying otsu_threshold ')
      
      %Showing image text
      figure(9)
      imshow(binary_image_text)
      title('Image text after applying otsu_threshold ')
      
      %3.2
      %Use signal.correlate2d and xcorr in python and matlab, respectively to correlate your binarize images.
      xcorr_image_A_and_text = xcorr2(binary_image_text,binary_image_A); %Corelating binarized images using function xcorr2
      intensity_image=mat2gray(xcorr_image_A_and_text);%Converting matrix to intensity image The returned matrix I contains values in the range 0.0 (black) to 1.0 (full
                                       % intensity or white)
      [row_max, col_max]= find(intensity_image ==max(max(intensity_image)));
      
      figure(10)
      imshow(intensity_image)
      title('Correlation of binarized images')
     
      figure(11)
      imshow(second_text_image)
      hold on;
      %we plot the image with the finded 'a' letters
      plot(col_max-6,row_max-6, 'g*','MarkerSize',6) % To adjust the center of letter a in the txt imag we subtract few values 
      title('Detected letter a on text image')
         
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%

%Convolution Function
function y = convolution (x,h)

x2 = h ; 
length_x = length(x); %Calculating length of the input signal x
length_h = length(h); %Calculating length of the input signal h

if length_x > length_h  %Applying Conditions
x2 = [x2 zeros(1,length_x - length_h)]; 
else 
x=[x zeros(1,length_h - length_x)]; 
end 

y=zeros(1,length_x + length_h-1);  %padding zeros

x=x(end:-1:1);  % Fliiping input signal x 

for i=1:length(y) 
if i<=length(x) 
y(i)=sum(x(1,length(x)-i+1:length(x)).*x2(1,1:i)); 
else 
j=i-length(x); 
y(i)=sum(x(1,1:length(x)-j).*x2(1,j+1:length(x2))); 
end 
end 
       


end



%Function 2 
function y1 =  Dirac(n,N) %Function Defination 

    if ((n<1)||(n>N))
            disp('Error : n should be inferior then N-1');  %Display error if n > N-1
            y1= 0;
    else
            s = zeros(1,N);  
            s(n) = 1 ;
            y1 = s;
           
    end
  
end


%Function 3 
function [convolution_matrix] = convolution2D(input_1,input_2) %Function convolution in 2D

flipped_signal = flipdim(flipdim(input_2,1),2); %Flipping  the signal

[row_size_in1,col_size_in1] = size(input_1); %Calculating size of input 1 signal
[row_size_in2,col_size_in2] = size(input_2); %Calculating size of input 2 signal
adding_row = 2*floor((row_size_in2/2));

adding_matrix = zeros(row_size_in1 + adding_row, col_size_in1 + adding_row); %Padding zeros
[ row_added_matrix, column_added_matrix] = size(adding_matrix); %Calculating size
adding_matrix(((adding_row/2) + 1):( row_added_matrix-( adding_row/2) ),( (adding_row/2)+1):( column_added_matrix-(adding_row/2) ))=input_1;
convolution_matrix=zeros(row_size_in1,col_size_in1); %Padding zeros

%Convolving two signal through double loop
for (i=((adding_row/2)+1):(row_added_matrix-(adding_row/2)))
    for (j=((adding_row/2)+1):(column_added_matrix-(adding_row/2)))
        
        convolution_matrix(i-adding_row/2,j-adding_row/2)=sum(sum(adding_matrix(i-adding_row/2:i+adding_row/2, j-adding_row/2: j+ adding_row/2).*flipped_signal));
                
    end
    
end
end
%%%%%********************************************************************%%%%