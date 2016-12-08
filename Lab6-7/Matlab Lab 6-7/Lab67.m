
%Name : Wajahat Akhtar
%Lab 6 and 7
function Lab67() %Function Defination 



% 1.1 Exercice 1 Filters response
% 1.1 - 3rd order Chebyshev Filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            % SECOND PART 3RD ORDER CHEBYSHEV FILTER
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*****************************************************************%  
%Low pass Chebyshev Filter
%*****************************************************************% 
[b,a] = cheby1(3,1,0.5,'low');  
[H, W]= freqz(b,a);
figure(1);
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low Pass - Chebyshev');

%*****************************************************************%  
%High pass Chebyshev Filter
%*****************************************************************% 
[b,a] = cheby1(3,1, 0.5, 'high');
[H, W]= freqz(b,a);
subplot(2,2,2);
plot(W/pi, abs(H));
title('High Pass - Chebyshev');

%*****************************************************************%  
%Band pass Chebyshev Filter
%*****************************************************************% 
[b,a] = cheby1(3,1, [0.2 0.4]);
[H, W]= freqz(b,a);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band-Pass - Chebyshev');

%*****************************************************************%  
%Band Stop Chebyshev Filter
%*****************************************************************% 
[b,a] = cheby1(3,1, [0.2 0.4] , 'stop');
[H, W]= freqz(b,a);
subplot(2,2,4);
plot(W/pi, abs(H)); 
title('Band-Stop - Chebyshev');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % SECOND PART 3RD ORDER BUTTERWORTH FILTER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3rd order  Butterworth Filter
%*****************************************************************%  
%Low pass Butterworth Filter
%*****************************************************************% 
[b,a] = butter(3,0.5,'low');
[H, W]= freqz(b,a);
figure(2);
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low Pass - Butterworth');

%*****************************************************************%  
%High pass Butterworth Filter
%*****************************************************************% 
[b,a] = butter(3,0.5, 'high');
[H, W]= freqz(b,a);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('High Pass - Butterworth');
%*****************************************************************%  
%Band pass Butterworth Filter
%*****************************************************************% 
[b,a] = butter(3,[0.2 0.4]);
[H, W]= freqz(b,a);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band-Pass - Butterworth');
%*****************************************************************%  
%Band Stop Butterworth Filter
%*****************************************************************%   
[b,a] = butter(3,[0.2 0.4] , 'stop');
[H, W]= freqz(b,a);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Band-Stop - Butterworth');
%*****************************************************************%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % LOW PASS CHEBYSHEV WITH VARYING ORDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2 - Low Pass Chebyshev with varying order = 3
%*****************************************************************%  
figure(3);
[b,a] = cheby1(3,1,0.5,'low');
[H, W]= freqz(b,a);
subplot(2,2,1);
plot(W/pi, abs(H));
title('Low Pass Chebyshev filter - Order:3');
%*****************************************************************%  
% Low Pass Chebyshev with varying order = 5
%*****************************************************************%  
[b,a] = cheby1(5,1,0.5,'low');
[H, W]= freqz(b,a);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('Low Pass Chebyshev filter - Order:5');
%*****************************************************************%  
% Low Pass Chebyshev with varying order = 10
%*****************************************************************%  
[b,a] = cheby1(10,1,0.5,'low');
[H, W]= freqz(b,a);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Low Pass Chebyshev filter - Order:10');
%*****************************************************************%  
% Low Pass Chebyshev with varying order = 20
%*****************************************************************%  
[b,a] = cheby1(20,1,0.5,'low');
[H, W]= freqz(b,a);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Low Pass Chebyshev filter - Order:20');

% After analysing the plots we can say that the by increasing the order of
% the filter ripples incresase but at the same time transitions become
% steeper which gives us better result So it depend whetehr we want
% steepnes in transition in our filter so we use filter which gives more
% transition and if we want to reduce riiples so we use another filter.







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                            % 1D FILTERING 
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.1 1D Filtering


%*****************************************************************%  
                %DIRAC FUNCTION INPUT SIGNAL
%*****************************************************************% 

k = 40 ;  % Given value of K
x_k = zeros(k,1);  % Initiallizing with zeros
x_k (20,1) = 1;    % Creating Dirac Function
figure(4);   
plot(x_k);  % Plotting Dirac Function
title('Dirac Function '); xlabel('d(k) '); ylabel('y(k)');

%*****************************************************************%  
%SETTING SCALING FACTROR
%*****************************************************************%  
scaling_Factor = 0.5; %scale of smoothing
Ts = 1;  %sampling freq
alpha = scaling_Factor * Ts ; %Alpha Factor Value 
a_signal = exp(-alpha) ;
figure(5);
stem (a_signal) ;
title('Expnential exp(-alpha) '); 

%*****************************************************************%  
%CALLING ANTI CAUSAL DERAVATIVE FUNCTION WITH DIRAC AS INPUT
%*****************************************************************%  
y_k = anti_causal_dervative (a_signal,alpha , x_k);
figure(6)
stem (y_k) ;
title('Anti - Causal Deravative ');

%*****************************************************************%  
%CALLING CAUSAL DERAVATIVE FUNCTION WITH DIRAC AS INPUT
%*****************************************************************%  
y2_k = causal_deravative (a_signal,alpha , x_k);
figure(7)
stem (y2_k) ;
title('Causal Deravative'); 

%*****************************************************************%  
%CALLING ANTI CAUSAL SMOOTHING FUNCTION WITH DIRAC AS INPUT
%*****************************************************************%  
yb_k = anticausal_smoothing (a_signal,alpha , x_k ,scaling_Factor );
figure(8)
stem (yb_k) ;
title('anti causal Smmothing '); 

%*****************************************************************%  
%CALLING CAUSAL DERAVATIVE FUNCTION WITH DIRAC AS INPUT
%*****************************************************************%  
y2b_k = causal_smoothing (a_signal,alpha , x_k , scaling_Factor );
figure(9)
stem (y2b_k) ;
title('Causal Smmothing'); 


%*****************************************************************%  
                    % STEP FUNCTION INPUT SIGNAL
%*****************************************************************% 

%Step input Signal
k = 40 ;  % Given value of K
x2_k = zeros(k,1);   % Initiallizing with zeros
for i = 10 : 30  % Creating Dirac Function
x2_k (i) = 1; % Plotting Dirac Function
end
figure(10)
stem(x2_k)
title('Step Function '); xlabel('u[n] '); ylabel('y[n]');
%*****************************************************************%  
%CALLING ANTI CAUSAL DERAVATIVE FUNCTION WITH STEP AS INPUT
%*****************************************************************%  
y3_k = anti_causal_dervative (a_signal,alpha , x2_k);
figure(11)
stem (y3_k) ;
title('Anti - Causal Deravative ');
%*****************************************************************%  
%CALLING CAUSAL DERAVATIVE FUNCTION WITH STEP AS INPUT
%*****************************************************************%  
y4_k = causal_deravative (a_signal,alpha , x2_k);
figure(12)
stem (y4_k) ;
title('Causal Deravative ');
%*****************************************************************%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % EXERCISE : 3 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Canny-Deriche filtering
%Loading Image
%*****************************************************************%  
lena_image = imread('.\barbara.gif');
figure(13);
imshow(lena_image);  %Displaying Image
outImage1 = zeros(size(lena_image)); % Initiallizing Image with zeros
outImage2 = zeros(size(lena_image)); % Initiallizing Image with zeros

%*****************************************************************%  
% 3.2 : Applying the smoothing (derivative) filter along the columns (rows) of the images
%*****************************************************************%  
for i = 1:size(lena_image, 2)

    image_signal = lena_image(:,i);
    response_causal = causal_deravative (a_signal,alpha , image_signal  );
    response_antiCausal = anti_causal_dervative (a_signal,alpha , image_signal);
    response = response_causal + response_antiCausal;
    
    outImage1(:,i) = response;            
end
figure(14);
imshow (outImage1, []);

%*****************************************************************%  
% 3.3 : Applying the smoothing (derivative) filter along the rows (columns) of the images
%*****************************************************************%  
for i = 1:size(lena_image, 2)

    image_signal2 = lena_image(i,:);
    response_causal2 = causal_deravative (a_signal,alpha , image_signal2  );
    response_antiCausal2 = anti_causal_dervative (a_signal,alpha , image_signal2);
    response2 = response_causal2 + response_antiCausal2;
    
    outImage2(i,:) = response2;   
end

figure(15);
imshow (outImage2, []);
%*****************************************************************%  
                        % END 
%*****************************************************************%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
