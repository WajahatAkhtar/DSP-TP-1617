function Lab1()  %Displaying all the functions
   
%%%% Question : 1.1 
    y1 = Dirac(0,20);  %Dirac Function 
    figure(1); stem(y1);
    title('Dirac Function')  %Plotting Dirac Function 
    xlabel('n')
    ylabel('delta(k-n)')

%%%% Question : 1.2   
    y2 = step(10,20);  %Dirac Function with Step
    figure(2); stem(y2);
    title('Dirac Function with step n = 10')   %Plotting Dirac Function with Step
    xlabel('n')
    ylabel('step(k-n)')
    
%%%% Question : 1.3
    y3 = ramp(2,10,20);  %Ramp Function 
    figure(3); stem(y3);
    title('Ramp Function with a=2 and n=10')   %Plotting Ramp Function 
    xlabel('n')
    ylabel('P(k-n)')
    
%%%% Question : 1.4
    y4 = geo(2,10,20);  % Geometric Function
    figure(4); stem(y4);
    title('Geometric Function')   %Plotting Geometric Function
    xlabel('n')
    ylabel('G(k-n)')
    
%%%% Question : 1.5
    y5 = box(3,10,20);  % BOX Function
    figure(5); stem(y5);
    title('BOX Function')   %Plotting BOX Function 
    xlabel('n')
    ylabel('B(k-n)')
    
%%%% Question : 1.6
    y6_a = sinfunc(1,10,100);  % Sin Function
    figure(6); stem(y6_a);
    title('Sin Function with period 1 of frequency 10hz and sampling Frequncy 100')   %Plotting Sin Function 
    xlabel('n')
    ylabel('sin(2*pi*f*t)')
    
    y6_b = sinfunc(2,10,1000);  % Sin Function
    figure(7); stem(y6_b);
    title('Sin Function with period 2 of frequency 10hz and sampling Frequncy 1000')   %Plotting Sin Function
    xlabel('n')
    ylabel('sin(2*pi*f*t)')
    
    y6_c = sinfunc(2,10,30);  % Sin Function
    figure(8); stem(y6_c);
    title('Sin Function with period 2 of frequency 10hz and sampling Frequncy 30')   %Plotting Sin Function
    xlabel('n')
    ylabel('sin(2*pi*f*t)')
    
    
    
%%%%%%%%************ Exercice : 2  Random signals *****************%%%%%%%%


%%% Question 2.1 : Generate an observation xn 1000 of the normal/gaussian random process N.
%%%%%%----For xn=1000---%%%%%%
 rand1 = randn(1,1000);  %%Using randn matlab function for generating arrays of random numbers whose elements are normally distributed with mean 0, variance , and standard deviation
 standard_deviation = 1; 
 mean = 0;
 figure(9);  hist(rand1,100) %%Plotting Histogram of random signal
 title('Histogram with xn=1000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 Theoritical_dist = normpdf(rand1,mean,standard_deviation);
 figure(10);  stem(rand1,Theoritical_dist);%% Plotting Theoritical distribution of the signal
 title('Theoretical distribution with xn=1000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
%%%%%%----For xn=10000---%%%%%%
 rand2 = randn(1,10000);  %%Using randn matlab function for generating arrays of random numbers whose elements are normally distributed with mean 0, variance , and standard deviation
 figure(11);  hist(rand2,100)
 title('Histogram with xn=10000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 Theoritical_dist = normpdf(rand2,mean,standard_deviation);
 figure(12); stem(rand2,Theoritical_dist);
 title('Theoretical distribution with xn=10000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 %%Comments on above question
 % As Gaussian noise is statistical noise having a probability density 
 % function (PDF) equal to that of the normal distribution,So after
 % inceasing the number of samples to 10000 we can see a clear change in
 % the fgures the we have now more samples in the same region with more
 % implicitly means more information with less noise. 
 

 
 
 %%% Question 2.2 : Generate an observation xu 1000 of the uniform law of the random process U and an observation xu.
 %%%%%%----For xu=1000---%%%%%%
 rand3=rand(1,1000);
 figure(13); hist(rand3,100) %%Plotting Histogram of random signal
 title('Histogram with xu=1000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 Theoritical_dist = unifpdf(rand3,mean,standard_deviation);
 figure(14); stem(rand3,Theoritical_dist);
 title('Theoretical distribution with xu=1000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 %%%%%%----For xu=10000---%%%%%%
 rand4=rand(1,10000);
 figure(15);  hist(rand4,100) %%Plotting Histogram of random signal
 title('Histogram with xu=10000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 Theoritical_dist = unifpdf(rand4,mean,standard_deviation);
 figure(16); stem(rand4,Theoritical_dist);
 title('Theoretical distribution with xu=10000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 %%Comments on above question
 %There is only change in there intensities with uniform distribution
 
 

 
%%% Question 2.3 :Compute the autocorrelation of the two observations and plot them. 
 [r1,lags1] = xcorr(rand1); %Cross-correlation
 figure(17);
 plot(lags1,r1)
 title('Autocorrelation xn=1000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 %%Comments on above question
  %We can see in the figure 17 that it has a peak in the middle which tells
  %that it is a stochastic process with more variation on the sides.

 [r2,lags2] = xcorr(rand2);
 figure(18);
 plot(lags2,r2)
 title('Autocorrelation xn=10000 normal/gaussian random process')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
  %%Comments on above question
  %by increasing number of sample to 10000 We can see the same result with
  %less variation in the sides.
 
 [r3,lags3] = xcorr(rand3);
 figure(19);
 plot(lags3,r3)
 title('Autocorrelation xu=1000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
  %We can see in the figure 19 that it has huge variation in the middle not 
  %a sigle peak as it was in case of random signal above which tells it is 
  %a not stochastic process/white noise with more variation on the sides.
 
 [r4,lags4] = xcorr(rand4);
 figure(20);
 plot(lags4,r4)
 title('Autocorrelation xu=10000 uniform law of the random process U')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
  %by increasing number of sample to 10000 We can see the same result with
  %less variation in the sides.
 
  
  
 
 %%% Question 2.4 :Generate three binary random signals s1; s2; s3 
  s1 = round(rand(1,50));  % Generating binary random signal 1
  s2 = round(rand(1,50));  % Generating binary random signal 2
  s3 = round(rand(1,50));  % Generating binary random signal 3
  s=[zeros(1,15) s1 zeros(1,15) s2 zeros(1,15) s3 zeros(1,15)]; % Generating Whole signal with shift of 15
 
 figure(21);
 plot(s); %Plotting whole signal
 title('Whole Signal with shift of 15')
 xlabel('Array of Random numbers')
 ylabel('Intensity at each index)')
 
 figure(22);
 [r1,lags1] = xcorr(s,s1); %cross correlation of whole signal and s1
 plot(lags1,r1) %Plotting cross correlation of s1
 title('Cross-correlation of binary random signal s1 with the Whole Signal')
 xlabel('Binary random signal 1')
 ylabel('Intensity at each index)')
 
 figure(23);
 [r1,lags1] = xcorr(s,s2); %cross correlation of whole signal and s2
 plot(lags1,r1) %Plotting cross correlation of s2
 title('Cross-correlation of  binary random signal s2 with the Whole Signal')
 xlabel('Binary random signal 2')
 ylabel('Intensity at each index)')
 
 
 figure(24);
 [r1,lags1] = xcorr(s,s3); %cross correlation of whole signal and s3
 plot(lags1,r1) %Plotting cross correlation of s3
 title('Cross-correlation of  binary random signal s3 with the Whole Signal')
 xlabel('Binary random signal 3')
 ylabel('Intensity at each index)')
 
 %%Comments
 %We can clearly see the peaks for each binary random signal when we do
 %cross-correlation of that signal which provides us a method to find
 %signal patterns.
 
end
%%-----------------------------------------------------------------------%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------Defination of all the functions called above----------%%%%%%%%%%

%%%********Exercice : 1 Deterministic signals*********************%%%%%%%%
%Question 1.1 : Dirac Function 
function y1 =  Dirac(n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1');  %Display error if n > N-1
    else
            s = zeros(N,1);  
            s(n+1) = 1 ;
            y1 = s;
           
    end
  
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%


%Question 1.2 : Dirac Function with step 10
function y2 =  step(n,N)  %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1');  %Display error if n > N-1
    else
            s = zeros(N,1);  
            for i = n+1:N
                s(i) = 1 ;
            end 
            y2 = s;
           
    end
  
end 
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%
    


%Question 1.3 : Ramp Function
function y3 =  ramp(a,n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1');  %Display error in n > N-1
    else
            
            s = zeros(1,N);   
            for j = n:N
                s(j) = a*(j-n);
            end
            y3 = s;
                  
    end
  
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%



%Question 1.4 :  Geometric function
function y4 =  geo(a,n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1'); %Display error in n > N-1
    else
            
            s = zeros(1,N);  
            for j = n:N
                s(j) = a^(j-n);
            end
            y4 = s;
                  
    end
  
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%


%Question 1.5 :  BOX function
function y5 =  box(a,n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1'); %Display error in n > N-1
    else
            
            s = zeros(1,N);  
            for j = (n-a):(n+a)
                s(j) = 1 ;
            end
            y5 = s;
                  
    end
  
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%




%Question 1.6 :  Sin function
function y6 =  sinfunc(n,f,Ts) %Function Defination 

         t = [0:1/Ts:n/f];
         y6 = sin(2*pi*f*t);
     
           
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%




















