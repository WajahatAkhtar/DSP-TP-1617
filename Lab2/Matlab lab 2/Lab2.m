%Lab 2 Digital Signal Processing
function y1 =  Lab2(f,fs) %Function Defination 


 %*****************************************************************%   
%Part 1 Reminder
         t1 = [0:1/fs:10]; %
         y1 = sin(2*pi*f*t1);
         figure(1)
         subplot(1,2,1),plot(t1,y1);
         
         t2 = [0:1/fs:20];
         y2 = sin(2*pi*(f/fs)*t2);
         subplot(1,2,2),plot(t2,y2);
        
%*****************************************************************%

%*****************************************************************%
        x = step(4,20);
         for k= 1:19
            y(k) = x(k)/2+ (x(k+1))/2;%Making signal non causal because k = k+1 it contains future values
         end
        figure(4)
        stem(y)
  
        %Making signal causal by shifting it to past by putting k = k-1 and
        %changing the index becuase matlab doesnt consider 0
        for k= 2:20
            y(k) = x(k)/2+ (x(k-1))/2;
        end
        figure(5)
        stem(y)
%*****************************************************************%          
        %Stability  
        %%% Part 1
        y5 = x; %x is a step 
        for i = 2:1:20
            y5(i) = y5(i-1)+x(i);
        end
        figure(6)
        stem(y5)
         %%Comments :
        
        
        %%%Part 2
         d = Dirac(4,20);
         y6 = d;
         for i = 2:1:20
            y6(i) = y6(i-1)+d(i);
        end
        figure(7)
        stem(y6)
         %%Comments :
         
        %%%Part 2
         d = Dirac(4,20);
         y7 = d;
         for i = 2:1:20
            y7(i) = d(i)+ 2*(y7(i-1));
        end
        figure(8)
        stem(y7)
         %%Comments :
        
        %%%Part 2
         d = Dirac(4,20);
         y8 = d;
         for i = 2:1:20
            y8(i) = d(i)+ (y8(i-1)/3);
        end
        figure(9)
        stem(y8)
        %%Comments : same impulse response with opposite direction of
        %%impulse
        
 %*****************************************************************%          
         
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%




   
  
 