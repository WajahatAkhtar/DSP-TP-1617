function Lab1()   %Displaying all the functions
    y1 = Dirac(0,20); %Dirac Function 
    figure(1); stem(y1);
    title('Dirac Function')  % Plotting Dirac Function 
    
    
    y2 = step(10,20); %Dirac Function with Step
    figure(2); stem(y2);
    title('Dirac Function with step n = 10')   % Plotting Dirac Function with Step
    
    
    
    y3 = ramp(2,10,20); %Dirac Function with Step
    figure(3); stem(y3);
    title('Ramp Function from a=2 and n=10')   % Plotting Dirac Function with Step
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 1 : Plotting Dirac Function 
function y1 =  Dirac(n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1'); %Display error in n > N-1
    else
            s = zeros(N,1);  % 
            s(n+1) = 1 ;
            y1 = s;
           
    end
  
end
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 2 :  Plotting Dirac Function with step 10
function y2 =  step(n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1'); %Display error in n > N-1
    else
            s = zeros(N,1);  % 
            for i = n+1:N
            s(i) = 1 ;
            end
            y2 = s;
           
    end
  
end
 
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%
    


%Question 2 :  Plotting Dirac Function with step 10
function y3 =  ramp(a,n,N) %Function Defination 

    if n >(N-1)
            Disp('Error : n should be inferior then N-1'); %Display error in n > N-1
    else
            
            s = zeros(1,N);  % 
            
            for j = n+1:N
            s(j) = 1;
            end
            
            for j = n+2:N
            s(j) = s(j-1)+s(j);
            end
            
            
            y3 = s;
            
           
    end
  
end
 
%%%%%%%%%----------------------------------------------------%%%%%%%%%%%%%%




