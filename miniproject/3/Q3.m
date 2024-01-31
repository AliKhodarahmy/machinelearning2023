clc
clear all
close all
%% GRADIEN DECENT TRAINING ALGORITHEM.
disp(' GRADIEN DECENT TRAINING ALGORITHEM.');
disp('  This programm develops as a course project for Fuzzy systems course')
disp('   by Habibollah Naeimi.')
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 1st Part: Parameter Setting.
%%  Parameters Initiating.
disp(' Parameters Initiating...');

M = 30;            
Alpha = 0.3;
Q = 100;                               
epsilon = 0;
InpuNum = 1;                       
y_Bar = zeros(M,1);                  
x_Bar = zeros(M,InpuNum);
Sigma = zeros(M,InpuNum);

disp(' Part 1: DONE!');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');

%% 2rd Part: Sampling.
%%  Calculating Samples.
Data_Pairs_Num = 300;      

SAMPLES_Num = 500;               
a1 = 0.3;
a2 = 0.6;

y_aprx = zeros(1,SAMPLES_Num);       
SAMPLES = zeros(SAMPLES_Num,InpuNum+1);
iNSamples = 0.2:0.01:0.21;
for k=2:SAMPLES_Num+33+InpuNum
    iNSamples(k + 1) = a1*iNSamples(k) + a2*iNSamples(k - 1) + 0.6*sin(pi*k) + 0.3*sin(3*pi*k) + 0.1*sin(5*pi*k);
end

iNSamples = iNSamples(2:end);
    
for i=1:SAMPLES_Num
    SAMPLES(i,:) = iNSamples(i:i+InpuNum);
end

Pairs = SAMPLES(1:Data_Pairs_Num,:);
disp(' Complete Sampling.');
%%  Primary Parameter Fixing Using Online initial Parameter Choosing.

x_Bar = Pairs(1:M,1:InpuNum);         
y_Bar = Pairs(1:M,end);
Sigma = repmat(((max(x_Bar)-min(x_Bar))/M),M,1);

disp(' Initial Parameters are Reasdy!');
disp(' ');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 3rd Part: Parameters Updating. 
%%  Updating.

z = zeros(M,1);       
iN_z = zeros(1,InpuNum);

 for p=1:size(Pairs,1)
   
     for l=1:M     
         for i=1:InpuNum
             iN_z(i) = exp(-(((Pairs(p,i)-x_Bar(l,i))/Sigma(l,i))^2));
         end
         z(l) = prod(iN_z);
     end
            
     b = sum(z);         
     a = sum(y_Bar.*z);    
     f = a/b;          

            for q=1:Q
                for l=1:M        
                    y_Bar(l) = y_Bar(l)-Alpha*(f-Pairs(p,end))/b*z(l);
                    for i=1:InpuNum                    
                        x_Bar(l,i) = x_Bar(l,i)-Alpha*(f-Pairs(p,end))/b*(y_Bar(l)-f)*z(l)*(2*(Pairs(p,i)-x_Bar(l,i))/(Sigma(l,i)^2));
                        Sigma(l,i) = Sigma(l,i)-Alpha*(f-Pairs(p,end))/b*(y_Bar(l)-f)*z(l)*(2*((Pairs(p,i)-x_Bar(l,i))^2)/(Sigma(l,i)^3));
                    end
                end
                
                if (f-Pairs(p,end))<epsilon 
                    break;
                end
            end
 end

disp(' Parameters Updated.');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 4th Part: Results.
%%  Results.
 
f = zeros(1,SAMPLES_Num);          
f(1:2) = SAMPLES(1:2,end);

for k=3:SAMPLES_Num
    
    for l=1:M                   
        for i=1:InpuNum
            iN_z(i) = exp(-(((SAMPLES(k,i)-x_Bar(l,i))/Sigma(l,i))^2));                    
        end
        z(l) = prod(iN_z);
    end
           
    b = sum(z);            
    a = sum(y_Bar.*z);               
    f(k) = a/b;                   
    y_aprx(k) = f(k);                  
             
end
%%  Plotting.

figure;
plot(SAMPLES(:,end));
hold on
plot(y_aprx,'r');
legend('Real Value','Predicted Value');

%%  Errors.
Error = SAMPLES(:,end)-y_aprx';
disp('Mean Square Error is:');
MSE = mse(Error)
disp('Mean Absolute Error is:');
MAE = mae(Error)

