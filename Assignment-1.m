
% Sourena Morteza Ghasemi(40171622)
% Amirreza Azadnia (40198570)


clc;
clear;
%EKF Assignment 1


%constants:
Q = [0.0025, 0, 0 ; 0 , 0.0025, 0; 0, 0, 0.025];
alpha_2 = 0.1;
alpha_1 = 0.1;
alpha_4 = 0.1;
alpha_3 = 0.1;
C=eye(3);  % In This Assignment we have three out put:X,Y,teta
dt=0.01;   % our sampling time
t=0;       %initial time
sikma_bar_prev=Q; %first guess of sikma
q=0;
m=0;
x_prev=0;      %from sensor
y_prev=0;      %from sensor
teta_prev=0;   %from sensor
v_t=0;         %from controller
w_t1=0.5;         %from controller  
H_t=C;       %becaus of the Assignemnt!

time=150;   %simulation time



%sensorTime;  %sensor time
%commandTime; %command Time
% Define file names
sensorFileName = 'Sensor(x,y,Theta)gravel_raw1.data';
commandFileName = 'Command(V,W)gravel_raw1.data';
commandData = dlmread(commandFileName);
commandTime = commandData(:,1);
sensorData = dlmread(sensorFileName);
sensorTime = sensorData(:,1);

initialsensortime=sensorTime(1,1);
initalcommandtime=commandTime(1,1);
% Read 'Sensor(x,y,Theta).data'

Xmatrix = sensorData(:,2);
Ymatrix = sensorData(:,3);
tetamatrix = sensorData(:,4);

% Read 'Command(V,W).data'

v_tmatrix = commandData(:,2);
w_tmatrix = commandData(:,3);


sensorTime1=initialsensortime;
commandTime1=initalcommandtime;





MuBarXArray = [];
MuBarYArray = [];
MuBarThetaArray = [];
x = [];
y = [];
teta = [];
V= [];
W= [];
timespan=[];
DET=[];

while t<= time   %simulation time
t=t+dt;
    
    
    if (t <= sensorTime1) && (t <= commandTime1)
    w_t=w_t1*v_t;
    f1=x_prev-(v_t/w_t)*sin(teta_prev)+(v_t/w_t)*sin(teta_prev+w_t*dt);
    f2=y_prev+(v_t/w_t)*cos(teta_prev)-(v_t/w_t)*cos(teta_prev+w_t*dt);
    f3=teta_prev+w_t*dt;
    f=[f1;f2;f3];
    
G=[1,0,-v_t*cos(teta_prev)/w_t+v_t*cos(teta_prev+w_t*dt)/w_t;0,1,-v_t*sin(teta_prev)/w_t+v_t*sin(teta_prev+w_t*dt)/w_t;0,0,1];




V_matrix=[ sin(teta_prev + w_t/100)/w_t - sin(teta_prev)/w_t  (v_t*sin(teta_prev))/w_t^2 + (v_t*cos(teta_prev + w_t/100))/(100*w_t) - (v_t*sin(teta_prev + w_t/100))/w_t^2;
 cos(teta_prev)/w_t - cos(teta_prev + w_t/100)/w_t  (v_t*cos(teta_prev + w_t/100))/w_t^2 - (v_t*cos(teta_prev))/w_t^2 + (v_t*sin(teta_prev + w_t/100))/(100*w_t);
                                                 0                                                                                                        1/100];
M_t=[alpha_3*v_t*v_t+alpha_4*w_t*w_t,0;0,alpha_2*v_t*v_t+alpha_1*w_t*w_t];

sikma_bar=G*sikma_bar_prev*G'+V_matrix*M_t*V_matrix';  

P =  sikma_bar + Q;

% Calculate the inverse

z=[x_prev;y_prev;teta_prev];


lambda = 1e-6;  % Small positive constant
Pinv = inv(P + lambda * eye(size(P)));
     % Pinv=inv(P);
      K= sikma_bar  *Pinv;
      
disp(M_t)


%//////end of prediction step
%Update start

MuBarT=f+K*(z-f);

 

%sikma_bar_prev=(eye(3)-K*H_t)*sikma_bar;
%update end


if det(sikma_bar) > 1e-30
    sikma_bar_prev = (eye(3) - K * eye(3)) * sikma_bar;
    DET(m+1)=det(sikma_bar_prev);
    else
   sikma_bar = eye(3);
end


 % 
     %disp(v_t)
    %disp(w_t)

elseif  t>= sensorTime1
 %update mesurments inputs
 %update mesurment time
m=m+1; 
x_prev = Xmatrix(m, 1);
y_prev = Ymatrix(m ,1);
teta_prev = tetamatrix(m, 1);
sensorTime1 = sensorTime(m ,1);

elseif  t>= commandTime1
%update control inputs
%update control time
q=q+1;     
v_t = v_tmatrix(q, 1);
w_t1 = w_tmatrix(q, 1);


commandTime1 = commandTime(q, 1);
    end

MuBarXArray(m+1) = MuBarT(1);
     MuBarYArray(m+1) = MuBarT(2);
      MuBarThetaArray(m+1) = MuBarT(3);

   x(m+1)=x_prev;
   y(m+1)=y_prev;
   teta(m+1)=teta_prev;
   V(m+1)= v_t;
   W(m+1) = w_t;
   timespan(m+1)=sensorTime1;

   
end


lowerLimit = 1.5708;  % Lower limit
upperLimit = 7.8540;  % Upper limit
MuBarThetaArray(MuBarThetaArray < lowerLimit) = MuBarThetaArray(MuBarThetaArray < lowerLimit) + 2*pi;
MuBarThetaArray(MuBarThetaArray > upperLimit) = MuBarThetaArray(MuBarThetaArray > upperLimit) - 2*pi;
 

teta_prev(teta_prev < lowerLimit) = teta_prev(teta_prev < lowerLimit) + 2*pi;
teta_prev(teta_prev > upperLimit) = teta_prev(teta_prev > upperLimit) - 2*pi;
 


figure;
plot(timespan, MuBarYArray, 'LineWidth', 1)
hold on
 plot(timespan, y,'r');
 title('mu_y against t and y against t');
 xlabel(' t');
 ylabel(' mu_y, y');
 % Create a legend
legend('mu_y/t', 'y/t');
grid on;    
    
% Plot mu_x against t and x against t
figure;
plot(timespan, MuBarXArray, 'LineWidth', 1)
hold on
 plot(timespan, x,'r');
 title('mu_x against t and x against t');
 xlabel(' t');
 ylabel(' mu_x, x');
 % Create a legend
legend('mu_x/t', 'x/t');
grid on;    
  
% Plot mu_x against mu_y and x against y
figure;
plot(MuBarXArray, MuBarYArray, 'LineWidth', 1)
hold on
 plot(x, y,'r');
 title('mu_x against mu_y and x against y');
 xlabel('mu_y, y');
 ylabel(' mu_x, x');
 % Create a legend
legend('mu_x/mu_y', 'x/y');
grid on;

% Plot mu_tetha against t and tetha against t
figure;
plot(timespan, MuBarThetaArray, 'LineWidth', 1)
hold on
 plot(timespan, teta ,'r');
 title('mu_tetha against t and tetha against t');
 xlabel('t');
 ylabel('mu_tetha, tetha');
xlim([70 80])
 % Create a legend
legend('mu_tetha/t', 'tetha/t');
grid on;

% Plot v_t against t 
figure;
plot(timespan, V, 'LineWidth', 1)
 title('v_t against t');
 xlabel('t');
 ylabel('v_t');
 % Create a legend
grid on;

% Plot v_t against t 
figure;
plot(timespan, W, 'LineWidth', 1)
 title('w_t against t');
 xlabel('t');
 ylabel('w_t');
 % Create a legend
grid on;


figure;
plot(timespan, DET, 'LineWidth', 1)
 title('det_sigma against t');
xlabel('t');
 ylabel('det_Sigma');
% xlim([50 55])
 % Create a legend
grid on;