% This code allows for the prediction of UAV states during a blackout 

g= 9.81;

load ('axControlledWindStep.mat')
load ('ayControlledWindStep.mat')
load ('azControlledWindStep.mat')
load ('pControlledWindStep.mat')
load ('qControlledWindStep.mat')
load ('rControlledWindStep.mat')
load('Float.mat')

load ('xControlledWindStep.mat')
load ('yControlledWindStep.mat')
load ('zControlledWindStep.mat')
load ('uControlledWindStep.mat')
load ('vControlledWindStep.mat')
load ('wControlledWindStep.mat')
load ('rollControlledWindStep.mat')
load ('pitchControlledWindStep.mat')
load ('yawControlledWindStep.mat')

%Observation vector z =[V , beta, alpha, xe,ye,h]
load ('VelocityControlledWindStep.mat')
load ('alphaControlledWindStep.mat')
load ('betaControlledWindStep.mat')


load ('Da.mat')
load ('De.mat')
load ('Dr.mat')
load ('Dt.mat')

load ('xeBefore.mat')
load ('yeBefore.mat')
load ('zeBefore.mat')
load ('uBefore.mat')

xeBef=  xeBefore(2,:);
yeBef=  yeBefore(2,:);
zeBef=  zeBefore(2,:);
uBef=  uBefore(2,:);

de_veri = De(2,:);
da_veri = Da(2,:);
dr_veri = Dr(2,:);
dT_veri = Dt(2,:);
Float_veri = Float(2,:);

Runs = 50; % Monte Carlo Runs
MSExUKF= zeros(1,Runs);
MSEyUKF= zeros(1,Runs);
MSEzUKF= zeros(1,Runs);

MSEuUKF= zeros(1,Runs);
MSEvUKF= zeros(1,Runs);
MSEwUKF= zeros(1,Runs);
 
MSErollUKF= zeros(1,Runs);
MSEpitchUKF= zeros(1,Runs);
MSEyawUKF= zeros(1,Runs);

%

MSExCKF= zeros(1,Runs);
MSEyCKF= zeros(1,Runs);
MSEzCKF= zeros(1,Runs);

MSEuCKF= zeros(1,Runs);
MSEvCKF= zeros(1,Runs);
MSEwCKF= zeros(1,Runs);

MSErollCKF= zeros(1,Runs);
MSEpitchCKF= zeros(1,Runs);
MSEyawCKF= zeros(1,Runs);

ax_meas = axControlledWindStep(2,:);
ay_meas = ayControlledWindStep(2,:) ;
az_meas = azControlledWindStep(2,:) ;
p_meas = pControlledWindStep(2,:) ;
q_meas = qControlledWindStep(2,:) ; 
r_meas = rControlledWindStep(2,:) ;

%Measured data, u

Velocity = VelocityControlledWindStep;
alpha = alphaControlledWindStep  ;
beta = betaControlledWindStep ;
x = xControlledWindStep;
y = yControlledWindStep ;
z = zControlledWindStep ;
u = uControlledWindStep;
v = vControlledWindStep;
w = wControlledWindStep;
roll = rollControlledWindStep;
pitch = pitchControlledWindStep;
yaw = yawControlledWindStep;


z_obse = [Velocity(2,:) ; alpha(2,:) ; beta(2,:) ; x(2,:) ; y(2,:) ; z(2,:)];

%Trim conditions which are x= [u; v; w; Roll_A; Pitch_A; Yaw_A; x; y; z]
x_true = [u(2,:) ; v(2,:) ; w(2,:) ; roll(2,:) ; pitch(2,:) ; yaw(2,:) ; xControlledWindStep(2,:) ; yControlledWindStep(2,:) ; zControlledWindStep(2,:)]; %No sensor Noise
x_sensor = [u(2,:) ; v(2,:) ; w(2,:) ; roll(2,:) ; pitch(2,:) ; yaw(2,:) ; x(2,:) ; y(2,:) ; z(2,:)]; %With Sensor noise
xo = [u(2,1) ; v(2,1) ; w(2,1) ; roll(2,1) ; pitch(2,1) ; yaw(2,1) ; x(2,1) ; y(2,1) ; z(2,1)];

for n = 1:1:Runs %Monte Carlo
    
time = axControlledWindStep(1,:);
N = length(time);
XXX = Velocity(2,:);
%random errors
W1yEr = (0.65e-4).*rand(1,1);
W2yEr = (3.13e-6).*rand(1,1);
W3yEr = (9e-4).*rand(1,1) ;

W4yEr= (0.05).*rand(1,1) ;
W5yEr= (0.05).*rand(1,1);
W6yEr= (0.02).*rand(1,1);

V1Er= (0.0158).*rand(1,1);
V2Er= (1e-4).*rand(1,1);
V3Er= (2.5e-2).*rand(1,1);
V4yEr= (0.0000001).*rand(1,1);
V5Er= (0.0000004).*rand(1,1);
V6Er= (0.2).*rand(1,1) ;
V7Er= (2).*rand(1,1);
V8Er= (1).*rand(1,1) ;
V9Er= (2).*rand(1,1) ;

%xo_dot = [udotWind(2,1) ; vdotWind(2,1) ; wdotWind(2,1) ; pdotWind(2,1) ; qdotWind(2,1) ; rdotWind(2,1) ;  xdotWind(2,1) ; ydotWind(2,1) ; zdotWind(2,1)];
xo_dot = [-1.1701e-8 ; 6.256e-20 ; 9.5279e-12 ; 5.0037e-19 ; -1.44995e-20 ; -1.2725e-19 ; 45.0741 ; 7.5205e-17 ; -3.3172e-9];

errors = [0.01369 ; 0.01369 ; 0.01369 ;2e-9 ; 1e-4 ;2.6e-1 ; 2 ; 1 ; 2] ;  %Remember p,q,r is in rad
P = diag(errors); 
yerrors = [ W1yEr ;W2yEr ;W3yEr ;W4yEr ;W5yEr ;W6yEr];
xerrors =[V1Er; V2Er; V3Er ;V4yEr ;V5Er ;V6Er ; V7Er ; V8Er ; V9Er];

R = diag (yerrors);
Q = diag (xerrors);
%Scaling Parameters
alpha = 1;          % Can be between 1 and 1e-14
beta = 2;           % Optimum for Gaussian distribution
k =  0;             % Because L>3 
L = length(xo);
Lamda = (alpha^2)*(L+ k) - L;
eta = sqrt(L+Lamda) ; 

% Weight vectors for steps 0 and 1
neu_m = zeros(1,(2*L + 1));
neu_c = zeros(1,(2*L + 1));
neu_o_m = Lamda/(L + Lamda);
neu_o_c= (Lamda/(L + Lamda)) + 1 - alpha^2 + beta;
neu_m (1) = neu_o_m;
neu_c (1) = neu_o_c;

% Weight vectors for steps 2 onwards
for i = 2: (2*L+1)
neu_m (i) = 1/(2*(L + Lamda));
neu_c (i) = neu_m (i);
end
neu_m = neu_m';
neu_c= neu_c';

%initialising estimate state vector for all u length
x = zeros (L,N) ; 
x (:,1) = xo ; 
u_dot = zeros(1,L);
v_dot = zeros(1,L);
w_dot = zeros(1,L);
roll_dot = zeros(1,L);
pitch_dot = zeros(1,L);
yaw_dot = zeros(1,L);
x_dot = zeros(1,L);
y_dot = zeros(1,L);
h_dot = zeros(1,L);
V = zeros(1,L);
alpha = zeros(1,L);
beta =zeros(1,L);
% X_int =zeros(L,2*L+1);
x_dot_kminus1 = [xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot xo_dot]; 
X_int= [xo xo xo xo xo xo xo xo xo xo xo xo xo xo xo xo xo xo xo];

for k= 2:N
   
if XXX(k) == 0  
    z_obse(:,k) = y_hat_k_kminus1;
    %x_true(1,k) = NaN;
    x_true(7,k) = NaN;
    x_true(8,k) = NaN ;
    x_true(9,k) = NaN;
end
    
%Determining sqrt of P
P_k_kminus_1 = P;
sqrtP = chol(P_k_kminus_1,'lower');              %can use sqrtm if this doesnt work

%Generating Sigma Points

X_Sigma = [x(:,k-1) , (x(:,k-1)+ (eta*sqrtP)) , (x(:,k-1) - (eta*sqrtP))];

u_state = X_Sigma (1,:);
v_state = X_Sigma (2,:);
w_state = X_Sigma (3,:);
Roll_state = X_Sigma (4,:);
Pitch_state = X_Sigma (5,:);
Yaw_state = X_Sigma (6,:);

%Prediction Transformation, will calculate x_dot  

%translation kinematics
for i = 1: 2*L + 1
u_dot (i) = (r_meas(k-1)*v_state (i)) - (q_meas(k-1)*w_state(i)) - (g*sin(Pitch_state(i)))+ (ax_meas(k-1));
v_dot (i) = (p_meas(k-1)*w_state (i)) - (r_meas(k-1)*u_state(i)) + (g*cos(Pitch_state(i))*sin(Roll_state(i)))+ (ay_meas(k-1));
w_dot (i)  = (q_meas(k-1)*u_state (i)) - (p_meas(k-1)*v_state(i)) + (g*cos(Pitch_state(i))*cos(Roll_state(i)))+ (az_meas(k-1));

%rotational kinematics
roll_dot (i)  = p_meas(k-1) + ((q_meas(k-1)*sin(Roll_state(i))*tan(Pitch_state(i))) + (r_meas(k-1)*cos(Roll_state(i))*tan(Pitch_state(i))));
pitch_dot (i)  = ((q_meas(k-1)*cos(Roll_state(i))) - (r_meas(k-1)*sin(Roll_state(i))));
yaw_dot(i)  = ((q_meas(k-1)*sin(Roll_state(i))) + (r_meas(k-1)*cos(Roll_state(i))))/ cos(Pitch_state(i));


%Position kinematics
x_dot(i)  = u_state(i)*cos(Pitch_state(i))*cos(Yaw_state(i))  + v_state(i)*((sin(Roll_state(i))*sin(Pitch_state(i))*cos(Yaw_state(i)))- (cos(Roll_state(i))*sin(Yaw_state(i)))) + w_state(i)*((cos(Roll_state(i))*sin(Pitch_state(i))*cos(Yaw_state(i))) +  (sin(Roll_state(i))*sin(Yaw_state(i)))) ;

y_dot (i) = u_state(i)*cos(Pitch_state(i))*sin(Yaw_state(i))  + v_state(i)*((sin(Roll_state(i))*sin(Pitch_state(i))*sin(Yaw_state(i))) + (cos(Roll_state(i))*cos(Yaw_state(i)))) + w_state(i)*((cos(Roll_state(i))*sin(Pitch_state(i))*sin(Yaw_state(i))) -  (sin(Roll_state(i))*cos(Yaw_state(i)))) ;

h_dot (i) = u_state(i)*sin(Pitch_state(i)) - v_state(i)*sin(Roll_state(i))*cos(Pitch_state(i)) - w_state(i)*cos(Roll_state(i))*cos(Pitch_state(i));

end
%Saving for integration

chi_k_kminus1_dot = [u_dot ; v_dot; w_dot; roll_dot; pitch_dot; yaw_dot; x_dot; y_dot; h_dot];

for j = 1:2*L+1   %goes from 1 to 19 ie through all columns
    X_int(:,j) = X_int(:,j) + 0.5*(time(k) - time(k-1))*(x_dot_kminus1(:,j) + chi_k_kminus1_dot(:,j)); 
end

chi_k_kminus1 = X_int;
x_dot_kminus1 = chi_k_kminus1_dot;
%Mean of predicted state
x_hat_sum = 0;
for i = 1: 2*L+1
    x_hat_sum  = x_hat_sum  + neu_m(i)*chi_k_kminus1(:,i);
end
x_hat_k_kminus1 =x_hat_sum; % chi_k_kminus1*neu_m;

%Getting the covariance of predicted state
P_k_minus_1 = Q;
for j = 1: (2*L+ 1)
P_k_minus_1 = P_k_minus_1 + neu_c(j)*(chi_k_kminus1(:,j)- x_hat_k_kminus1)*transpose(( chi_k_kminus1(:,j)- x_hat_k_kminus1));
end

u_state = chi_k_kminus1 (1,:);
v_state = chi_k_kminus1 (2,:);
w_state = chi_k_kminus1 (3,:);
Roll_state = chi_k_kminus1 (4,:);
Pitch_state = chi_k_kminus1 (5,:);
Yaw_state = chi_k_kminus1 (6,:);
x_state = chi_k_kminus1 (7,:);
y_state = chi_k_kminus1 (8,:);
z_state = chi_k_kminus1 (9,:); 

%Observation Update y =[V , beta, alpha, xe,ye,h]
for w = 1: 2*L + 1
V (w)= sqrt(u_state(w)^2 + v_state(w)^2 + w_state(w)^2) ; 
alpha (w) = atan (w_state(w)/u_state(w));
beta (w) = asin(v_state(w)/(sqrt(u_state(w)^2 + v_state(w)^2 + w_state(w)^2)));
end

Psi_k_kminus_1 = [V; alpha; beta; x_state; y_state; z_state];
Ly = size(Psi_k_kminus_1);
%The predicted output
y_hat_k_kminus1 = Psi_k_kminus_1*neu_m ; 

P_k_yy = R;
P_k_xy = zeros(L,Ly(1,1));

for j = 1: 2*L +1
P_k_yy = P_k_yy + neu_c(j)*( Psi_k_kminus_1(:,j)- y_hat_k_kminus1)*(( Psi_k_kminus_1(:,j)- y_hat_k_kminus1))';
P_k_xy = P_k_xy +  neu_c(j)*( chi_k_kminus1(:,j)- x_hat_k_kminus1)*(( Psi_k_kminus_1(:,j)- y_hat_k_kminus1))';
end

Kk = P_k_xy/P_k_yy;
%Measurment update
x(:,k) = x_hat_k_kminus1 + Kk*(z_obse(:,k) - y_hat_k_kminus1);
P = P_k_minus_1 - Kk*P_k_yy*Kk' ;
end
x_UKF = x;
x_true(7,:); x_true(8,:); x_true(9,:);


MSExUKF(n)=(sum( (1/N)*(xeBef - x(7,1:k)).^2))^0.5;
MSEyUKF(n)= (sum((1/N)*(yeBef - x(8,1:k)).^2))^0.5;
MSEzUKF(n)= (sum((1/N)*(zeBef - x(9,1:k)).^2))^0.5;

MSEuUKF(n)= (sum((1/N)*(uBef - x(1,1:k)).^2))^0.5;
MSEvUKF(n)= (sum((1/N)*(x_true(2,:) - x(2,1:k)).^2))^0.5;
MSEwUKF(n)= (sum((1/N)*(x_true(3,:) - x(3,1:k)).^2))^0.5;

MSErollUKF(n)= (sum((1/N)*(x_true(4,:) - x(4,1:k)).^2))^0.5;
MSEpitchUKF(n)= (sum((1/N)*(x_true(5,:) - x(5,1:k)).^2))^0.5;
MSEyawUKF(n)= (sum((1/N)*(x_true(6,:) - x(6,1:k)).^2))^0.5;

%random errors
W1yEr = (0.65e-4).*rand(1,1);
W2yEr = (3.13e-6).*rand(1,1);
W3yEr = (9e-4).*rand(1,1) ;

W4yEr= (0.05).*rand(1,1) ;
W5yEr= (0.05).*rand(1,1);
W6yEr= (0.02).*rand(1,1);

V1Er= (1).*rand(1,1);
V2Er= (0.01369).*rand(1,1);
V3Er= (2.5e-2).*rand(1,1);
V4yEr= (2e-8).*rand(1,1);
V5Er= (1e-8).*rand(1,1);
V6Er= (3e-5).*rand(1,1) ;
V7Er= (2).*rand(1,1);
V8Er= (1).*rand(1,1) ;
V9Er= (2.2).*rand(1,1) ;

x = xControlledWindStep;
y = yControlledWindStep ;
z = zControlledWindStep ;
Velocity = VelocityControlledWindStep;
alpha = alphaControlledWindStep  ;
beta = betaControlledWindStep ;

z_obse = [Velocity(2,:) ; alpha(2,:) ; beta(2,:) ; x(2,:) ; y(2,:) ; z(2,:)];
%L = length(xo);
xL = length (xo);
yL = length (z_obse (:,1));
w_error = zeros(9,1) ;
v_error = zeros (6,1);
x_a = [ xo ; w_error ; v_error];

errors = [0.01369 ; 0.01369 ; 0.01369 ;2e-9 ; 1e-4 ;2.6e-1 ; 2 ; 1 ; 2] ;  %Remember p,q,r is in rad
P = diag(errors); 
yerrors = [ W1yEr ;W2yEr ;W3yEr ;W4yEr ;W5yEr ;W6yEr];
xerrors =[V1Er; V2Er; V3Er ;V4yEr ;V5Er ;V6Er ; V7Er ; V8Er ; V9Er];

xo_dot = [-1.1701e-8 ; 6.256e-20 ; 9.5279e-12 ; 5.0037e-19 ; -1.44995e-20 ; -1.2725e-19 ; 45.0741 ; 7.5205e-17 ; -3.3172e-9];

error_argumented = [errors; xerrors; yerrors];
Pa = diag(error_argumented);

R = diag (yerrors);
Q = diag (xerrors);
L = length (x_a) ; 

u_dot = zeros(1,L);
v_dot = zeros(1,L);
w_dot = zeros(1,L);                                                                                                                                                                                                                 
roll_dot = zeros(1,L);
pitch_dot = zeros(1,L);
yaw_dot = zeros(1,L);
x_dot = zeros(1,L);
y_dot = zeros(1,L);
h_dot = zeros(1,L);
V = zeros(1,L);
alpha = zeros(1,L);
beta =zeros(1,L);

u_doty = zeros(1,L);
v_doty = zeros(1,L);
w_doty = zeros(1,L);                                                                                                                                                                                                                 
roll_doty = zeros(1,L);
pitch_doty = zeros(1,L);
yaw_doty = zeros(1,L);
x_doty = zeros(1,L);
y_doty = zeros(1,L);
h_doty = zeros(1,L);

X_int =zeros(xL,2*L);
x_dot_kminus1 = zeros(xL,2*L);

for i = 1: 2*L
x_dot_kminus1 (:,i) = xo_dot;
X_int (:,i) = xo;
end

y_dot_kminus1= x_dot_kminus1;
y_int =X_int;

x = zeros (L,N) ; 
x(:,1) = x_a;
e = zeros (L,2*L);
D = ones(L,1);
e(:,1:L) = diag(D);
e(:, L+1:2*L) = - diag(D);
eta = sqrt(L)*e;
xhatkminus1_kminus1 = x(:,1);

for k = 2:N

if XXX(k) == 0  
    z_obse(:,k) = zhatk_kminus1;
end
    
P_kminus1_kminus_1 = Pa;
Skminus1_kminus1 = chol(P_kminus1_kminus_1,'lower');  
X_kminus1_kminus1 = zeros (L,2*L);
%Cubature sigma points
    for i = 1:2*L
        X_kminus1_kminus1(:,i) = Skminus1_kminus1*eta(:,i) + xhatkminus1_kminus1 ;
    end  
omega = 1/(2*L);

u_state =  X_kminus1_kminus1(1,:);
v_state = X_kminus1_kminus1 (2,:);
w_state = X_kminus1_kminus1 (3,:);
Roll_state = X_kminus1_kminus1 (4,:);
Pitch_state = X_kminus1_kminus1 (5,:);
Yaw_state = X_kminus1_kminus1 (6,:);
x_state = X_kminus1_kminus1 (7,:);
y_state = X_kminus1_kminus1 (8,:);
z_state = X_kminus1_kminus1 (9,:); 
w_noise = X_kminus1_kminus1  (xL+1:2*xL,:);
v_noise = X_kminus1_kminus1  (2*xL+1:L,:);

for i = 1: 2*L
u_dot (i) = (r_meas(k-1)*v_state (i)) - (q_meas(k-1)*w_state(i)) - (g*sin(Pitch_state(i))) + (ax_meas(k-1)); 
v_dot (i) = (p_meas(k-1)*w_state (i)) - (r_meas(k-1)*u_state(i)) + (g*cos(Pitch_state(i))*sin(Roll_state(i)))+ (ay_meas(k-1));
w_dot (i)  = (q_meas(k-1)*u_state (i)) - (p_meas(k-1)*v_state(i)) + (g*cos(Pitch_state(i))*cos(Roll_state(i)))+ (az_meas(k-1));

%rotational kinematics
roll_dot (i)  = p_meas(k-1) + ((q_meas(k-1)*sin(Roll_state(i))*tan(Pitch_state(i))) + (r_meas(k-1)*cos(Roll_state(i))*tan(Pitch_state(i))));
pitch_dot (i)  = ((q_meas(k-1)*cos(Roll_state(i))) - (r_meas(k-1)*sin(Roll_state(i))));
yaw_dot(i)  = ((q_meas(k-1)*sin(Roll_state(i))) + (r_meas(k-1)*cos(Roll_state(i))))/ cos(Pitch_state(i));

%Position kinematics
x_dot(i)  = u_state(i)*cos(Pitch_state(i))*cos(Yaw_state(i))  + v_state(i)*((sin(Roll_state(i))*sin(Pitch_state(i))*cos(Yaw_state(i)))- (cos(Roll_state(i))*sin(Yaw_state(i)))) + w_state(i)*((cos(Roll_state(i))*sin(Pitch_state(i))*cos(Yaw_state(i))) +  (sin(Roll_state(i))*sin(Yaw_state(i)))) ;

y_dot (i) = u_state(i)*cos(Pitch_state(i))*sin(Yaw_state(i))  + v_state(i)*((sin(Roll_state(i))*sin(Pitch_state(i))*sin(Yaw_state(i))) + (cos(Roll_state(i))*cos(Yaw_state(i)))) + w_state(i)*((cos(Roll_state(i))*sin(Pitch_state(i))*sin(Yaw_state(i))) -  (sin(Roll_state(i))*cos(Yaw_state(i)))) ;

h_dot (i) = u_state(i)*sin(Pitch_state(i)) - v_state(i)*sin(Roll_state(i))*cos(Pitch_state(i)) - w_state(i)*cos(Roll_state(i))*cos(Pitch_state(i));

end

Xstar_k_kminus1_dot = [u_dot ; v_dot; w_dot; roll_dot; pitch_dot; yaw_dot; x_dot; y_dot; h_dot];

for j = 1:2*L   %goes from 1 to 18 ie through all columns
    X_int(:,j) =  X_int(:,j) + 0.5*(time(k) - time(k-1))*(x_dot_kminus1(:,j) + Xstar_k_kminus1_dot(:,j)); 
end
Xstar_k_kminus1 = X_int + w_noise;
x_dot_kminus1 = Xstar_k_kminus1_dot;
%Propagating the sigma points

Xsum = 0;
XProduct_star = 0;

for i = 1:2*L
    XProduct_star = XProduct_star + omega*Xstar_k_kminus1(:,i)*Xstar_k_kminus1(:,i)';
    Xsum = Xsum + omega*Xstar_k_kminus1(:,i);
end

xhatk_kminus1 = Xsum;

% Estimate the predicted state error covariance
P_k_kminus_1 = XProduct_star - xhatk_kminus1*xhatk_kminus1' ;

%Observation Update
Sk_k_1 = chol(P_k_kminus_1,'lower');
X_k_kminus1 = zeros (xL,2*L);
for i = 1:2*L
        X_k_kminus1(:,i) = Sk_k_1*eta(1:xL,i) + xhatk_kminus1;
end  

u_statey =  X_k_kminus1(1,:);
v_statey = X_k_kminus1 (2,:);
w_statey = X_k_kminus1 (3,:);
Roll_statey = X_k_kminus1 (4,:);
Pitch_statey = X_k_kminus1 (5,:);
Yaw_statey = X_k_kminus1 (6,:);


for i = 1: 2*L
u_doty (i) = (r_meas(k)*v_statey (i)) - (q_meas(k)*w_statey(i)) - (g*sin(Pitch_statey(i)))+ (ax_meas(k));
v_doty (i) = (p_meas(k)*w_statey (i)) - (r_meas(k)*u_statey(i)) + (g*cos(Pitch_statey(i))*sin(Roll_statey(i)))+ (ay_meas(k));
w_doty (i)  = (q_meas(k)*u_statey (i)) - (p_meas(k)*v_statey(i)) + (g*cos(Pitch_statey(i))*cos(Roll_statey(i)))+ (az_meas(k));

%rotational kinematics
roll_doty (i)  = p_meas(k) + ((q_meas(k)*sin(Roll_statey(i))*tan(Pitch_statey(i))) + (r_meas(k)*cos(Roll_statey(i))*tan(Pitch_statey(i))));
pitch_doty (i)  = ((q_meas(k)*cos(Roll_statey(i))) - (r_meas(k)*sin(Roll_statey(i))));
yaw_doty(i)  = ((q_meas(k)*sin(Roll_statey(i))) + (r_meas(k)*cos(Roll_statey(i))))/ cos(Pitch_statey(i));

%Position kinematics
x_doty(i)  = u_statey(i)*cos(Pitch_statey(i))*cos(Yaw_statey(i))  + v_statey(i)*((sin(Roll_statey(i))*sin(Pitch_statey(i))*cos(Yaw_statey(i)))- (cos(Roll_statey(i))*sin(Yaw_statey(i)))) + w_statey(i)*((cos(Roll_statey(i))*sin(Pitch_statey(i))*cos(Yaw_statey(i))) +  (sin(Roll_statey(i))*sin(Yaw_statey(i)))) ;

y_doty (i) = u_statey(i)*cos(Pitch_statey(i))*sin(Yaw_statey(i))  + v_statey(i)*((sin(Roll_statey(i))*sin(Pitch_statey(i))*sin(Yaw_statey(i))) + (cos(Roll_statey(i))*cos(Yaw_statey(i)))) + w_statey(i)*((cos(Roll_statey(i))*sin(Pitch_statey(i))*sin(Yaw_statey(i))) -  (sin(Roll_statey(i))*cos(Yaw_statey(i)))) ;

h_doty (i) = u_statey(i)*sin(Pitch_statey(i)) - v_statey(i)*sin(Roll_statey(i))*cos(Pitch_statey(i)) - w_statey(i)*cos(Roll_statey(i))*cos(Pitch_statey(i));

end

Ystar_k_kminus1_dot = [u_doty ; v_doty; w_doty; roll_doty; pitch_doty; yaw_doty; x_doty; y_doty; h_doty];

for j = 1:2*L   %goes from 1 to 18 ie through all columns
    y_int(:,j) =  y_int(:,j) + 0.5*(time(k) - time(k-1))*(y_dot_kminus1(:,j) + Ystar_k_kminus1_dot(:,j)); 
end
Ystar_k_kminus1 = y_int ;
y_dot_kminus1 = Ystar_k_kminus1_dot ;

u_statey =  Ystar_k_kminus1(1,:);
v_statey = Ystar_k_kminus1 (2,:);
w_statey = Ystar_k_kminus1 (3,:);
Roll_statey = Ystar_k_kminus1 (4,:);
Pitch_statey = Ystar_k_kminus1 (5,:);
Yaw_statey = Ystar_k_kminus1 (6,:);
x_statey = Ystar_k_kminus1 (7,:);
y_statey = Ystar_k_kminus1 (8,:);
z_statey = Ystar_k_kminus1 (9,:);

for w = 1: 2*L
V (w)= sqrt(u_statey(w)^2 + v_statey(w)^2 + w_statey(w)^2) ; 
alpha (w) = atan (w_statey(w)/u_statey(w));
beta (w) = asin(v_statey(w)/(sqrt(u_statey(w)^2 + v_statey(w)^2 + w_statey(w)^2)));

end

Z_k_kminus_1 = [V; alpha; beta; x_statey; y_statey; z_statey];

Z_k_kminus_1 = Z_k_kminus_1 + v_noise;

%Calculating the covariences
Zsum = 0;
Pzz_sum = 0;
Pxz_sum = 0;

for i = 1:2*L
    Pzz_sum = Pzz_sum + Z_k_kminus_1(:,i)*Z_k_kminus_1(:,i)';
    Pxz_sum = Pxz_sum + Xstar_k_kminus1(:,i)*Z_k_kminus_1(:,i)';
    Zsum = Zsum + Z_k_kminus_1(:,i);
end
zhatk_kminus1 = omega*Zsum;
Pzz = omega*Pzz_sum - zhatk_kminus1*zhatk_kminus1'; % +R;
Pxz = omega*Pxz_sum  - xhatk_kminus1*zhatk_kminus1' ;

%Calculating the kalman gain
Kk = Pxz/Pzz; 

x(1:xL,k) = xhatk_kminus1 + Kk*(z_obse(:,k) - zhatk_kminus1);
P = P_k_kminus_1 - Kk*Pzz*Kk';
xhatkminus1_kminus1 = xhatk_kminus1;

%re_argumenting
x_a = [ x(1:xL,k) ; w_error ; v_error];
xhatkminus1_kminus1 = [ xhatk_kminus1 ; w_error ; v_error];
x(:,k) = x_a;
Pa = zeros(L,L);
Pa(1:xL,1:xL) = P;
Pa(xL+1:2*xL , xL+1 : 2*xL) = Q ; 
Pa(2*xL+1:L , 2*xL+1:L) = R ; 


end
x_true(7,:);
x_true(8,:);
x_true(9,:);


MSExCKF(n)=(sum( (1/N)*(xeBef - x(7,1:k)).^2))^0.5;
MSEyCKF(n)= (sum((1/N)*(yeBef - x(8,1:k)).^2))^0.5;
MSEzCKF(n)= (sum((1/N)*(zeBef - x(9,1:k)).^2))^0.5;

MSEuCKF(n)= (sum((1/N)*(uBef - x(1,1:k)).^2))^0.5;
MSEvCKF(n)= (sum((1/N)*(x_true(2,:) - x(2,1:k)).^2))^0.5;
MSEwCKF(n)= (sum((1/N)*(x_true(3,:) - x(3,1:k)).^2))^0.5;

MSErollCKF(n)= (sum((1/N)*(x_true(4,:) - x(4,1:k)).^2))^0.5;
MSEpitchCKF(n)= (sum((1/N)*(x_true(5,:) - x(5,1:k)).^2))^0.5;
MSEyawCKF(n)= (sum((1/N)*(x_true(6,:) - x(6,1:k)).^2))^0.5;

end

RUNSS= 1:1:Runs;
figure (23)
subplot(3,3,1)
 plot(RUNSS, MSEuUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEuCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE u (ms^{-1})'); grid off; legend('UKF','CKF');
 %legend boxoff


subplot(3,3,2)
 plot(RUNSS, MSEvUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEvCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE v (ms^{-1})'); grid off; %legend('UKF','CKF');

subplot(3,3,3)
 plot(RUNSS, MSEwUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEwCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE w (ms^{-1})'); grid off; %legend('UKF','CKF');

subplot(3,3,4)
 plot(RUNSS, MSErollUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSErollCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE \Phi (^{o})'); grid off; %legend('UKF','CKF');
 
subplot(3,3,5)
 plot(RUNSS, MSEpitchUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEpitchCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE \Theta (^{o})'); grid off; 

subplot(3,3,6)
 plot(RUNSS, MSEyawUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEyawCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE \Psi (^{o})'); grid off; %legend('UKF','CKF');
 
RUNSS= 1:1:Runs; 
subplot(3,3,7)
 plot(RUNSS, MSExUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSExCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE x (m)'); grid off; %legend('UKF','CKF');

subplot(3,3,8)
 plot(RUNSS, MSEyUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEyCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE y (m)'); grid off; %legend('UKF','CKF');
 
subplot(3,3,9)
 plot(RUNSS, MSEzUKF,'k','MarkerSize', 0.1)
 hold on
 plot(RUNSS, MSEzCKF,'b','MarkerSize', 0.1)
 hold off
 xlabel('Run'); ylabel('RMSE h (m)'); grid off;% legend('UKF','CKF');
%  
 
 x_CKF = x;
figure(11)
subplot(3,3,1)
 plot(time(1:k), x_CKF(1,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(1,1:k),'b-.')
 hold on
 plot(time(1:k), uBef,'k--')
 xlabel('Time (s)'); ylabel('u (ms^{-1})'); grid off; legend('CKF','UKF','Measured');
 hold off

subplot(3,3,2)
 plot(time(1:k), x_CKF(2,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(2,1:k),'b-.')
 hold on
 plot(time(1:k), x_true(2,:),'k--')
 xlabel('Time (s)'); ylabel('v (ms^{-1})'); grid off;% legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,3)
 plot(time(1:k), x_CKF(3,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(3,1:k),'b-.')
 hold on
 plot(time(1:k), x_true(3,:),'k--')
 xlabel('Time (s)'); ylabel('w (ms^{-1})'); grid off; %legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,4)
 plot(time(1:k), x_CKF(4,1:k)*180/pi,'g--')
 hold on
 plot(time(1:k), x_UKF(4,1:k)*180/pi,'b-.')
 hold on
 plot(time(1:k), x_true(4,:)*180/pi,'k--')
 xlabel('Time (s)'); ylabel('\Phi (^{o})'); grid off; %legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,5)
 plot(time(1:k), x_CKF(5,1:k)*180/pi,'g--')
 hold on
 plot(time(1:k), x_UKF(5,1:k)*180/pi,'b-.')
 hold on
 plot(time(1:k), x_true(5,:)*180/pi,'k--')
 xlabel('Time (s)'); ylabel('\Theta (^{o})'); grid off; %legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,6)
 plot(time(1:k), x_CKF(6,1:k)*180/pi,'g--')
 hold on
 plot(time(1:k), x_UKF(6,1:k)*180/pi,'b-.')
 hold on
 plot(time(1:k), x_true(6,:)*180/pi,'k--')
 xlabel('Time (s)'); ylabel('\Psi (^{o})'); grid off; %legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,7)
 plot(time(1:k), x_CKF(7,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(7,1:k),'b-.')
 hold on
 plot(time(1:k), x_true(7,:),'k--')
 xlabel('Time (s)'); ylabel('x (m)'); grid off; legend('CKF','UKF','Measured'); box off; legend boxoff
 hold off
 
subplot(3,3,8)
 plot(time(1:k), x_CKF(8,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(8,1:k),'b-.')
 hold on
 plot(time(1:k), x_true(8,:),'k--')
 xlabel('Time (s)'); ylabel('y (m)'); grid off; %legend('CKF','UKF','Measured');
 hold off
 
subplot(3,3,9)
 plot(time(1:k), x_CKF(9,1:k),'g--')
 hold on
 plot(time(1:k), x_UKF(9,1:k),'b-.')
 hold on
 plot(time(1:k), x_true(9,:),'k--')
 xlabel('Time (s)'); ylabel('h (m)'); grid off; %legend('CKF','UKF','Measured');
 hold off 

 figure(50)
plot3(x_CKF(8,1:k),x_CKF(7,1:k),x_CKF(9,1:k),'g--');
hold on
plot3(x_UKF(8,1:k),x_UKF(7,1:k),x_UKF(9,1:k),'b-.');
hold on
plot3(x_true(8,:),x_true(7,:),x_true(9,:),'k--');
xlabel('Y distance (m)'); ylabel('X distance (m)');zlabel('Altitude(m)'); grid off; legend('CKF','UKF','True'); legend boxoff;
hold off
