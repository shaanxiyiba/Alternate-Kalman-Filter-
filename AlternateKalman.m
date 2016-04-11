SIGNOISE=300.; 
X=0.;
Y=0.; 
XH=1000.; 
YH=2000.; 
XDH=0.; 
YDH=0.; 
XR1=1000000.; 
YR1=20000.*3280.; 
XR2=50000000.; 
YR2=20000.*3280.; 
ORDER=4; TS=1.; 
TF=20000.; 
PHIS=0.; 
T=0.; 
S=0.; 
H=.01; 
PHI=zeros(ORDER,ORDER); 
P=zeros(ORDER,ORDER); 
IDNP=eye(ORDER); 
Q=zeros(ORDER,ORDER); 
P(1,1)=1000.^2; 
P(2,2)=100.^2; 
P(3,3)=2000.^2; 
P(4,4)=100.^2; 
RMAT(1,1)=SIGNOISE^2; 
RMAT(1,2)=0.; 
RMAT(2,1)=0.; 
RMAT(2,2)=SIGNOISE^2; 
TS2=TS*TS; 
TS3=TS2*TS;
Q(1,1)=PHIS*TS3/3.; 
Q(1,2)=PHIS*TS2/2.; 
Q(2,1)=Q(1,2); 
Q(2,2)=PHIS*TS; 
Q(3,3)=PHIS*TS3/3.; 
Q(3,4)=PHIS*TS2/2.; 
Q(4,3)=Q(3,4); 
Q(4,4)=PHIS*TS; 
count=0; 
while T<=TF
    XR1OLD=XR1; 
    XR2OLD=XR2; 
    XOLD=X; 
    YOLD=Y; 
    XR1D=-14600.; 
    XR2D=-14600.; 
    XD=100.; 
    YD=0.; 
    XR1=XR1+H*XR1D; 
    XR2=XR2+H*XR2D; 
    X=X+H*XD; 
    Y=Y+H*YD; 
    T=T+H; 
    XR1D=-14600.; 
    XR2D=-14600.; 
    XD=100.; 
    YD=0.; 
    XR1=.5*(XR1OLD+XR1+H*XR1D); 
    XR2=.5*(XR2OLD+XR2+H*XR2D); 
    X=.5*(XOLD+X+H*XD); 
    Y=.5*(YOLD+Y+H*YD); 
    S=S+H;
    if S>=(TS-.00001) 
        S=0.; 
        XB=XH+XDH*TS; 
        YB=YH+YDH*TS; 
        R1B=sqrt((XR1-XB)^2+(YR1-YB)^2); 
        R2B=sqrt((XR2-XB)^2+(YR2-YB)^2); 
        HMAT(1,1)=-(XR1-XB)/R1B; 
        HMAT(1,2)=0.; 
        HMAT(1,3)=-(YR1-YB)/R1B; 
        HMAT(1,4)=0.; 
        HMAT(2,1)=-(XR2-XB)/R2B; 
        HMAT(2,2)=0.; 
        HMAT(2,3)=-(YR2-YB)/R2B; 
        HMAT(2,4)=0.; 
        HT=HMAT'; 
        PHI(1,1)=1.; 
        PHI(1,2)=TS; 
        PHI(2,2)=1.; 
        PHI(3,3)=1.; 
        PHI(3,4)=TS; 
        PHI(4,4)=1.; 
        PHIT=PHI'; 
        PHIP=PHI*P; 
        PHIPPHIT=PHIP*PHIT; 
        M=PHIPPHIT+Q; 
        HM=HMAT*M; 
        HMHT=HM*HT; 
        HMHTR=HMHT+RMAT; 
        HMHTRINV=inv(HMHTR); 
        MHT=M*HT; 
        GAIN=MHT*HMHTRINV; 
        KH=GAIN*HMAT; 
        IKH=IDNP-KH; 
        P=IKH*M; 
        R1NOISE=SIGNOISE*randn; 
        R2NOISE=SIGNOISE*randn; 
        R1=sqrt((XR1-X)^2+(YR1-Y)^2); 
        R2=sqrt((XR2-X)^2+(YR2-Y)^2);
        RES1=R1+R1NOISE-R1B; 
        RES2=R2+R2NOISE-R2B; 
        
        XH=XB+GAIN(1,1)*RES1+GAIN(1,2)*RES2; 
        XDH=XDH+GAIN(2,1)*RES1+GAIN(2,2)*RES2; 
        
        YH=YB+GAIN(3,1)*RES1+GAIN(3,2)*RES2; 
        YDH=YDH+GAIN(4,1)*RES1+GAIN(4,2)*RES2; 
        
        ERRX=X-XH; SP11=sqrt(P(1,1)); 
        ERRXD=XD-XDH; 
        
        SP22=sqrt(P(2,2)); 
        ERRY=Y-YH; 
        
        SP33=sqrt(P(3,3)); 
        ERRYD=YD-YDH; 
        
        SP44=sqrt(P(4,4)); 
        SP11P=-SP11; 
        SP22P=-SP22; 
        SP33P=-SP33; 
        SP44P=-SP44; 
        
        count=count+1; 
        
        %D = speed
        %H = location
        
        ArrayT(count)=T; 
        ArrayX(count)=X; %actual location
        ArrayXH(count)=XH; %estimated location
        ArrayXD(count)=XD; %constant speed
        ArrayXDH(count)=XDH; %estimated speed
        
        
        
        ArrayY(count)=Y; %reciever altitude
        ArrayYH(count)=YH; %reciever altitude estimate
        ArrayYD(count)=YD; %receiver altitude velocity
        ArrayYDH(count)=YDH; %reciever altitude velocity estimate
        
        
        ArrayERRX(count)=ERRX; %location error
        ArraySP11(count)=SP11; 
        ArraySP11P(count)=SP11P; 
        ArrayERRXD(count)=ERRXD; %speed error
        ArraySP22(count)=SP22; 
        ArraySP22P(count)=SP22P; 
        
        
        ArrayERRY(count)=ERRY; %reciever altitude error
        ArraySP33(count)=SP33; 
        ArraySP33P(count)=SP33P; 
        ArrayERRYD(count)=ERRYD;%receiver altitude velocity error 
        ArraySP44(count)=SP44; 
        ArraySP44P(count)=SP44P; 
    end
end

figure
subplot(2,1,1)
plot(ArrayX)
hold on
grid on
plot(ArrayXH,'red--');
xlabel('Time in seconds');
ylabel('Receiver Location');
title('Location Estimation using Alternative approach of Kalman filtering');
legend('Actual Location','Estimated Location');

subplot(2,1,2)
plot(ArrayXD)
hold on
grid on
plot(ArrayXDH,'red--');
xlabel('Time in seconds');
ylabel('Receiver Velocity');
title('Velocity Estimation using Alternative approach of Kalman filtering');
legend('Actual Velocity','Estimated Velocity');

figure
subplot(2,1,1)
plot(ArrayY)
hold on
grid on
plot(ArrayYH,'red--');
xlabel('Time in seconds');
ylabel('Receiver Altitude');
title('Altitude Estimation using Alternative approach of Kalman filtering');
legend('Actual Altitude','Estimated Altitude');

subplot(2,1,2)
plot(ArrayYD)
hold on
grid on
plot(ArrayYDH,'red--');
xlabel('Time in seconds');
ylabel('Receiver Altitude Velocity');
title('Altitude Velocity Estimation using Alternative approach of Kalman filtering');
legend('Actual Altitude Velocity','Estimated Altitude Velocity');

figure
subplot(2,1,1)
grid on
plot(ArrayERRX,'black--')
xlabel('Time in seconds');
ylabel('Receiver Location Error');
title('Location Error Estimation using Alternative approach of Kalman filtering');
legend('Estimated Error');

subplot(2,1,2)
grid on
plot(ArrayERRXD,'black--')
xlabel('Time in seconds');
ylabel('Receiver Velocity Error');
title('Velocity Error Estimation using Alternative approach of Kalman filtering');
legend('Estimated Error-Velocity');

figure
subplot(2,1,1)
grid on
plot(ArrayERRY,'black--')
xlabel('Time in seconds');
ylabel('Receiver Altitude Error');
title('Altitude Error Estimation using Alternative approach of Kalman filtering');
legend('Estimated Error-Altitude');

subplot(2,1,2)
grid on
plot(ArrayERRY,'black--')
xlabel('Time in seconds');
ylabel('Receiver Altitude Velocity Error');
title('Altitude Velocity Error Estimation using Alternative approach of Kalman filtering');
legend('Estimated Error-Altitude Velocity');
