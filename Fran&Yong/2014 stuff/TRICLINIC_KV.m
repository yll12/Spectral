%parameters
clear all
clear workspace


h=1e-3; %plate thickness in m

rho=1500;  %kg/m3;


% %in N/m2
e(1,1) = 2.18e8;  
e(1,2) = 7.65e7; 
e(1,3) = 1.64e7; 
e(1,4) = -3.60e6;
e(1,5) = 6.88e5;
e(1,6) = 1.16e8; 

e(2,2) = 7.11e7;    
e(2,3) = 1.92e7; 
e(2,4) = -7.71e5;
e(2,5) = 2.15e6;
e(2,6) = 5.00e7;

e(3,3) = 4.22e7;      
e(3,4) = -9.644e5;
e(3,5) = 6.27e5;
e(3,6) = -3.07e6;

e(4,4) = 1.11e7;     
e(4,5) = 2.89e6;
e(4,6) = -1.15e6;

e(5,5) = 1.36e7;      
e(5,6) = 1.48e6;

e(6,6) = 9.35e7;


% % %in N/m2 %with 0.01 tolerance
% e(1,1) = 0;  
% e(1,2) = 0; 
% e(1,3) = 0; 
% e(1,4) = 0;
% e(1,5) = 0;
% e(1,6) = 0; 
% 
% e(2,2) = 0;  
% e(2,3) = 0; 
% e(2,4) = 0;
% e(2,5) = 0;
% e(2,6) = 0;
% 
% e(3,3) = 0;     
% e(3,4) = 0;
% e(3,5) = 0;
% e(3,6) = 0;
% 
% e(4,4) = 0;   
% e(4,5) = 0;
% e(4,6) = 0;
% 
% e(5,5) = 0;   
% e(5,6) = 0;
% 
% e(6,6) = 0;

%the matrices have been obtained from the orthogonal matrix: material 1, by
%means of a rotation with:beta=0.5764 and delta=0.3256

% %in N/m2
% c(1,1) = 74.29e9   -1i*(w(1,m)/(2*pi*2e6))*e(1,1);  
% c(1,2) = 28.94e9   -1i*(w(1,m)/(2*pi*2e6))*e(1,2);
% c(1,3) = 5.86e9     -1i*(w(1,m)/(2*pi*2e6))*e(1,3); 
% c(1,4) = 0.20e9     -1i*(w(1,m)/(2*pi*2e6))*e(1,4);
% c(1,5) = -0.11e9    -1i*(w(1,m)/(2*pi*2e6))*e(1,5);
% c(1,6) = 37.19e9   -1i*(w(1,m)/(2*pi*2e6))*e(1,6); 

 
%integration points
N=120;%checked for 100 and 120

[y,D]=chebdif(N,2);
y=((y+1)/2)*h; %y
y=y/h; %y hat
D1=((1/2)^-1)*D(:,:,1);
D2=((1/2)^-2)*D(:,:,2);

% withdamping
wmin = 0;
wmax = 112;
w = wmin:(wmax-wmin)/8.3e3:wmax;
w = w.^4;

c66=5.9472e10; 

for m=400:50:(size(w,2)-2000); %from 5000 onwards something is wrong

C(2,2) = 25.69e9   -1i*(w(1,m)/(2*pi*2e6))*e(2,2);
C(2,3) = 5.65e9    -1i*(w(1,m)/(2*pi*2e6))*e(2,3); 
C(2,4) = 9.28e7    -1i*(w(1,m)/(2*pi*2e6))*e(2,4);
C(2,5) = -8.01e7   -1i*(w(1,m)/(2*pi*2e6))*e(2,5);
C(2,6) = 17.52e9   -1i*(w(1,m)/(2*pi*2e6))*e(2,6);    

C(3,3) =   12.11e9  -1i*(w(1,m)/(2*pi*2e6))*e(3,3);     
C(3,4) =   1.33e7   -1i*(w(1,m)/(2*pi*2e6))*e(3,4);   
C(3,5) =   -0.86e7  -1i*(w(1,m)/(2*pi*2e6))*e(3,5);
C(3,6) =   0.22e9   -1i*(w(1,m)/(2*pi*2e6))*e(3,6);
  
C(4,4) =  4.18e9  -1i*(w(1,m)/(2*pi*2e6))*e(4,4);   
C(4,5) =  1.31e9  -1i*(w(1,m)/(2*pi*2e6))*e(4,5);
C(4,6) =  9.49e7  -1i*(w(1,m)/(2*pi*2e6))*e(4,6);   
    
C(5,5) =  5.35e9  -1i*(w(1,m)/(2*pi*2e6))*e(5,5);   
C(5,6) =  -7.05e7 -1i*(w(1,m)/(2*pi*2e6))*e(5,6);

C(6,6) =  28.29e9 -1i*(w(1,m)/(2*pi*2e6))*e(6,6);

   
k1=C(2,2)/C(6,6);%kappa1
k2=C(2,3)/C(6,6);%kappa2
k3=C(2,4)/C(6,6);%kappa3
k4=C(2,5)/C(6,6);%kappa4
k5=C(2,6)/C(6,6);%kappa5

k6=C(3,3)/C(6,6);%kappa6
k7=C(3,4)/C(6,6);%kappa7
k8=C(3,5)/C(6,6);%kappa8
k9=C(3,6)/C(6,6);%kappa9

k10=C(4,4)/C(6,6);%kappa10
k11=C(4,5)/C(6,6);%kappa11
k12=C(4,6)/C(6,6);%kappa12

k13=C(5,5)/C(6,6);%kappa13
k14=C(5,6)/C(6,6);%kappa14

what= w(1,m)*(h/(sqrt(C(6,6)/rho))); %w hat

%Q matrices

Z=zeros(N);

%2nd degree in K_hat
Q2=[-eye(N)*k13,      -eye(N)*k11,     -eye(N)*k8;
        -eye(N)*k11,      -eye(N)*k10,     -eye(N)*k7;
        -eye(N)*k8,         -eye(N)*k7,      -eye(N)*k6];

S2=zeros(3*N);
  
Q2(1,:)=S2(1,:);
Q2(N,:)=S2(N,:);
Q2((N+1),:)=S2((N+1),:);
Q2((2*N),:)=S2((2*N),:);
Q2((2*N+1),:)=S2((2*N+1),:);
Q2((3*N),:)=S2((3*N),:);

%1st degree in K_hat

Q1=[1i*(k14+k14)*D1,    1i*(k12+k4)*D1,      1i*(k9+k11)*D1;
        1i*(k12+k4)*D1,      1i*(k3+k3)*D1,        1i*(k2+k10)*D1;
        1i*(k9+k11)*D1,      1i*(k2+k10)*D1,      1i*(k7+k7)*D1];

    
S1=[1i*k4*eye(N),      1i*k3*eye(N),       1i*k2*eye(N);
        1i*k11*eye(N),    1i*k10*eye(N),     1i*k7*eye(N);
        1i*k14*eye(N),    1i*k12*eye(N),     1i*k9*eye(N)];
    

Q1(1,:)=S1(1,:);
Q1(N,:)=S1(N,:);
Q1((N+1),:)=S1((N+1),:);
Q1((2*N),:)=S1((2*N),:);
Q1((2*N+1),:)=S1((2*N+1),:);
Q1((3*N),:)=S1((3*N),:);


%0th degree in K_hat
L0u =     D2  +(what)^2*eye(N);
L0v =k1*D2  +(what)^2*eye(N);
L0w=k10*D2+(what)^2*eye(N);

Q0=[L0u,        k5*D2,    k12*D2;
        k5*D2,    L0v,        k3*D2;
        k12*D2,  k3*D2,    L0w];
    
S0=[ k5*D1,      k1*D1,   k3*D1;
         k12*D1,    k3*D1,   k10*D1;
         D1,           k5*D1,    k12*D1];

Q0(1,:)=S0(1,:);
Q0(N,:)=S0(N,:);
Q0((N+1),:)=S0((N+1),:);
Q0((2*N),:)=S0((2*N),:);
Q0((2*N+1),:)=S0((2*N+1),:);
Q0((3*N),:)=S0((3*N),:);


%companion matrices
M1 = [-Q1, -Q0; eye(3*N), zeros(3*N)];
M2 = [Q2, zeros(3*N); zeros(3*N), eye(3*N)];

[V, k] = eig(M1, M2); 

khat =diag(k);

kh=sort(khat);

kh=kh(find(isfinite(real(kh))==1));

kh=kh(find(real(kh)>=0));

kh=kh(find(imag(kh)>-1e-9));


a=real(kh)/h;
b=imag(kh)/h;

n=1;

% for j=1:1:3*N  
for j=1:1:size(a,1)       
    %FILTER
    if abs(b(j,1)/a(j,1))<=0.5;
        A(n,1)=abs(a(j,1));
        B(n,1)=abs(b(j,1));
        n=n+1;
    else
        n=n;
    end
    
end

% figure(1)
% hold on
% plot( A(1:size(A,1),1), B(1:size(A,1),1), 'red *')
% xlabel('Re(k)'), ylabel('Im(k)')
% grid on
% hold on

% figure(5)
% hold on
% plot((w(1,m)*(h/(pi*sqrt(c66/rho)))), ((1/pi)*mods(1:size(mods,1),1)),   'red o')
% axis([0 6 0 6])
% xlabel('w'), ylabel('Re(k)')
% grid on
% hold on

for j=1:1:size(A,1)
    
f(j,m)=(w(1,m))/(2e6*pi); %in MHz-mm
vp(j,m)=(w(1,m)*1e-3)/(A(j,1)); %in mm/us   

figure(1)
hold on
plot( f(j,m)*1e6 , vp(j,m)*1e3, 'red o'  );
axis([0 5*1e6 0 10*1e3])
axis([0 8*1e6 0 12*1e3])
xlabel('Freq.-Thick. (MHz-mm)')
ylabel('Phase Vel. (mm/us)')
grid on
hold on




% figure(2)
% hold on
% plot((w(1,m)*(h/(pi*sqrt(c66/rho)))), ((1/pi)*abs(real(kh(1:N,1)))), 'red o')
% axis([0 15 0 15])
% xlabel('w'), ylabel('Re(k)')
% grid on
% hold on

% figure(2)
% hold on
% plot((w(1,m))/(2e3*pi), abs(imag(kh(1:N,1))),  'blue o')
% axis([0 4000 0 100])
% xlabel('w'), ylabel('Im(k)')
% grid on
% hold on

% figure(3)
% hold on
% plot((w(1,m)*(h/(pi*sqrt(5.9472e10/rho)))), A((1:size(A,1)),1)*(h/(pi)),  'red')
% axis([0 8 0 8])
% xlabel('w')
% ylabel('Re(k)')
% grid on
% hold on

figure(2)
hold on
plot( (w(1,m)*1e6)/(2e6*pi) , abs(B(j,1)), 'red o'  );
axis([0 5*1e6 0 150])
axis([0 8*1e6 0 200])
xlabel('Freq. (MHz)')
ylabel('Attenuation. (Np/m)')
title('SAFE (red circles) vs. Spectral Kelvin-Voigt (blue lines): Lamb Modes Attenuation, Material 1')
grid on
hold on


end
end

