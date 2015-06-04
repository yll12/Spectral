%CODE 5
%Anisotropic Plate 

%case:

%TETRAGONAL7

%Propagation Parallel to the fibres (alpha=0)

%Half Thickness
h=50e-3;

rho=8938.4; %kg/m3;
% %in N/m2

c(1,1) =  235.54e9;
c(1,2) =  88.64e9;
c(1,3) =  88.78e9;
c(1,4) =  19.05e9;
c(1,5) =  -2.97e9;
c(1,6) =  -6.54e9;

c(2,2) =  215.85e9; 
c(2,3) =  108.43e9;
c(2,4) =  -14.14e9; 
c(2,5) =  2.35e9;
c(2,6) =  22.16e9; 

c(3,3) =  215.77e9;
c(3,4) =  -4.90e9;
c(3,5) =  0.61e9;
c(3,6) =  -15.57e9;

c(4,4) =  58.24e9;      
c(4,5) =  -17.68e9;  
c(4,6) =  -7.02e9; 

c(5,5) = 39.69e9;
c(5,6) = 11.66e9; 

c(6,6) = 44.02e9; 

alpha=0; %any angle
C=c;
%[C]=rotated_C_matrix(alpha, c);

ka= C(5,5)/C(6,6);
kb= C(5,6)/C(6,6);
kc= C(4,5)/C(6,6);
kd= C(4,6)/C(6,6);
ke= C(2,5)/C(6,6);
kf= C(2,6)/C(6,6);
kg= C(3,5)/C(6,6);
kh= C(3,6)/C(6,6);
ki= C(4,4)/C(6,6);
kj= C(2,4)/C(6,6);
kk= C(2,2)/C(6,6);
kl= C(3,4)/C(6,6);
km= C(2,3)/C(6,6);
kn= C(3,3)/C(6,6);

kmin = 0;
kmax = 18;
k = kmin:(kmax-kmin)/4.2e3:kmax;
k = k.^3;

%%G
N=50; %Integration points %Integration points    
    
% generate 
% generate Chebyshev differentiation matrices
[y,D]=chebdif(N,2);
y1=((y+1)/2)*h; %y1

D1=(1/2)^-1*D(:,:,1);
D2=(1/2)^-2*D(:,:,2);


O=zeros(N);


% for m=150:20:(size(k,2)-1500);
% for 
m=1000;

    
%Differential Operator

L11= -(h*k(1,(m+1)))^2*ka*eye(N) +2*kb*1i*(h*k(1,(m+1)))*D1      +D2; 
L12= -(h*k(1,(m+1)))^2*kc*eye(N) +(kd+ke)*1i*(h*k(1,(m+1)))*D1 +kf*D2; 
L13= -(h*k(1,(m+1)))^2*kg*eye(N) +(kh+kc)*1i*(h*k(1,(m+1)))*D1 +kd*D2; 

% L21= L12; 
L22= -(h*k(1,(m+1)))^2*ki*eye(N) +2*kj*1i*(h*k(1,(m+1)))*D1       +kk*D2; 
L23= -(h*k(1,(m+1)))^2*kl*eye(N) +(km+ki)*1i*(h*k(1,(m+1)))*D1 +kj*D2; 

% L31= L13; 
% L32= L23; 
L33= -(h*k(1,(m+1)))^2*kn*eye(N) +2*kl*1i*(h*k(1,(m+1)))*D1        +ki*D2; 

Lp=[L11, L12, L13; L12, L22, L23; L13, L23, L33];

L=Lp;

%BC's

%T2
S11= kf*h*D1  + 1i*((h^2)*k(1,(m+1)))*ke*eye(N);
S12= kk*h*D1 + 1i*((h^2)*k(1,(m+1)))*kj*eye(N);
S13= kj*h*D1  + 1i*((h^2)*k(1,(m+1)))*km*eye(N);

%T4
S21= kd*h*D1 + 1i*((h^2)*k(1,(m+1)))*kc*eye(N);
S22= kj*h*D1  + 1i*((h^2)*k(1,(m+1)))*ki*eye(N);
S23= ki*h*D1  + 1i*((h^2)*k(1,(m+1)))*kl*eye(N);

%T6
S31= h*D1       + 1i*((h^2)*k(1,(m+1)))*kb*eye(N);
S32= kf*h*D1  + 1i*((h^2)*k(1,(m+1)))*kd*eye(N);
S33= kd*h*D1 + 1i*((h^2)*k(1,(m+1)))*kh*eye(N);

S=[S31, S32, S33; S11, S12, S13; S21, S22, S23];

%introducing the BC's in the problem

L(1,:)=S(1,:);
L(N,:)=S(N,:);

L((N+1),:)=S((N+1),:);
L((2*N),:)=S((2*N),:);

L(((2*N)+1),:)=S(((2*N)+1),:);
L((3*N),:)=S((3*N),:);


%Matrix M on the RHS
Mp=-((h^2)*rho/C(6,6))*eye(N); Mp(1,1)=0; Mp(N,N)=0;

M=[Mp, O, O; O, Mp, O; O, O, Mp];


[P,E]=eig(L,M);

[w,ii]=sort(real(sqrt(diag(E))));
p1= P(:, ii);

p(:,(m+1))=w(:,1); 
q=1;
for n=1:3*N
  
    if p(n,(m+1))==0;
        q=q;
    else
        W(q,(m+1))=p(n,(m+1));
        P1(:,q)=p1(:,n);
        q=q+1;
    end
end    

for j=1:1:12;
       
       fd(j,(m+1))=(W(j,(m+1))*h)/(2e3*pi);            %in MHz-mm
       vp(j,(m+1))=(W(j,(m+1))*1e-3)/(k(1,(m+1)));     %in mm/us    
       
       figure (1)
       plot( fd(j,(m+1)) , vp(j,(m+1)), 'red *' );%plots real part of wavenumber, defined as in Graff (8.1.13)
       axis([0 5 0 10])
       hold on
       grid on
end    

% end



% figure (1)
% 
% % for m=150:5:(size(k,2)-1500);
%  for j=1:1:12 %range of modes plotted 
%  
% if j==12;        
% plot( fd(j,(150:size(k,2)-1500)) , vp(j,(150:size(k,2)-1500)), '-b','linewidth', 2);%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
% else
% plot( fd(j,(150:size(k,2)-1500)) , vp(j,(150:size(k,2)-1500)), '-b','linewidth', 2);%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
% end
% 
%  end
% % end
% axis([0 5 0 10])
% 
% grid on


m=j;

Ux=zeros(N,m);
Uy=zeros(N,m);
Uz=zeros(N,m);

nux=zeros(1,m);
nuy=zeros(1,m);
nuz=zeros(1,m);

NUx=zeros(N,m);
NUy=zeros(N,m);
NUz=zeros(N,m);



%The mode shape on the left of the plot corresponds to the highest mode for
%that frequency, as we move to the right of the plot we reach, on the right
%hand side the lowest mode for that frequency.
for i=1:1:m;

%Solid
Ux(1:N,i)= real(P1(1:1:N,i))-imag(P1(1:1:N,i));
Uy(1:N,i)= real(P1((N+1):1:2*N,i))-imag(P1((N+1):1:2*N,i));
Uz(1:N,i)= real(P1((2*N+1):1:3*N,i))-imag(P1((2*N+1):1:3*N,i));


nux(1,i)=max(abs(Ux(:,i)));
nuy(1,i)=max(abs(Uy(:,i)));
nuz(1,i)=max(abs(Uz(:,i)));


NUx(:,i)= Ux(:,i)/nux(1,i);
NUy(:,i)= Uy(:,i)/nuy(1,i);
NUz(:,i)= Uz(:,i)/nuz(1,i);


figure(3)
%Solid
subplot(3,m,i)
plot( NUx(1:N,i), y1(:,1), 'red')
hold on
axis([-1 1 0 0.05])
grid on

if i==m;
title('Ur displacement Profiles')
end

figure(3)
%steel
subplot(3,m,(i+m))
plot( NUy(1:N,i), y1(:,1), 'black')
axis([-1 1 0 0.05])
hold on
grid on

if i==m;
title('U0 displacement Profiles')
end

figure(3)
%steel
subplot(3,m,(i+(2*m)))
plot( NUz(1:N,i), y1(:,1), 'blue')
axis([-1 1 0 0.05])
hold on
grid on

if i==m;
title('Uz displacement Profiles')
end

end

