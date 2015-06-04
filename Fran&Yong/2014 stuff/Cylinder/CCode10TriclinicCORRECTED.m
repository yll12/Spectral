%CCode10 Triclinic (Parallel)

%The stifness matrix has already been change accordingly to atch
%cylindrical coordinates when deriving the EOM in mapple. So we simply use
%the entries given in cartesian coordinates to define the stifness matrix.

%radio interno (a) y externo (b) del cilindro
a=5; 
b=5.005; %m


%Physical Parameters

rho=8938.4;  %kg/m3;

%in N/m2
c(1,1) = 182.1e9;
c(1,2) = 119.1e9;
c(1,3) = 111.8e9; 
c(1,4) = 26.2e9;
c(1,5) = -17.3e9;
c(1,6) = 0;
c(2,2) = 167.7e9;
c(2,3) = 126.2e9;
c(2,4) = -26.2e9;
c(2,5) = 17.3e9;
c(2,6) = 0;
c(3,3) = 174.9e9; 
c(3,4) = 0;
c(3,5) = 4.237e9;
c(3,6) = 0;
c(4,4) = 92.1e9;      
c(4,5) = 0;
c(4,6) = 17.3e9;
c(5,5) = 53.5e9;   
c(5,6) = 26.2e9;
c(6,6) = 67.6e9;

alpha=pi/9.3; %angle in radians

%[C]=rotated_C_matrix(alpha, c);

C = c;

%this generates a quasi triclinic matrix from a monoclinic one to compare with the plate case
%this matrix is in the cartesian system, but the change to cylindrical has already been taken 
%into acount when getting the expression for the equations from mappleN=70; 
    
Vt=sqrt(67.6e9/rho);
 
N=90;

% generate Chebyshev differentiation matrices
[x,D]=chebdif(N,2);
h=b-a;
r=(h*x+a+b)/2;
D1=(2/h)*D(:,:,1);
D2=(2/h)^2*D(:,:,2);

O=zeros(N);

kmin = 0;
kmax = 35;
k = kmin:(kmax-kmin)/4.2e3:kmax;
k = k.^3;
xi=k.*2*h; %range of x to be plotted

% for m=0:1:(size(k,2)-1);

m=0;
    
%Differential Operator

L11= C(2,2)*(D2+diag(r.^-1)*D1)  -C(1,1)*diag(r.^-2)*eye(N) +1i*k(1,(m+1))*C(2,4)*(2*D1+diag(r.^-1)*eye(N)) -k(1,(m+1))^2*C(4,4)*eye(N); 
L12= 1i*k(1,(m+1))*C(2,5)*(D1+diag(r.^-1)*eye(N))-1i*k(1,(m+1))*C(1,5)*diag(r.^-1)*eye(N) -k(1,(m+1))^2*C(4,5)*eye(N) +C(2,6)*D2 -C(1,6)*diag(r.^-1)*(D1-diag(r.^-1)*eye(N)) +1i*k(1,(m+1))*C(4,6)*(D1-diag(r.^-1)*eye(N));
L13= 1i*k(1,(m+1))*C(2,3)*(D1+diag(r.^-1)*eye(N)) -1i*k(1,(m+1))*C(1,3)*diag(r.^-1)*eye(N) -k(1,(m+1))^2*C(3,4)*eye(N) +C(2,4)*(D2+diag(r.^-1)*D1) -C(1,4)*diag(r.^-1)*D1 +1i*k(1,(m+1))*C(4,4)*D1;

L21= 1i*k(1,(m+1))*C(2,5)*D1 +C(2,6)*(D2+2*diag(r.^-1)*D1) +1i*k(1,(m+1))*C(1,5)*diag(r.^-1)*eye(N) +C(1,6)*(diag(r.^-1)*D1+diag(r.^-2)*eye(N)) -k(1,(m+1))^2*C(4,5)*eye(N) +1i*k(1,(m+1))*C(4,6)*(D1+2*diag(r.^-1)*eye(N));
L22= -k(1,(m+1))^2*C(5,5)*eye(N) +1i*k(1,(m+1))*C(5,6)*(2*D1+diag(r.^-1)*eye(N)) +C(6,6)*(D2+diag(r.^-1)*D1-diag(r.^-2)*eye(N));
L23= -k(1,(m+1))^2*C(3,5)*eye(N) +1i*k(1,(m+1))*C(3,6)*(D1+2*diag(r.^-1)*eye(N)) +1i*k(1,(m+1))*C(4,5)*D1 +C(4,6)*(D2+2*diag(r.^-1)*D1);

L31= 1i*k(1,(m+1))*C(2,3)*D1 +C(2,4)*(D2+diag(r.^-1)*D1) +1i*k(1,(m+1))*C(1,3)*diag(r.^-1)*eye(N) +C(1,4)*(diag(r.^-1)*D1) -k(1,(m+1))^2*C(3,4)*eye(N) +1i*k(1,(m+1))*C(4,4)*(D1+diag(r.^-1)*eye(N));
L32= -k(1,(m+1))^2*C(3,5)*eye(N) +1i*k(1,(m+1))*C(4,5)*(D1+diag(r.^-1)*eye(N)) +1i*k(1,(m+1))*C(3,6)*(D1-diag(r.^-1)*eye(N)) +C(4,6)*D2;
L33= -k(1,(m+1))^2*C(3,3)*eye(N) +1i*k(1,(m+1))*C(3,4)*(D1+diag(r.^-1)*eye(N)) +1i*k(1,(m+1))*C(3,4)*D1 +C(4,4)*(D2+diag(r.^-1)*D1);

Lp=[L11, L12, L13; L21, L22, L23; L31, L32, L33];

L=Lp;

%BC's

S11= C(2,2)*D1 +C(1,2)*diag(r.^-1)*eye(N)+1i*k(1,(m+1))*C(2,4)*eye(N); 
S12= 1i*k(1,(m+1))*C(2,5)*eye(N)+C(2,6)*(D1-diag(r.^-1)*eye(N));
S13= 1i*k(1,(m+1))*C(2,3)*eye(N)+C(2,4)*D1;

S21= C(2,4)*D1 +C(1,4)*diag(r.^-1)*eye(N)+1i*k(1,(m+1))*C(4,4)*eye(N);
S22= 1i*k(1,(m+1))*C(4,5)*eye(N)+C(4,6)*(D1-diag(r.^-1)*eye(N));
S23= 1i*k(1,(m+1))*C(3,4)*eye(N)+C(4,4)*D1;

S31= C(2,6)*D1 +C(1,6)*diag(r.^-1)*eye(N)+1i*k(1,(m+1))*C(4,6)*eye(N);
S32= 1i*k(1,(m+1))*C(5,6)*eye(N)+C(6,6)*(D1-diag(r.^-1)*eye(N));
S33= 1i*k(1,(m+1))*C(3,6)*eye(N)+C(4,6)*D1;

S=[S11, S12, S13; S21, S22, S23; S31, S32, S33];



%introducing the BC's in the problem
L(1,:)=S(1,:);
L(N,:)=S(N,:);
L((N+1),:)=S((N+1),:);
L((2*N),:)=S((2*N),:);
L(((2*N)+1),:)=S(((2*N)+1),:);
L((3*N),:)=S((3*N),:);


%Matrix M on the RHS
Mp=-rho*eye(N); Mp(1,1)=0; Mp(N,N)=0;

M=[Mp, O, O; O, Mp, O; O, O, Mp];


[P,E]=eig(L,M);
w=sort(real(sqrt(diag(E))));

for j=1:1:3*N
        W(j,(m+1))=w(j,1);        
end

% end

% hold on
% figure (1)
% for m=0:1:(size(k,2)-500);
%     for j=1:1:30 %range of modes plotted 
% plot( ((h*k(1,(m+1))/pi)) , ((h*W(j,(m+1)))/(pi*Vt)), 'blue' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% xlabel('2h/kpi')
% ylabel('2hw/Vpi')
% title('Triclinic Cylinder vs. Plate')
% 
% grid on