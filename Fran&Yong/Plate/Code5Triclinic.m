%CODE 5
%Anisotropic Plate 

%case:

%TETRAGONAL7

%Propagation Parallel to the fibres (alpha=0)

%Half Thickness
h=2.5e-3;

% %TETRAGONAL6
% %Physical Parameters
% rho=7280;  %kg/m3;

% %in N/m2
% %in N/m2
% C(1,1) = 4.53e10;
% C(1,2) = 4.00e10; 
% C(1,3) = 4.15e10; 
% C(1,4) = 0;
% C(1,5) = 0;
% C(1,6) = 0; 
% C(2,2) = 4.53e10;   
% C(2,3) = 4.15e10; 
% C(2,5) = 0;
% C(2,6) = 0;
% C(3,3) = 4.51e10; 
% C(3,4) = 0;
% C(3,5) = 0;
% C(3,6) = 0;
% C(4,4) = 0.65e10;    
% C(4,5) = 0;
% C(4,6) = 0;
% C(5,5) = 0.65e10;   
% C(5,6) = 0;
% C(6,6) = 1.21e10;

%Physical Parameters
%TETRAGONAL7 
% rho=6950;  %kg/m3;
% 
% %in N/m2
% %in N/m2
% C(1,1) = 10.92e10; 
% C(1,2) = 6.83e10; 
% C(1,3) = 5.28e10; 
% C(1,4) = 0;
% C(1,5) = 0;
% C(1,6) = 1.36e10; 
% C(2,2) = 10.92e10;  
% C(2,3) = 5.28e10; 
% C(2,4) = 0;
% C(2,5) = 0;
% C(2,6) = -1.36e10;
% C(3,3) = 9.17e10; 
% C(3,4) = 0;
% C(3,5) = 0;
% C(3,6) = 0;
% C(4,4) = 2.67e10;   
% C(4,5) = 0;
% C(4,6) = 0;
% C(5,5) = 2.67e10;
% C(5,6) = 0;
% C(6,6) = 3.37e10;

%ORTHOROMBIC
% 
% rho=5300;  %kg/m3;
% 
% 
% %in N/m2
% %in N/m2
% C(1,1) = 23.9e10; 
% C(1,2) = 10.4e10; 
% C(1,3) = 5e10; 
% C(1,4) = 0;
% C(1,5) = 0;
% C(1,6) = 0; 
% C(2,2) = 24.7e10;   
% C(2,3) = 5.2e10;  
% C(2,5) = 0;
% C(2,6) = 0;
% C(3,3) = 13.5e10;  
% C(3,4) = 0;
% C(3,5) = 0;
% C(3,6) = 0;
% C(4,4) = 6.5e10;   
% C(4,5) = 0;
% C(4,6) = 0;
% C(5,5) = 6.6e10;   
% C(5,6) = 0;
% C(6,6) = 7.6e10;

%TRIGONAL

% % TELLURIUM
% 
% rho=6250;  %kg/m3;
% 
% %in N/m2

% C(1,1) = 3.27e10; 
% C(1,2) = 0.86e10; 
% C(1,3) = 2.49e10; 
% C(1,4) = 1.24e10;
% C(1,5) = 0;
% C(1,6) = 0; 
% C(2,2) = 3.27e10;   
% C(2,3) = 2.49e10;  
% C(2,4) = -1.24e10;
% C(2,5) = 0;
% C(2,6) = 0;
% C(3,3) = 7.22e10;   
% C(3,4) = 0;
% C(3,5) = 0;
% C(3,6) = 0;
% C(4,4) = 3.12e10;   
% C(4,5) = 0;
% C(4,6) = 0;
% C(5,5) = 3.12e10;   
% C(5,6) = 1.24e10;
% C(6,6) = (C(1,1)-C(1,2))/2;

% Matrix generated by rotating the orthorombic matrix of the SA case
% angles: beta=0.5764 and delta=0.3256

rho=8938.4; %kg/m3;
% 
% 
% %in N/m2

C(1,1) = 20.787e10;  
C(1,2) = 10.906e10; 
C(1,3) = 9.341e10; 
C(1,4) = 1.6574e10;
C(1,5) = -1.615e10;
C(1,6) = -2.3188e10; 
C(2,2) = 16.77e10;  
C(2,3) = 13.624e10; 
C(2,4) = -2.4719e10;
C(2,5) = 1.128e10;
C(2,6) = 0.86831e10;
C(3,3) = 18.591e10;     
C(3,4) = 0.81453e10;
C(3,5) = 0.82076e10;
C(3,6) = 1.4505e10;
C(4,4) = 10.023e10;   
C(4,5) = 1.4505e10;
C(4,6) = 0.58388e10;
C(5,5) = 3.5110e10;   
C(5,6) = 1.6574e10;
C(6,6) = 5.9472e10;

Vt=sqrt(5.9472e10/rho); 

N=90; %Integration points %Integration points    
    
% generate Chebyshev differentiation matrices
[x,D]=chebdif(N,2);
x=x*h;
D1=(h^-1)*D(:,:,1);
D2=(h^-2)*D(:,:,2);

O=zeros(N);

kmin = 0;
kmax = 35;
k = kmin:(kmax-kmin)/4.2e3:kmax;
k = k.^3;
xi=(kmax^3)*2*h; %range of x to be plotted

% for m=0:100:(size(k,2)-1);
% m=1300;
%Differential Operator
% for m=0:1:10;
m=0;
L11= -k(1,(m+1))^2*C(5,5)*eye(N) +2*C(5,6)*1i*k(1,(m+1))*D1            +C(6,6)*D2; 
L12= -k(1,(m+1))^2*C(4,5)*eye(N) +(C(4,6)+C(2,5))*1i*k(1,(m+1))*D1 +C(2,6)*D2; 
L13= -k(1,(m+1))^2*C(3,5)*eye(N) +(C(3,6)+C(4,5))*1i*k(1,(m+1))*D1 +C(4,6)*D2; 

% L21= L12; 
L22= -k(1,(m+1))^2*C(4,4)*eye(N) +2*C(2,4)*1i*k(1,(m+1))*D1        +C(2,2)*D2; 
L23= -k(1,(m+1))^2*C(3,4)*eye(N) +(C(2,3)+C(4,4))*1i*k(1,(m+1))*D1 +C(2,4)*D2; 

% L31= L13; 
% L32= L23; 
L33= -k(1,(m+1))^2*C(3,3)*eye(N) +2*C(3,4)*1i*k(1,(m+1))*D1        +C(4,4)*D2; 

Lp=[L11, L12, L13; L12, L22, L23; L13, L23, L33];

L=Lp;

%BC's

S11= C(2,6)*D1 + 1i*k(1,(m+1))*C(2,5)*eye(N);
S12= C(2,2)*D1 + 1i*k(1,(m+1))*C(2,4)*eye(N);
S13= C(2,4)*D1 + 1i*k(1,(m+1))*C(2,3)*eye(N);

S21= C(4,6)*D1 + 1i*k(1,(m+1))*C(4,5)*eye(N);
S22= C(2,4)*D1 + 1i*k(1,(m+1))*C(4,4)*eye(N);
S23= C(4,4)*D1 + 1i*k(1,(m+1))*C(3,4)*eye(N);

S31= C(6,6)*D1 + 1i*k(1,(m+1))*C(5,6)*eye(N);
S32= C(2,6)*D1 + 1i*k(1,(m+1))*C(4,6)*eye(N);
S33= C(4,6)*D1 + 1i*k(1,(m+1))*C(3,6)*eye(N);

S=[S31, S32, S33; S11, S12, S13; S21, S22, S23];

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

% w=w(find(imag(w)==0));

%for j=1:1:3*N;
%       W(j,(m+1))=w(j,1); 
%        fd(j,(m+1))=(w(j,1)*2*h)/(2e3*pi); %in MHz-mm
%        vp(j,(m+1))=(w(j,1)*1e-3)/(k(1,(m+1)));     %in mm/us       
%end

% %storing data
%     if N==30
%         F=W; %data storaging matrix N==30
%     elseif  N==40
%         G=W; %data storaging matrix N==40
%     elseif  N==50
%         H=W; %data storaging matrix N==50
%     elseif  N==60
%         I=W; %data storaging matrix  N==60
%     else
%         J=W;  %data storaging matrix  N==70
%            
%     end
% 
    p(:,(m+1))=w(:,1); 
  q=1;
for n=1:2*N
  
    if p(n,(m+1))==0;
        q=q;
    else
        W(q,(m+1))=p(n,(m+1));
        q=q+1;
    end
end

    
% end

% figure (2)
% 
% for m=0:1:(size(k,2)-500);
%     for j=2:1:25%range of modes plotted 
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*F(j,(m+1)))/(pi*Vt)), 'red' );        
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*G(j,(m+1)))/(pi*Vt)), 'yellow' );
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*H(j,(m+1)))/(pi*Vt)), 'green' );
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*I(j,(m+1)))/(pi*Vt)), 'blue' );
% %plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% % axis([0 12 0 8])
% grid on

figure (3)
% for m=0:100:(size(k,2)-500);
    for j=1:1:30 %range of modes plotted 
        plot( (2*h*k(1,(m+1))/pi) , ((2*h*bigw5cppn90(j,(m+1)))/(pi*Vt)), 'red *' );%plots real part of wavenumber, defined as in Graff (8.1.13)
    hold on
    end
% end
axis([0 10 0 10])
% axis([0 12 0 8])
grid on


% figure (3)
% 
% for m=1000:1:(size(k,2)-1);
%     for j=2:1:20 %range of modes plotted 
% plot( fd(j,(m+1)) , vp(j,(m+1)), 'black' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% axis([0 20 0 20])
% grid on

