%CCODE 15

%SV & L Modes

%Indium (Tetragonal)
%radio interno (a) y externo (b) del cilindro
a=5; 
b=5.005; %m
%Physical Parameters
rho=5300;  %kg/m3;

%in N/m2
C(1,1) = 23.9e9; 
C(1,2) = 10.4e9; 
C(1,3) = 5e9; 
C(2,2) = 24.7e9;  
C(2,3) = 5.2e9; 
C(3,3) = 13.5e9; 
C(4,4) = 6.5e9;   
C(5,5) = 6.6e9;   
C(6,6) = 7.6e9;

N=90;


 Vt=sqrt(C(6,6)/rho);      
    
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

for m=0:1:(size(k,2)-1);

 %Differential Operator

L11= C(1,1)*(D2+diag(r.^-1)*D1)-k(1,(m+1))^2*C(5,5)*eye(N)-C(2,2)*diag(r.^-2);
L12=  1i*k(1,(m+1))*( (C(1,3)+C(5,5))*D1 + (C(1,3)-C(2,3))*diag(r.^-1) );

L21=  1i*k(1,(m+1))*( C(1,3)*D1+ C(2,3)*diag(r.^-1) +C(5,5)*(D1+diag(r.^-1)) );
L22=  C(5,5)*(D2+diag(r.^-1)*D1)-k(1,(m+1))^2*C(3,3)*eye(N);

Lp=[L11, L12; L21, L22];

L=Lp;

%BC's

S11= (C(1,1)*D1+C(1,2)*diag(r.^-1));
S12= (C(1,3)*1i*k(1,(m+1))*eye(N));
S21= (1i*k(1,(m+1))*C(5,5)*eye(N));
S22= (C(5,5)*D1);


S=[S11, S12; S21, S22];


%introducing the BC's in the problem
L(1,:)=S(1,:);
L(N,:)=S(N,:);
L((N+1),:)=S((N+1),:);
L((2*N),:)=S((2*N),:);

%Matrix M on the RHS
Mp=-rho*eye(N); Mp(1,1)=0; Mp(N,N)=0;

M=[Mp, O; O, Mp];


[P,E]=eig(L,M);
w=sort(real(sqrt(diag(E))));

for j=1:1:2*N
        W(j,(m+1))=w(j,1); 
        fd(j,(m+1))=(w(j,1)*h)/(2e3*pi); %in MHz-mm
        vp(j,(m+1))=(w(j,1)*1e-3)/(k(1,(m+1))); %in mm/us        
end

end


% figure (3)
% for m=0:1:(size(k,2)-500);
%     for j=1:1:30 %range of modes plotted 
% plot( ((h*k(1,(m+1))/pi)) , ((h*W(j,(m+1)))/(pi*Vt)), 'blue' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% 
% grid on
% 

figure (3)

for m=0:1:(size(k,2)-1);
    for j=1:1:30 %range of modes plotted 
plot( fd(j,(m+1)) , vp(j,(m+1)), 'black' );%plots real part of wavenumber, defined as in Graff (8.1.13)
hold on
    end
end
axis([0 10 0 6])
xlabel('Freq.-Thick. (MHz-mm)')
ylabel('Phase Vel. (m/ms)')
title('Orthorombic Cylinder')
% title('Orthorombic DISPERSE vs. Spectral Method')

grid on
