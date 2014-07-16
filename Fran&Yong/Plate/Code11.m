%CODE 11

%Anisotropic Plate 

%case:
%Orthorombic:
%Propagation along the fibres (alpha=0)
%SH modes

%Half Thickness
h=5e-3;
%Barium Sodium Niobate

%Physical Parameters


rho=5300;  %kg/m3;

% CARTESIAN
%in N/m2
C(1,1) = 23.9e10; 
C(1,2) = 10.4e10; 
C(1,3) = 5e10; 
C(2,2) = 24.7e10; 
C(2,3) = 5.2e10; 
C(3,3) = 13.5e10; 
C(4,4) = 6.5e10;    
C(5,5) = 6.6e10;   
C(6,6) = 7.6e10;

% rho=8938.4;    %kg/m3;
% 
% %in N/m2
% C(1,1) = 182.1e9; 
% C(1,2) = 119.1e9; 
% C(1,3) = 111.8e9; 
% C(2,2) = 167.7e9;  
% C(2,3) = 126.2e9; 
% C(3,3) = 174.9e9; 
% C(4,4) = 92.1e9;   
% C(5,5) = 53.5e9;   
% C(6,6) = 67.7e9;

Vt=sqrt(C(6,6)/rho);

kmin = 0;
kmax = 18;
k = kmin:(kmax-kmin)/4.2e3:kmax;
k = k.^3;
xi=k.*2*h; %range of x to be plotted


for N=70:5:70; %Integration points %Integration points    
    
% generate Chebyshev differentiation matrices
[x,D]=chebdif(N,2);
x=(x+1)*(h/2);
D1=((h/2)^-1)*D(:,:,1);
D2=((h/2)^-2)*D(:,:,2);

O=zeros(N);

for m=0:50:(size(k,2)-1);

    
%Differential Operator
Lp=-k(1,(m+1))^2*C(5,5)*eye(N)+C(6,6)*D2;

L=Lp;

%BC's
S=C(6,6)*D1;

%introducing the BC's in the problem
L(1,:)=S(1,:);
L(N,:)=S(N,:);

%Matrix M on the RHS
M=-rho*eye(N); M(1,1)=0; M(N,N)=0;

[P,E]=eig(L,M);
w=sort(real(sqrt(diag(E))));

p=zeros(N, size(k,2));

p(:,(m+1))=w(:,1); 
  q=1;
for n=1:N
  
    if p(n,(m+1))==0;
        q=q;
    else
        W(q,(m+1))=p(n,(m+1));
        q=q+1;
    end
end    
end

%storing data
    if N==70
        H=W; %data storaging matrix N==60
    elseif N==75
        J=W;  %data storaging matrix  N==70  
    
    end
    
end

% %  restoring for plotting
% for m=0:1:(size(k,2)-1);
% for j=1:1:N
%   
%         fd(j,(m+1))=(H(j,1)*2*h)/(2e3*pi); %in MHz-mm
%         vp(j,(m+1))=(H(j,1)*1e-3)/(k(1,(m+1))); %in mm/us
% end
% 
% end


% %check with different N
% figure (1)
% 
% for m=0:50:(size(k,2)-500);
%     for j=1:1:25%range of modes plotted 
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*H(j,(m+1)))/(pi*Vt)), 'red o' );        
% plot( (2*h*k(1,(m+1))/pi) , ((2*h*J(j,(m+1)))/(pi*Vt)), 'blue *' );
% %plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% % axis([0 12 0 8])
% grid on
% 
% hold on
figure(2)
for m=0:50:(size(k,2)-500);
    for j=1:1:25 %range of modes plotted 
plot( (h*k(1,(m+1))/pi) , ((h*H(j,(m+1)))/(pi*Vt)), 'black o' );%plots real part of wavenumber, defined as in Graff (8.1.13)
hold on
    end
end
axis([0 8 0 8])
grid on


% figure (3)
% 
% for m=1000:1:(size(k,2)-1);
%     for j=2:1:20 %range of modes plotted 
% plot( fd(j,(m+1)) , vp(j,(m+1)), 'black' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% xlabel('Freq.-Thick. (MHz-mm)')
% ylabel('Phase Vel. (m/ms)')
% title('Orthorombic DISPERSE vs. Spectral Method')
% axis([0 20 0 20])
% grid on

