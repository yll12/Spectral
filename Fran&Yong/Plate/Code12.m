%CODE 12

%Anisotropic Plate 

%case:
%Orthorombic:
%Propagation along the fibres (alpha=0)
%SV&L modes



%Barium Sodium Niobate
%Half Thickness
h=2.5e-3;

%Physical Parameters
rho=5300;  %kg/m3;

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

Vt=sqrt(C(6,6)/rho); 

% for N=70:5:75; %Integration points %Integration points    
    
N=90;
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
xi=k.*2*h; %range of x to be plotted

% for m=0:1:(size(k,2)-1);
m=0;
%Differential Operator

L11=D2*C(2,2)-k(1,(m+1))^2*C(4,4)*eye(N);
L12=D1*1i*k(1,(m+1))*(C(2,3)+C(4,4));
%L12=L21;
L22=-k(1,(m+1))^2*C(3,3)*eye(N)+D2*C(4,4);

Lp=[L11, L12; L12, L22];

L=Lp;

%BC's

S11=C(2,2)*D1;
S12=C(2,3)*1i*k(1,(m+1))*eye(N);
S21=C(4,4)*1i*k(1,(m+1))*eye(N);
S22=C(4,4)*D1;


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

p=zeros((2*N), size(k,2));

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

%storing data
%     if N==70
%         H=W; %data storaging matrix N==60
%     elseif N==75
%         J=W;  %data storaging matrix  N==70  
%     
%     end
    
% end

% %  restoring for plotting
% for m=0:1:(size(k,2)-1);
% for j=1:1:2*N
%   
%         fd(j,(m+1))=(H(j,1)*2*h)/(2e3*pi); %in MHz-mm
%         vp(j,(m+1))=(H(j,1)*1e-3)/(k(1,(m+1))); %in mm/us
% end
% 
% end


% figure (2)
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

figure (2)
hold on
% for m=0:1:(size(k,2)-500);
    for j=1:1:30 %range of modes plotted 
plot( (2*h*k(1,(m+1))/pi) , ((2*h*bigw12cppn90(j,(m+1)))/(pi*Vt)), 'red *' );%plots real part of wavenumber, defined as in Graff (8.1.13)
hold on
    end
% end
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
% axis([0 20 0 20])
% grid on

% %DSIPERSE COMPARISON
% for m=0:1:(size(k,2)-500);
%     for j=1:1:40 %range of modes plotted 
% plot( (2*h*k(1,(m+1))) , ((H(j,(m+1)))/1e6), 'black' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
%     end
% end
% axis([0 10 0 7])
% grid on
