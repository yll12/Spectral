clear all
clear workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%AXES CONFIGURATION%%%%%%%%%%%%%%%%%

 %Plane {x,z} === Plane of the plate
 %{y} axis ORTHOGONAL to the plane of the plate 
 %PROPAGATION ALONG {z} axis
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLATE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%% 

h=50e-3;   %plate thickness in m

rho=8938.4;  %kg/m3;

% %Stiffness Cartesian Matrixin in N/m2
% 
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


%%%%%%%%%%%%%%DAMPING CONSTANTS MATRIX%%%%%%%%%%%%%%%%%%%%%%
%in N/m2
% e(1,1) = 2.18e8;  
% e(1,2) = 7.65e7; 
% e(1,3) = 1.64e7; 
% e(1,4) = -3.60e6;
% e(1,5) = 6.88e5;
% e(1,6) = 1.16e8; 
% 
% e(2,2) = 7.11e7;    
% e(2,3) = 1.92e7; 
% e(2,4) = -7.71e5;
% e(2,5) = 2.15e6;
% e(2,6) = 5.00e7;
% 
% e(3,3) = 4.22e7;      
% e(3,4) = -9.644e5;
% e(3,5) = 6.27e5;
% e(3,6) = -3.07e6;
% 
% e(4,4) = 1.11e7;     
% e(4,5) = 2.89e6;
% e(4,6) = -1.15e6;
% 
% e(5,5) = 1.36e7;      
% e(5,6) = 1.48e6;
% 
% e(6,6) = 9.35e7;


% % % %in N/m2 %with 0.01 tolerance
e(1,1) = 0;  
e(1,2) = 0; 
e(1,3) = 0; 
e(1,4) = 0;
e(1,5) = 0;
e(1,6) = 0; 

e(2,2) = 0;  
e(2,3) = 0; 
e(2,4) = 0;
e(2,5) = 0;
e(2,6) = 0;

e(3,3) = 0;     
e(3,4) = 0;
e(3,5) = 0;
e(3,6) = 0;

e(4,4) = 0;   
e(4,5) = 0;
e(4,6) = 0;

e(5,5) = 0;   
e(5,6) = 0;

e(6,6) = 0;

%%%%%%%%%%%%%VELOCITY TO RESCALE THE PLOT%%%%%%%%%%%%%%%%%
% Vt=sqrt(/rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%HYSTERIC MODEL%%%%%%%%%%%%%%%%%%%%%
%Expression for the C omplex Stiffness Matrix which will enter the
%equations

%C=c-1i*e;
%Where 
%c is the cartesian stiffness matrix

%e in the damping constants matrix

%%%%%%%%%%%%%%%%%%CODE PARAMETERS%%%%%%%%%%%%%%%%%%%%%
%the matrices have been obtained from the orthogonal matrix: material 1, by
%means of a rotation with:beta=0.5764 and delta=0.3256 

% %in N/m2
% c(1,1) = 74.29e9   -1i*e(1,1);  
% c(1,2) = 28.94e9   -1i*e(1,2);
% c(1,3) = 5.86e9    -1i*e(1,3); 
% c(1,4) = 0.20e9    -1i*e(1,4);
% c(1,5) = -0.11e9   -1i*e(1,5);
% c(1,6) = 37.19e9   -1i*e(1,6); 

 
%integration points
N=83;%cheked for 80 and 120 MUST BE ODD FOR THE PLOT AND NORMSLIZATION OF PHASE

[y,D]=chebdif(N,2);
y=((y+1)/2)*h; %y
y=y/h; %y hat
D1=((1/2)^-1)*D(:,:,1);
D2=((1/2)^-2)*D(:,:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%ANGULAR FREQUENCY RANGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wmin = 0;
% wmax = 100;
% w = wmin:(wmax-wmin)/8.3e3:wmax;
% w = w.^4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%ANGULAR FREQUENCY RANGE SAFE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wmin = 1e3*2*pi; %fmin=1kHz
wmax = 0.1e6 *2*pi; %fmax=0.1MHz
step = 1e3*2*pi; %step=1kHz
w = wmin:step:wmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for m=400:25:(size(w,2)-5900); %from 5000 onwards something is wrong
% for m=1:1:size(w,2); %equally spaced w SAFE
for m=40;
C(:,:) = c(:,:)  -1i*e(:,:);

   
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

k15=C(6,6)/C(6,6);%kappa15

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
L0u = k15*D2  +(what)^2*eye(N);
L0v =k1*D2  +(what)^2*eye(N);
L0w=k10*D2+(what)^2*eye(N);

Q0=[L0u,        k5*D2,    k12*D2;
        k5*D2,    L0v,        k3*D2;
        k12*D2,  k3*D2,    L0w];
    
S0=[ k5*D1,      k1*D1,   k3*D1;
         k12*D1,    k3*D1,   k10*D1;
         k15*D1,           k5*D1,    k12*D1];

Q0(1,:)=S0(1,:);
Q0(N,:)=S0(N,:);
Q0((N+1),:)=S0((N+1),:);
Q0((2*N),:)=S0((2*N),:);
Q0((2*N+1),:)=S0((2*N+1),:);
Q0((3*N),:)=S0((3*N),:);


%companion matrices
M1 = [-Q1, -Q0; eye(3*N), zeros(3*N)];
M2 = [Q2, zeros(3*N); zeros(3*N), eye(3*N)];

[p, W] = eig(M1, M2); 

% khat =diag(k);
% kh=sort(khat);
% kh=kh(find(isfinite(real(kh))==1));
% kh=kh(find(real(kh)>=0));
% kh=kh(find(imag(kh)>-1e-9));
% 

[kh, ii]=sort(diag(W));

p = p(:, ii);

% khref=kh(find(isfinite(real(kh))==1));
% khref=khref(find(real(khref)>=0));
% khref=khref(find(imag(khref)>-1e-12));
clear KH;
clear P;

t=1;
for s=1:1:size(kh,1)
if isfinite(real(kh(s,1)))==1 && real(kh(s,1))>=0 && imag(kh(s,1))>-1e-6;
    KH(t,1)=kh(s,1);
    P(:,t) =p(:,s);
    t=t+1;
else
    t=t;
end
end    

%The filtring algorithm works khref and KH agree. Now we need to clasify
%the eigenvalues in propagating and non-propagating branches taking care of
%keeping all of them for further study. So we will not discard any now.

a=real(KH)/h;
b=imag(KH)/h;

q=1;
i=1;

clear A;
clear B;
clear NA;
clear NB;
clear PropMod;
clear NoPropMod;

for s=1:1:size(a,1)
    %FILTER FOR PROPAGATING MODES AND NON-PROPAGATING

    if abs(b(s,1)/a(s,1))<=0.01 
%         && a(s,1)<=2e3    
        A(q,1)=a(s,1);
        B(q,1)=b(s,1);
        PropMod(:,q)=P(:,s);
        q=q+1;
        i=i;
    else
        NA(i,1)=a(s,1);
        NB(i,1)=b(s,1);
        NoPropMod(:,i)=P(:,s);
        i=i+1;
        q=q;
    end
    
end


for j=1:1:size(A,1)
    
%Frequency 
fd(j,m)=(w(1,m)*h)/(2e3*pi); %in MHz
%Phase Velocity
vp(j,m)=(w(1,m)*1e-3)/(A(j,1)); %in mm/us   


%%%%%%%%%%%%%%%%%%%%%%%%%ADIMESIONAL PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% plot( ((h*A(s,1)/pi)) , ((h*w(1,m))/(pi*Vt)), 'black o' );%plots real part of wavenumber, defined as in Graff (8.1.13)
% hold on
% axis([0 8 0 4])
% xlabel('2·h·k / \pi')
% ylabel(' 2·h·\omega / \pi·V_{66} ' )
% title('Flexural Modes n=0 (lines) emph{vs}. Torsional and Longitudinal (circles) ')
% grid on   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%MODES PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%%SET UP WHITE BACKGROUND FOR FIGURE%%%%%%%%%%%%%%%%%%%%%%%%%
Modes=figure(1);
set(Modes    ,...  
'Color'    , 'white');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%PLOT SCM RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCM = plot( fd(j,m) , vp(j,m));
set( SCM                 ,... 
'LineWidth'    , 1       ,...
'LineStyle'    , 'none'  ,... 
'Marker'       , 'o'     ,...
'MarkerSize'   , 6       ,... 
'Color'        , 'b'     );
hold on

%%%%%%%%%%%%%%%%%%%%%%%%PLOT SAFE RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Here you must modify the SAFE plot line below so that you use your
%%%%%results. The first coordinate corresponds to the Freq-Thick product
%%%%%and the second coordinate is the value of the Imaginary Part of the
%%%%%Wavenumber in (Neper/m), that is, just plot the Imaginary Part of the
%%%%%Wavenumber in (1/m). 
% SAFE = plot( fd(j,m) , vp(j,m));
% set( SAFE                 ,... 
% 'LineWidth'    , 1       ,...
% 'LineStyle'    , 'none'  ,... 
% 'Marker'       , '*'     ,...
% 'MarkerSize'   , 5       ,... 
% 'Color'        , 'r'     );
% hold on

%%%%%%%%%%%%%%%Make an "I feel pretty" Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%      X AXIS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1label=xlabel('Freq.-Thick. (MHz-mm)');
set( X1label            ,...
'FontWeight' , 'bold'   ,... 
'FontAngle'  , 'default',...
'FontSize'   , 12       ,... 
'Color'      , 'k'      );

%%%%%%%%%%%%%%%%%%%%%%%%%      Y AXIS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1label=ylabel('Phase Velocity (m/ms)');
set( Y1label            ,...
'FontWeight' , 'bold'   ,... 
'FontAngle'  , 'default',...
'FontSize'   , 12       ,... 
'Color'      , 'k'      );

%%%%%%%%%%%%%%%%%%%%%%%%%%  AXES PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca                     , ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .00] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'XTick'       , 0:1:5     , ...
  'YTick'       , 0:1:10    , ...
  'LineWidth'   , 1         );

axis([0 5 0 10])
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%    FIGURE  TITLE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1Title=title({'Hysteric Model', '50mm Triclinic Plate'});
set( G1Title            , ...
'FontWeight' , 'bold'   ,... 
'FontAngle'  , 'default',...
'FontSize'   , 16       ,... 
'Color'      , 'k'      );



%%%%%%%%%%%%%%%INCLUDE A LITTLE LEGEND TO EXPLAIN THE GRAPH%%%%%%%%%%%%%%%%
% G1Legend = legend( ...
%   [SCM, SAFE]   , ...
%   'SCM'         , ...
%   'SAFE'        , ...
%   'location', 'NorthEast');


%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

%%%%%%%%%%%%%%%%%%%%%ATTENUATION PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%SET UP WHITE BACKGROUND FOR FIGURE%%%%%%%%%%%%%%%%%%%%%%%%%
Att=figure(2);
set(Att    ,...  
'Color'      , 'white');
hold on

%%%%%%%%%%%%%%%%%%%%%%%%PLOT SCM RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AttSCM = plot( (w(1,m)*h)/(2e3*pi) , B(j,1));
set( AttSCM                 ,... 
'LineWidth'    , 1       ,...
'LineStyle'    , 'none'  ,... 
'Marker'       , 'o'     ,...
'MarkerSize'   , 5       ,... 
'Color'        , 'b'     );
hold on

%%%%%%%%%%%%%%%%%%%%%%%%PLOT SAFE RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Here you must modify the SAFE plot line below so that you use your
%%%%%results. The first coordinate corresponds to the Freq-Thick product
%%%%%and the second coordinate is the value of the Imaginary Part of the
%%%%%Wavenumber in (Neper/m), that is, just plot the Imaginary Part of the
%%%%%Wavenumber in (1/m). 
AttSAFE = plot( (w(1,m)*h)/(2e3*pi) , B(j,1));
set( AttSAFE             ,... 
'LineWidth'    , 1       ,...
'LineStyle'    , 'none'  ,... 
'Marker'       , '*'     ,...
'MarkerSize'   , 5       ,... 
'Color'        , 'r'     );
hold on

%%%%%%%%%%%%%%%Make an "I feel pretty" Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%      X AXIS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X2label=xlabel('Freq.-Thick. (MHz-mm)');
set( X2label            , ...
'FontWeight' , 'bold'   ,... 
'FontAngle'  , 'default',...
'FontSize'   , 12       ,... 
'Color'      , 'k'      );

%%%%%%%%%%%%%%%%%%%%%%%%%      Y AXIS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y2label=ylabel('Attenuation (Np/m)');
set( Y2label            , ...
'FontWeight' , 'bold'   ,... 
'FontAngle'  , 'default',...
'FontSize'   , 12       ,... 
'Color'      , 'k'      );

%%%%%%%%%%%%%%%%%%%%%%%%%%  AXES PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca                     ,...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .00] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'XTick'       , 0:1:5     , ...
  'YTick'       , 0:0.2:1  , ...
  'LineWidth'   , 1         );

axis([0 5 0 1])
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%    FIGURE  TITLE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G2Title=title({'Hysteric Model', '50mm Triclinic Plate: Attenuation'});
set( G2Title             , ...
'FontWeight' , 'bold'    ,... 
'FontAngle'  , 'default' ,...
'FontSize'   , 16        ,... 
'Color'      , 'k'      );

%%%%%%%%%%%%%%%INCLUDE A LITTLE LEGEND TO EXPLAIN THE GRAPH%%%%%%%%%%%%%%%%
G2Legend = legend( ...
  [AttSCM, AttSAFE]   , ...
  'SCM'               , ...
  'SAFE'              , ...
  'location', 'NorthWest');
 
end
end


%Number of propagating modes
m=size(PropMod,2);

Ur=zeros(N,m);
U0=zeros(N,m);
Uz=zeros(N,m);

nur=zeros(1,m);
nu0=zeros(1,m);
nuz=zeros(1,m);

NUr=zeros(N,m);
NU0=zeros(N,m);
NUz=zeros(N,m);

yh=y.*h;

%The mode shape on the left of the plot corresponds to the highest mode for
%that frequency, as we move to the right of the plot we reach, on the right
%hand side the lowest mode for that frequency.
for i=1:1:m;

%Solid
Ux(1:N,i)= abs(PropMod((3*N+1):1:(4*N),i));
Uy(1:N,i)= abs(PropMod((4*N+1):1:(5*N),i));
Uz(1:N,i)= abs(PropMod((5*N+1):1:(6*N),i));


nux(1,i)=max(Ux(:,i));
nuy(1,i)=max(Uy(:,i));
nuz(1,i)=max(Uz(:,i));


NUx(:,i)= Ux(:,i)/nux(1,i);
NUy(:,i)= Uy(:,i)/nuy(1,i);
NUz(:,i)= Uz(:,i)/nuz(1,i);


figure(3)
%Solid
hold on
subplot(3,m,i)
plot( NUx(1:N,i), yh(:,1), 'k o','Markersize', 4)
hold on
grid on



figure(3)
%steel
hold on
subplot(3,m,(i+m))
plot( NUy(1:N,i), yh(:,1), 'k o','Markersize', 4)
hold on
grid on



figure(3)
%steel
hold on
subplot(3,m,(i+(2*m)))
plot( NUz(1:N,i), yh(:,1), 'k o','Markersize', 4)
hold on
grid on



end


fUx=zeros(N,m);
fUy=zeros(N,m);
fUz=zeros(N,m);

for i=1:1:m;

%Solid
fUx(1:N,i)= angle(PropMod(1:1:N,i));
fUy(1:N,i)= angle(PropMod((N+1):1:2*N,i));
fUz(1:N,i)= angle(PropMod((2*N+1):1:3*N,i));

CVphasex(1,i)=fUx( (((N-1)/2)+1) ,i); %centre value phase
CVphasey(1,i)=fUy( (((N-1)/2)+1) ,i); %centre value phase
CVphasez(1,i)=fUz( (((N-1)/2)+1) ,i);

FUx(1:N,i)=fUx(1:N,i)-CVphasex(1,i);
FUy(1:N,i)=fUy(1:N,i)-CVphasey(1,i);
FUz(1:N,i)=fUz(1:N,i)-CVphasez(1,i);

for j=1:1:N
  
    if FUx(j,i)<-pi;
        
       FUx(j,i)= 2*pi + FUx(j,i);
        
    elseif FUx(j,i)>pi;
    
       FUx(j,i)= FUx(j,i) - 2*pi;
       
    else
        
        FUx(j,i)= FUx(j,i);
        
    end 
    
    
    if FUy(j,i)<-pi;
        
       FUy(j,i)= 2*pi + FUy(j,i);
        
    elseif FUy(j,i)>pi;
    
       FUy(j,i)= FUy(j,i) - 2*pi;
       
    else
        
        FUy(j,i)= FUy(j,i);
        
    end
    
     if FUz(j,i)<-pi;
        
       FUz(j,i)= 2*pi + FUz(j,i);
        
    elseif FUz(j,i)>pi;
    
       FUz(j,i)= FUz(j,i) - 2*pi;
       
    else
        
        FUz(j,i)= FUz(j,i);
        
     end
    
end
    

end


%The mode shape on the left of the plot corresponds to the highest mode for
%that frequency, as we move to the right of the plot we reach, on the right
%hand side the lowest mode for that frequency.
for i=1:1:m;

figure(3)
%Solid
hold on
subplot(3,m,i)
plot( FUx(1:N,i), yh(:,1), '-r','linewidth', 2)
set(gca,'XTick',-3.1417:3.1417/2:3.1417)
xlabel('{\pi}' , 'FontWeight' , 'bold', 'FontSize' , 15)
set(gca,'XTickLabel',{'-1','-1/2','0','1/2','1'} )
set(gca,'YTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-3.1417 3.1417 0 yh(1,1)])
hold on


if i==3;
title({ 'b)', ' ' ,  'U_x',}, 'FontWeight' , 'bold', 'FontSize' , 12 )
end

figure(4)
%Solid
hold on
subplot(3,m,i)
plot( fUx(1:N,i), yh(:,1), 'black ')
hold on
axis([-3.1417 3.1417 0 yh(1,1)])
grid on

if i==m;
title('Ux displacement Profiles')
end



figure(3)
%steel
hold on
subplot(3,m,(i+m))
plot( FUy(1:N,i), yh(:,1), '-r','linewidth', 2)
set(gca,'XTick',-3.1417:3.1417/2:3.1417)
xlabel('{\pi}' , 'FontWeight' , 'bold', 'FontSize' , 15)
set(gca,'XTickLabel',{'-1','-1/2','0','1/2','1'} )
set(gca,'YTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-3.1417 3.1417 0 yh(1,1)])
hold on


if i==3;
title('U_y', 'FontWeight' , 'bold', 'FontSize' , 12 )
end

figure(4)
%steel
hold on
subplot(3,m,(i+m))
plot( fUy(1:N,i), yh(:,1), 'black ')
axis([-3.1417 3.1417 0 yh(1,1)])
hold on
grid on

if i==m;
title('Uy displacement Profiles')
end


figure(3)
%steel
hold on
subplot(3,m,(i+(2*m)))
plot( FUz(1:N,i), yh(:,1), '-r','linewidth', 2 )
set(gca,'XTick',-3.1417:3.1417/2:3.1417)
xlabel('{\pi}' , 'FontWeight' , 'bold', 'FontSize' , 15)
set(gca,'XTickLabel',{'-1','-1/2','0','1/2','1'} )
set(gca,'YTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-3.1417 3.1417 0 yh(1,1)])
hold on


if i==3;
title('U_z ', 'FontWeight' , 'bold', 'FontSize' , 12 )
end


figure(4)
%steel
hold on
subplot(3,m,(i+(2*m)))
plot( fUz(1:N,i), yh(:,1), 'black')
axis([-3.1417 3.1417 0 yh(1,1)])
hold on
grid on

if i==m;
title('Uz displacement Profiles')
end

end


%3D PLOTS

Z=yh;

for i=1:1:m; 
  
CFUx=cos(FUx);
SFUx=sin(FUx);

CFUy=cos(FUy);
SFUy=sin(FUy);

CFUz=cos(FUz);
SFUz=sin(FUz);

Xx(1:N,i)= NUx(1:N,i).*CFUx(1:N,i);
Yx(1:N,i)= NUx(1:N,i).*SFUx(1:N,i);

Xy(1:N,i)= NUy(1:N,i).*CFUy(1:N,i);
Yy(1:N,i)= NUy(1:N,i).*SFUy(1:N,i);

Xz(1:N,i)= NUz(1:N,i).*CFUz(1:N,i);
Yz(1:N,i)= NUz(1:N,i).*SFUz(1:N,i);

a=Xx(1:N,i);
b=Yx(1:N,i);

figure(8)
%Solid
subplot(3,m,i)
plot3( a, b, Z, '-b ','linewidth', 3 );
set(gca,'XTick',-1:1:1)
set(gca,'XTickLabel',{'-1','0','1'} )
set(gca,'YTick',-1:1:1)
set(gca,'YTickLabel',{'-1','0','1'} )
set(gca,'ZTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-1 1 -1 1 0 yh(1,1)])  
grid on
hold on

if i==3;
title({'c)'}, 'FontWeight' , 'bold', 'FontSize' , 12 )
xlabel('Re U_x' , 'FontWeight' , 'bold', 'FontSize' , 12)
ylabel('{Im U_x}' , 'FontWeight' , 'bold', 'FontSize' , 12)
end


x=Xy(1:N,i);
y=Yy(1:N,i);

figure(8)
%Solid
subplot(3,m,(i+m))
plot3( x, y, Z, '-b ','linewidth', 3 );
set(gca,'XTick',-1:1:1)
set(gca,'XTickLabel',{'-1','0','1'} )
set(gca,'YTick',-1:1:1)
set(gca,'YTickLabel',{'-1','0','1'} )
set(gca,'ZTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-1 1 -1 1 0 yh(1,1)])  
grid on
hold on

if i==3;
xlabel('Re U_y' , 'FontWeight' , 'bold', 'FontSize' , 12)
ylabel('{Im U_y}' , 'FontWeight' , 'bold', 'FontSize' , 12)
end

u=Xz(1:N,i);
v=Yz(1:N,i);


figure(8)
%steel
subplot(3,m,(i+(2*m)))
plot3( u, v, Z, '-b ','linewidth', 3 );
set(gca,'XTick',-1:1:1)
set(gca,'XTickLabel',{'-1','0','1'} )
set(gca,'YTick',-1:1:1)
set(gca,'YTickLabel',{'-1','0','1'} )
set(gca,'ZTick',0:0.025:yh(1,1))
set(gca, 'GridLineStyle',':');
axis([-1 1 -1 1 0 yh(1,1)])  
grid on
hold on

if i==3;
xlabel('Re U_z' , 'FontWeight' , 'bold', 'FontSize' , 12)
ylabel('{Im U_z}' , 'FontWeight' , 'bold', 'FontSize' , 12)
end



end




T=0:((2*pi)/1000):2*pi;

% T=0; %Select a value for the product time*frequency=T to see the profile at a certain time

Z(:,1)=yh(:,1);


i=1; %select mode to display, correspondign to those in figure 3. 1st mode (i=1) is the first from the left.
    
for j = 1:1:size(T,2);
        
        argx= T(1,j) - FUx(1:N,i);
        COSFUx=cos(argx);
        AA =  NUx(1:N,i).*COSFUx(1:N,1);
    
        argy= T(1,j) - FUy(1:N,i);
        COSFUy=cos(argy);
        BB =  NUy(1:N,i).*COSFUy(1:N,1);
        
        argz= T(1,j) - FUz(1:N,i);
        COSFUz=cos(argz);
        CC =  NUz(1:N,i).*COSFUz(1:N,1);
        
        
        figure(9)
        plot(AA,Z);
        grid on;
        axis([-1 1 0 yh(1,1)])
        title('Ux displacement Profiles')
      
        figure(10)
        plot(BB,Z);
        grid on;
        axis([-1 1 0 yh(1,1)])
        title('Uy displacement Profiles')
    
        figure(11)
        plot(CC,Z);
        grid on;
        axis([-1 1 0 yh(1,1)])
        title('Uz displacement Profiles')
        
        F(j) = getframe;      
       
end


 movie(F,1,100); 
