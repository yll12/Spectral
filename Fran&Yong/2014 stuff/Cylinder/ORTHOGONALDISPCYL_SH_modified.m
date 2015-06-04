a=5; 
b=5.005;
rho=5300;  

C(1,1) = 23.9e9; 
C(1,2) = 10.4e9; 
C(1,3) = 5e9; 
C(2,2) = 24.7e9;  
C(2,3) = 5.2e9; 
C(3,3) = 13.5e9; 
C(4,4) = 6.5e9;   
C(5,5) = 6.6e9;   
C(6,6) = 7.6e9;

N=100;
    
Vt=sqrt(C(6,6)/rho);  

[x,D]=chebdif(N,2);
h=b-a;
r=(h*x+a+b)/2;
D1=(2/h)*D(:,:,1);
D2=(2/h)^2*D(:,:,2);

O=zeros(N);

k=0;

m=0;

Lp=C(6,6)*(diag(r.^2)*D2+diag(r)*D1-eye(N))-C(4,4)*diag(r.^2)*(k(1,(m+1)))^2;

L=Lp;

S=C(6,6)*(diag(r.^2)*D1-diag(r));

L(1,:)=S(1,:);
L(N,:)=S(N,:);

Mp=-rho*diag(r.^2)*eye(N); Mp(1,1)=0; Mp(N,N)=0;

M=Mp;

[P,E]=eig(L,M);
w=sort(real(sqrt(diag(E))));
