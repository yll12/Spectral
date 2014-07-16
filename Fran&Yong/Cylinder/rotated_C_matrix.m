%computes the rotated C matrix from the stiffness matrix c of the material

function [C] = rotated_C_matrix(alpha, c)

a=alpha

if a==pi/2
    psi=0
    phi=1
else
psi=cos(alpha)
phi=sin(alpha)
end

C=zeros(6);

C(1,1)=psi^2*(c(1,1)*psi^2+c(1,3)*phi^2-2*c(1,5)*phi*psi)+phi^2*(c(1,3)*psi^2+c(3,3)*phi^2-2*c(3,5)*phi*psi)-2*phi*psi*(c(1,5)*psi^2+c(3,5)*phi^2-2*c(5,5)*phi*psi);
C(1,2)=c(1,2)*psi^2+c(2,3)*phi^2-2*c(2,5)*phi*psi; 
C(1,3)=psi^2*(c(1,1)*phi^2+c(1,3)*psi^2+2*c(1,5)*phi*psi)+phi^2*(c(1,3)*phi^2+c(3,3)*psi^2+2*c(3,5)*phi*psi)-2*phi*psi*(c(1,5)*phi^2+c(3,5)*psi^2+2*c(5,5)*phi*psi);  
C(1,4)=psi^2*(c(1,4)*psi+c(1,6)*phi)+phi^2*(c(3,4)*psi+c(3,6)*phi)-2*phi*psi*(c(4,5)*psi+c(5,6)*phi);
C(1,5)=psi^2*(c(1,1)*phi*psi-c(1,3)*phi*psi+c(1,5)*(-phi^2+psi^2))+phi^2*(c(1,3)*phi*psi-c(3,3)*phi*psi+c(3,5)*(-phi^2+psi^2))-2*phi*psi*(c(1,5)*phi*psi-c(3,5)*phi*psi+c(5,5)*(-phi^2+psi^2)); 
C(1,6)=psi^2*(-c(1,4)*phi+c(1,6)*psi)+phi^2*(-c(3,4)*phi+c(3,6)*psi)-2*phi*psi*(-c(4,5)*phi+c(5,6)*psi); 

C(2,2)=c(2,2);  
C(2,3)=c(1,2)*phi^2+c(2,3)*psi^2+2*c(2,5)*phi*psi;
C(2,4)=c(2,4)*psi+c(2,6)*phi;
C(2,5)=c(1,2)*phi*psi-c(2,3)*phi*psi+c(2,5)*(-phi^2+psi^2);
C(2,6)=-c(2,4)*phi+c(2,6)*psi;

C(3,3)=phi^2*(c(1,1)*phi^2+c(1,3)*psi^2+2*c(1,5)*phi*psi)+psi^2*(c(1,3)*phi^2+c(3,3)*psi^2+2*c(3,5)*phi*psi)+2*phi*psi*(c(1,5)*phi^2+c(3,5)*psi^2+2*c(5,5)*phi*psi);
C(3,4)=phi^2*(c(1,4)*psi+c(1,6)*phi)+psi^2*(c(3,4)*psi+c(3,6)*phi)+2*phi*psi*(c(4,5)*psi+c(5,6)*phi);
C(3,5)=phi^2*(c(1,1)*phi*psi-c(1,3)*phi*psi+c(1,5)*(-phi^2+psi^2))+psi^2*(c(1,3)*phi*psi-c(3,3)*phi*psi+c(3,5)*(-phi^2+psi^2))+2*phi*psi*(c(1,5)*phi*psi-c(3,5)*phi*psi+c(5,5)*(-phi^2+psi^2));
C(3,6)=phi^2*(-c(1,4)*phi+c(1,6)*psi)+psi^2*(-c(3,4)*phi+c(3,6)*psi)+2*phi*psi*(-c(4,5)*phi+c(5,6)*psi);

C(4,4)=psi*(c(4,4)*psi+c(4,6)*phi)+phi*(c(4,6)*psi+c(6,6)*phi);
C(4,5)=psi*(c(1,4)*phi*psi-c(3,4)*phi*psi+c(4,5)*(-phi^2+psi^2))+phi*(c(1,6)*phi*psi-c(3,6)*phi*psi+c(5,6)*(-phi^2+psi^2));
C(4,6)=psi*(-c(4,4)*phi+c(4,6)*psi)+phi*(-c(4,6)*phi+c(6,6)*psi);

C(5,5)=phi*psi*(c(1,1)*phi*psi-c(1,3)*phi*psi+c(1,5)*(-phi^2+psi^2))-phi*psi*(c(1,3)*phi*psi-c(3,3)*phi*psi+c(3,5)*(-phi^2+psi^2))+(-phi^2+psi^2)*(c(1,5)*phi*psi-c(3,5)*phi*psi+c(5,5)*(-phi^2+psi^2));
C(5,6)=phi*psi*(-c(1,4)*phi+c(1,6)*psi)-phi*psi*(-c(3,4)*phi+c(3,6)*psi)+(-phi^2+psi^2)*(-c(4,5)*phi+c(5,6)*psi);

C(6,6)=-phi*(-c(4,4)*phi+c(4,6)*psi)+psi*(-c(4,6)*phi+c(6,6)*psi);

C;
a;






