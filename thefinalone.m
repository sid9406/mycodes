j=input('write number of joints')
M=input('number of members');

cone=zeros(j);    %cone is the connectivity matrix
conn=[];          %member is connected between which joints

syms c d tempx tempy e f;
L=[];              %length of members
angl=[];           %angle with global X in counter clockwise
coor=[];      %takes the coordinates of joints

for n=1:j
    coor(n,1)=input('enter x coordinate');
    coor(n,2)=input('enter y coordinate');
end
   %input of coordinates are taken

for n=1:M
    conn(n,1)=input('enter smaller node');
    conn(n,2)=input('enter larger node');
end
    %input of how members are connected is taken as input
    
for n=1:M-1
    for m=1:M-n
        if conn(m,2)>conn(m+1,2)
            tempx=conn(m,1);
            tempy=conn(m,2);
            conn(m,1)=conn(m+1,1);
            conn(m,2)=conn(m+1,2);
            conn(m+1,1)=tempx;
            conn(m+1,2)=tempy;
        end
    end
end

for n=1:M-1
    for m=1:M-n
        if conn(m,1)>conn(m+1,1)
            tempx=conn(m,1);
            tempy=conn(m,2);
            conn(m,1)=conn(m+1,1);
            conn(m,2)=conn(m+1,2);
            conn(m+1,1)=tempx;
            conn(m+1,2)=tempy;
        end
    end
end
   %sorting of the conn is done
   
for n=1:M
    cone(conn(n,1),conn(n,2))=n;
end

   % cone matrix is filled
   
for n=1:M
    L(1,n)=((coor(conn(n,2),1)-coor(conn(n,1),1))^2+(coor(conn(n,2),2)-coor(conn(n,1),2))^2)^0.5;
    angl(n,1)=atan((coor(conn(n,2),2)-coor(conn(n,1),2))/(coor(conn(n,2),1)-coor(conn(n,1),1)));
end

%length and angles are calculated and stored 

syms ar coss sinn;
disp('1.fixed   2.roller  3.free 4.hinge 5.pin')
t=[];

% t denotes the joint condition

for n=1:j
    t(n,1)=input('type');
end


X=[];
Y=[];
Z=zeros(j);
a=1;
% X,Y,Z stores the numbering of degree of freedom


for n=1:j
    if t(n,1)==2    %i.e for roller
        X(n,1)=a;
        a=a+1;
        for i=1:M
            if conn(i,2)==n
                Z(conn(i,2),conn(i,1))=a;
                a=a+1;
            end
        end
        for i=1:M
            if conn(i,1)==n
                Z(conn(i,1),conn(i,2))=a;
                a=a+1;
            end
        end
    end
    
    if t(n,1)==3     % for free
        X(n,1)=a;
        a=a+1;
       Z(n,n)=a;
       a=a+1;
       Y(n,1)=a;
       a=a+1;
    end
    
    if t(n,1)==4      %i.e for hinge
        for i=1:M
            if conn(i,2)==n
                Z(conn(i,2),conn(i,1))=a;
                a=a+1;
            end
        end
        for i=1:M
            if conn(i,1)==n
                Z(conn(i,1),conn(i,2))=a;
                a=a+1;
            end
        end
    end
    
    if t(n,1)==5  %i.e  for pin
        X(n,1)=a;
        a=a+1;
        for i=1:M
            if conn(i,2)==n
                Z(conn(i,2),conn(i,1))=a;
                a=a+1;
            end
        end
        for i=1:M
            if conn(i,1)==n
                Z(conn(i,1),conn(i,2))=a;
                a=a+1;
            end
        end
        Y(n,1)=a;
        a=a+1;
    end 
end

 %numbering of unknown is done
 
 
u=a-1;   % u is the no.of unknown degree of freedom

for n=1:j
    if t(n,1)==1    % for fixed
        X(n,1)=a;
        a=a+1;
        Z(n,n)=a;
        a=a+1;
        Y(n,1)=a;
        a=a+1;
    end
    
    if t(n,1)==2    %for roller
        Y(n,1)=a;
        a=a+1;
    end
    
    if t(n,1)==4    %for hinge
        X(n,1)=a;
        a=a+1;
        Y(n,1)=a;
        a=a+1;
    end
end

%numbering of known is done

K=zeros(a-1);
%K is the stiffness matrix(GLOBAL)


xi=[-.90617984593 -.538469310 0 .5384693101 0.9061798459];
w=[.236926885 .47862867 .56888888 .47862867 .236926885];

%xi is the abscisa points on[-1,1]
%w is the weight

z=sym('z',[1 5]);
% z maps xi points from [-1,1] to [0,L]

for n=1:M   %loop for all elements
    
    syms E I area ent int choice;
    choice=input('enter choice 1.constant 2.varying');
    % choice 1 represents E,I,A are constant all over the member
    
    if choice==1
     area=input('area');
     ent=input('E');
     int=input('I');
     for i=1:5
        ar(1,i)=area;
        E(1,i)=ent;
        I(1,i)=int;
     end
    end
    
    if choice==2
         for i=1:5
        ar(1,i)=input('area');
        E(1,i)=input('enter E');
        I(1,i)=input('enter I');
         end
    end
       syms x;
       % x is the symbolic variable
       
        N=[(L(1,n)-x)/L(1,n) x/L(1,n)];
        B=diff(N,x);
       
    sum=0;
      
  for i=1:5
    z(1,i)=L(1,n)*(1+xi(1,i))/2;
    F=B.'*E(1,i)*B*ar(1,i);
    sum=sum+subs(F,x,z(1,i))*w(1,i);
  end
  
 k1=sum*L(1,n)/2;
 N=[1-3*x^2/L(1,n)^2+2*x^3/L(1,n)^3 x-2*x^2/L(1,n)+(x^3/L(1,n)^2) 3*x^2/L(1,n)^2-2*x^3/L(1,n)^3 -x^2/L(1,n)+x^3/L(1,n)^2];
 B=diff(N,x,2);
 
 sum=0;
 for i=1:5
     F=B.'*E(1,i)*B*I(1,i);
    sum=subs(F,x,z(1,i))*w(1,i)+sum;
 end
 
 k2=sum*L(1,n)/2;
 
 symbol=[k1(1,1) 0       0       k1(1,2)   0       0;
   0       k2(1,1) k2(1,2) 0         k2(1,3) k2(1,4);
   0       k2(2,1) k2(2,2) 0         k2(2,3) k2(2,4);
   k1(2,1)  0       0       k1(2,2)   0       0;
   0        k2(3,1) k2(3,2) 0         k2(3,3)  k2(3,4);
   0        k2(4,1) k2(4,2) 0         k2(4,3)  k2(4,4)];

   %symbol is the local stiffness matrix for the element
   
   
        coss=cos(angl(n,1));
        sinn=sin(angl(n,1));
        
        % coss and sinn are the cosines and sines of angle of element with
        % the GLOBAL: X
        
        
        T=[coss sinn 0 0 0 0;
           -sinn coss 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 coss sinn 0;
           0 0 0 -sinn coss 0;
           0 0 0 0 0 1];
       
       %T is the transformation matrix
       
       
        symbol=T.'*symbol*T;
        k=zeros(a-1);
        for A=1:j
            for b=1:j
                if cone(A,b)==n
                    c=A;
                    f=A;
                    d=b;
                    e=b;
                end
            end
        end
        if t(c,1)==1||t(c,1)==3
            d=c;
        end
        if t(e,1)==1||t(e,1)==3
            f=e;
        end
        G=[0 X(c,1) Y(c,1) Z(c,d) X(e,1) Y(e,1) Z(e,f);
           X(c,1) symbol(1,1) symbol(1,2) symbol(1,3) symbol(1,4) symbol(1,5) symbol(1,6);
           Y(c,1) symbol(2,1) symbol(2,2) symbol(2,3) symbol(2,4) symbol(2,5) symbol(2,6);
           Z(c,d) symbol(3,1) symbol(3,2) symbol(3,3) symbol(3,4) symbol(3,5) symbol(3,6);
           X(e,1) symbol(4,1) symbol(4,2) symbol(4,3) symbol(4,4) symbol(4,5) symbol(4,6);
           Y(e,1) symbol(5,1) symbol(5,2) symbol(5,3) symbol(5,4) symbol(5,5) symbol(5,6);
           Z(e,f) symbol(6,1) symbol(6,2) symbol(6,3) symbol(6,4) symbol(6,5) symbol(6,6)];
        for i=2:7
            for m=2:7
                k(G(i,1),G(1,m))=G(i,m);
            end
        end
        K=K+k;
end
    
    % GLOBAL stiffness matrix is calculated
    % now its time to input loading
    
Q=sym('q',[1 a-1]).';
D=sym('d',[1 a-1]).';

  %Q,D are the final loading and displacement vectors
  
q1=[];
qq=[];
% q1 is for distributed loading

for n=1:a-1
    q1(n,1)=0;
end

q2=[];  % q2 is for concentrated loading

R=[];
for n=1:M
    for A=1:j
            for b=1:j
                if cone(A,b)==n
                    c=A
                    f=A;
                    d=b;
                    e=b;
                end
            end
    end
    CC=c;
    DD=d;
        if t(c,1)==1||t(c,1)==3
            d=c;
        end
        if t(e,1)==1||t(e,1)==3
            f=e;
        end
    for i=1:a-1
        qq(i,1)=0;
    end
    sum=0;
    
    for o=1:5
       R(1,o)=input('enter w');     % R stores the distributed loading
      z(1,o)=L(1,n)*(1+xi(1,o))/2;
      sum=sum+((R(1,o)*(L(1,n)-z(1,o))^2*z(1,o))/L(1,n)^2)*w(1,o);
    end
    
    sum=sum*L(1,n)/2;
    qq(Z(c,d),1)=-1*sum;
    sum=0;
    for o=1:5
        sum=sum+((R(1,o)*z(1,o)^2*(L(1,n)-z(1,o)))/L(1,n)^2)*w(1,o);
    end
    sum=sum*L(1,n)/2;
    qq(Z(e,f),1)=sum;
    sum=0;
    for o=1:5
        sum=sum+R(1,o)*z(1,o)*w(1,o);
    end
    
    sum=sum*L(1,n)/2;
    
    sum=sum+qq(Z(c,d),1)+qq(Z(e,f),1);
    sum=sum/L(1,n);
    
    
    qq(Y(DD,1),1)=-1*sum*cos(angl(n,1));
    qq(X(DD,1),1)=sum*sin(angl(n,1));
    sum=0;
    for o=1:5
        sum=sum+R(1,o)*w(1,o);
    end
    sum=sum*L(1,n)/2;
    sum=sum-(qq(Y(DD,1),1)^2+qq(X(DD,1))^2)^0.5;
    
    qq(Y(CC,1),1)=-1*sum*cos(angl(n,1));
    qq(X(CC,1),1)=sum*sin(angl(n,1));
    q1=q1+qq
    
     sum1=0;
     sum2=0;
  for o=1:5              %for axial distributed loading
    R(1,o)=input('enter w');
    z(1,o)=L(1,n)*(1+xi(1,o))/2;
    sum1=sum1+(L(1,n)-z(1,o))*R(1,o)*w(1,o)/2;
    sum2=sum2+z(1,o)*R(1,o)*w(1,o)/2;
  end
 
  q1(Y(CC,1),1)=q1(Y(CC,1),1)+sum1*sin(angl(n,1))
  q1(X(CC,1),1)=q1(X(CC,1),1)+sum1*cos(angl(n,1))
  q1(Y(DD,1),1)=q1(Y(DD,1),1)+sum1*sin(angl(n,1))
  q1(X(DD,1),1)=q1(X(DD,1),1)+sum1*cos(angl(n,1))
      
end



    for n=1:u
    q2(n,1)=input('enter known loading at joints');  % enter known loading
    end
    
    for n=1:u
    q2(n,1)=q1(n,1)+q2(n,1);
    end
    
    for n=u+1:a-1
    D(n,1)=input('enter known displacements');
    end
    
    sub=zeros(u);
    for n=1:u
    for m=u+1:a-1
    sub(n,1)=sub(n,1)+K(n,m)*D(m,1);
    end
    end
    
    for n=1:u
    q2(n,1)=q2(n,1)-sub(n,1);
    end
    
    T=zeros(u);
    for n=1:u
    for m=1:u
    T(n,m)=K(n,m);
    end
    end
    
    d=[];
    d=pinv(T)*q2;
    
    for n=1:u
    D(n,1)=d(n,1);
    end
  
    %D is ready
    Q=K*D;
    
    for n=u+1:a-1
        Q(n,1)=Q(n,1)-q1(n,1);
    end
    
    % Q is ready
    
    for n=1:j
        % first printing reactions at node n
        
        disp('print Q for node')
        disp(n)
        vpa(Q(X(n,1)))
        for i=1:j
            if Z(n,i)~=0
        vpa(Q(Z(n,i)))
            end
        end
        vpa(Q(Y(n,1)))
        
        % now printing displacement at node n
        
        disp('print D')
        vpa(D(X(n,1)))
        for i=1:j
            if Z(n,i)~=0
        vpa(D(Z(n,i)))
            end
        end
        vpa(D(Y(n,1)))
    end
    
    
        
     
       
            

