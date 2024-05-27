clc
clear
mu=2
l1=0.8;
r1=[]
aoi=[]
p2=[];
aoibase=[]
disp('Correlated Arrivals')
   l2=1;
for p1=0.1:0.2:0.9
    k=1
for mu=0.9:0.1:1.5
D0=[-2,0;0,-0.5]
D1=[2*p1,2*(1-p1);0.5*(1-p1),0.5*p1]

%[D0,D1]=map1(a1,l1)
A0=mu*eye(size(D0))
A1=D0-mu*eye(size(D0))
A2=D1
g=A0+A1+A2;

x1=size(D0);

[G,R]=QBD_IS(A0,A1,A2)
B0=D0
pi=QBD_pi(B0,A2,R)
x=length(pi);

x2=x/x1(1);
for i=1:x2;
    P{i}=pi((i-1)*x1(1)+1:(i)*x1(1));
end
SS=x1(1);
e = ones(SS,1);
g(x1(1),1:x1(1))=e;
k1=eye(x1(1))
theta=inv(g)*k1(1:x1(1),x1(1))
%theta=transpose(stat(g))
sum(theta)
w=0;
for i=1:x2
w=w+P{i}*e;
end;

dr1=mu^2*(1-P{1}*e);
alpha=transpose(theta)*D1/l1;

nr1=P{1}*transpose(alpha)*D1*(eye(x1(1))-R-D0)^(-2)*(eye(x1(1))-R)^(-2)*e

aoi(k)=l1*(1/(mu*l1)+(transpose(nr1)*e)/mu^2+transpose(theta)*inv(-D0)*e/l1)
%correlation Coefficient
r1(k)=(1*(transpose(theta)*inv(-D0)*D1)*inv(-D0)*e/l1-1/l1^2)/(2*transpose(theta)*inv(-D0)*e/l1-1/l1^2)
rho=l1/mu
aoibase(k)=(1/mu)*(1+1/rho+rho^2/(1-rho))
k=k+1;
end
p2(l2)=plot((0.9:0.1:1.5),aoi)
hold on
l2=l2+1;

end
legend([p2(1),p2(2),p2(3),p2(4),p2(5)],'p=0.1','p=0.3','p=0.5','p=0.7','p=0.9')
xlabel('Sevice rate \mu')
ylabel('AAoI')
