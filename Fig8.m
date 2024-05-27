clc
clear
mu=0.6
r1=[]
u1=[];
aoi=[]
aoibase=[]
peakaoi=[]
mp=[];
disp('Arrival')
a1=3
    for l2=3:7
          u1(1:l2)=1/l2;
          k=1
    for l1=0.1:0.02:0.55
        [D0,D1]=map1(a1,l1)
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
    l11=l1*0.5/l2;
    alpha=transpose(theta)*0.5*(l2)*D1/l11;

    nr1=P{1}*transpose(alpha)*0.5*D1*(eye(x1(1))-R-D0)^(-2)*(eye(x1(1))-R)^(-2)*e

    aoi(k)=l11*(1/(mu*l11)+(transpose(nr1)*e)/mu^2+transpose(theta)*inv(-D0)*e/l11)
    %correlation Coefficient
    r1(k)=(1*(transpose(theta)*inv(-D0)*0.5*(l2)*D1)*inv(-D0)*e/l11-1/l11^2)/(2*transpose(theta)*inv(-D0)*e/l11-1/l11^2)
    rho=l1/mu
    aoibase(k)=(1/mu)*(1+1/rho+rho^2/(1-rho))
    peakaoi(k)=(1/l11)+(P{1}*R*inv(eye(x1(1))-R)^(-1)*e)/l11;
    k=k+1;
  end
   x3=0.1:0.02:0.55
   mp(l2)=plot(x3,peakaoi)
   hold on
    end
legend([mp(7),mp(6),mp(5),mp(4),mp(3)],'N=7','N=6','N=5','N=4','N=3')
xlabel('Arrival rate \lambda_1')
ylabel('Average Peak AoI for source1')

 %legend(mp,'N=7','N=6','N=5','N=4,'N=3')