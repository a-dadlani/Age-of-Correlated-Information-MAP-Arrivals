clc
clear
for i = 1 : 5
    for j = 1 : 3
        AA{i,j} = zeros(1,1);
    end;
end;

file = fopen('para1.txt','w+');
Sl=[47,44,77,48,42,75,54,48,83,53,48,79,56,53,82];
Su=[52,47,80,52,46,79,58,52,87,57,52,82,60,55,86];
 sl=[6,4,10,4,2,8,7,4,11,6,4,10,8,8,12];
 su=[8,6,12,8,6,12,11,8,15,10,8,11,12,10,15];

% for la1 = 3 : 1 :5
%     for la2 = 6 : 1 : 8
%         for mu = 0.2 : 0.1 : 0.4
%             for t1 = 2 : 1 : 5
                % for t2 = 0.5 : 0.1:1 


i1=5; i2=5;
j1=3; j2=3;
% S1=32; S2=36;
% s1=1; s2=5;
% for i = S1 : S2
% for j = s1 : s2
%  XX(i-S1+1,j-s1+1) = 99999;
%         end;
%     end;
z=14;

p=0.7; q=1-p; la2=0.45;N=10;

for ar=i1:i2
    fprintf('\n%d  \n',ar);
for le=j1:j2
    fprintf('\n'); xInc = 1;
    z=z+1;
   S1=Sl(z);
   S2=Su(z);
   s1=sl(z);
   s2=su(z);
for S = S1 : S2
%     fprintf('%d\t',S);
    yInc = 1;
for s = s1 : s2
 Q = S-s; 
 t1 =1;
 t2 =2;
 %t1=0.5,1,1.5
%t2=1.5,2,2.5:
    
la1=3.5;
[D0,D1,n] = map(ar,la1);

mu=0.4;
%mu=0.3,0.4,0.5
[T,a,m] = pht(le,mu);

nm=n*m;
Nnm=(N+1)*nm;
sNnm=(s+1)*Nnm;
SNnm=sNnm+Q*(N+1)*n;
Nn=(N+1)*n;

I1=eye(N+1,N+1);
I2=eye(n,n);
I3=eye(m,m);
e=ones(SNnm,1);
e1=ones(Nnm,1);
e2=ones(Nn,1);
e3=ones(m,1);
e4=ones(nm,1);
e5=ones(n,1);

T0=-T*e3;

I12=kron(I1,I2);
I23=kron(I2,I3);
t1I12=t1*I12;
t1I12a=kron(t1I12, a);
I1D1I3=kron(I1,kron(D1,I3));
I1qD1I3=kron(I1,kron(q*D1,I3));
pD1I3=kron(p*D1,I3);
D1a=kron(D1,a);
I1D1a=kron(I1, D1a);
I1I2T0=kron(I12,T0);

F2=zeros(SNnm, SNnm);
F0=zeros(SNnm, SNnm);
F1=zeros(SNnm, SNnm);
G0=zeros(SNnm, SNnm);
%%%%%%% Matrix F2 %%%%%%%

F2(sNnm+1:sNnm+Nn, s*Nnm+1:sNnm)=t1I12a;

for l=s+2:S
        F2(sNnm+(l-s-1)*Nn+1:sNnm+(l-s-1)*Nn+Nn, sNnm+(l-s-2)*Nn+1:sNnm+(l-s-2)*Nn+Nn)=t1I12;
end

%%%%%%%%%%%% Matrix F0 %%%%%%%%%

F0(1:Nnm, 1:Nnm)=I1D1I3;
for l=1:s
    F0(l*Nnm+1:l*Nnm+Nnm, l*Nnm+1:l*Nnm+Nnm)=I1qD1I3;
end

%%%%%%%%%%%%% Matrix F1 %%%%%%%%%%%
 for k=0:(N-1)
F1(k*nm+1:k*nm+nm, (k+1)*nm+1:(k+1)*nm+nm)=(N-k)*la2*I23;
end;
for k=0:N
    F1(k*nm+1:k*nm+nm,k*nm+1:k*nm+nm)=kron(D0-((N-k)*la2*I2), I3)+kron(I2,T);
end;
for l=1:s
    for k=0:N
    F1(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, l*Nnm+k*nm+1:l*Nnm+k*nm+nm)=kron(D0-((N-k)*la2*I2+k*t2*I2),I3)+kron(I2,T);
    end;
end;
for l=s+1:S
   for k=0:N
       F1(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n, sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n)=D0-(((N-k)*la2+k*t2+t1)*I2);
   end;
end;
for l=1:s
    for k=0:N
    F1(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, (l-1)*Nnm+k*nm+1:(l-1)*Nnm+k*nm+nm)=pD1I3+((N-k)*la2*I23);
    
    end;
end;
for k= 0:N
   F1(sNnm+k*n+1:sNnm+k*n+n, s*Nnm+k*nm+1:s*Nnm+k*nm+nm)=D1a+kron((N-k)*la2*I2,a); 
end;
for l= s+2:S
    for k=0:N
      F1(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n, sNnm+(l-s-2)*Nn+k*n+1:sNnm+(l-s-2)*Nn+k*n+n)=D1+(N-k)*la2*I2;
    end;
end;

for l=1:s
    for k=1:N
        F1(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, (l-1)*Nnm+(k-1)*nm+1:(l-1)*Nnm+(k-1)*nm+nm)=k*t2*I23;
    end;
end;
for k=1:N
        F1(sNnm+k*n+1:sNnm+k*n+n, s*Nnm+(k-1)*nm+1:s*Nnm+(k-1)*nm+nm)=kron(k*t2*I2,a);
end;
for l=s+2:S
    for k=1:N
        F1(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n, sNnm+(l-s-2)*Nn+(k-1)*n+1:sNnm+(l-s-2)*Nn+(k-1)*n+n )=k*t2*I2;
    end;
end;
for l=0:s
    F1(l*Nnm+1:l*Nnm+Nnm, sNnm+(l+Q-s-1)*Nn+1:sNnm+(l+Q-s-1)*Nn+Nn)=I1I2T0;
end;


%%%%%%%% Matrix G0 %%%%%%%%


 for k=0:(N-1)
G0(k*nm+1:k*nm+nm, (k+1)*nm+1:(k+1)*nm+nm)=(N-k)*la2*I23;
end;
for k=0:N
    G0(k*nm+1:k*nm+nm,k*nm+1:k*nm+nm)=kron(D0-((N-k)*la2*I2), I3)+kron(I2,T);
end;
for l=1:s
    for k=0:N
    G0(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, l*Nnm+k*nm+1:l*Nnm+k*nm+nm)=kron(D0-((N-k)*la2*I2+k*t2*I2),I3)+kron(I2,T);
    end;
end;
for l=s+1:S
   for k=0:N
       G0(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n,sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n)=D0-(((N-k)*la2+k*t2)*I2);
   end;
end;
for l=1:s
    for k=0:N
    G0(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, (l-1)*Nnm+k*nm+1:(l-1)*Nnm+k*nm+nm)=pD1I3+((N-k)*la2*I23);
    
    end;
end;
for k= 0:N
   G0(sNnm+k*n+1:sNnm+k*n+n, s*Nnm+k*nm+1:s*Nnm+k*nm+nm)=D1a+kron((N-k)*la2*I2,a); 
end;
for l= s+2:S
    for k=0:N
      G0(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n, sNnm+(l-s-2)*Nn+k*n+1:sNnm+(l-s-2)*Nn+k*n+n)=D1+(N-k)*la2*I2;
    end;
end;

for l=1:s
    for k=1:N
        G0(l*Nnm+k*nm+1:l*Nnm+k*nm+nm, (l-1)*Nnm+(k-1)*nm+1:(l-1)*Nnm+(k-1)*nm+nm)=k*t2*I23;
    end;
end;
for k=1:N
        G0(sNnm+k*n+1:sNnm+k*n+n, s*Nnm+(k-1)*nm+1:s*Nnm+(k-1)*nm+nm)=kron(k*t2*I2,a);
end;
for l=s+2:S
    for k=1:N
        G0(sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n, sNnm+(l-s-2)*Nn+(k-1)*n+1:sNnm+(l-s-2)*Nn+(k-1)*n+n )=k*t2*I2;
    end;
end;
for l=0:s
    G0(l*Nnm+1:l*Nnm+Nnm, sNnm+(l+Q-s-1)*Nn+1:sNnm+(l+Q-s-1)*Nn+Nn)=I1I2T0;
end;
rr=G0+F0;
rrr=rr*e;
F=F1+F2+F0;
sum=F*e;
sum1=0;
for i= 1:SNnm
    sum1=sum1+sum(i,1);
end;


[G,R,U]=QBD_IS(F2,F1,F0);
pi=QBD_pi(F2,G0,R);
x=length(pi);
x1=SNnm;
x2=x/x1;
for i=1:x2;
    P{i}=pi((i-1)*x1+1:(i)*x1);
end
SS=SNnm;
e = ones(SS,1);
w=0;
for i=1:x2
w=w+P{i}*e;
end;




%%%%%%% System performance Measures %%%%%%%


mi = 0;
parfor i = 1 : x2
    for j = 0 : s
        
        mi=mi+j*(P{i}(1, j*Nnm+1:j*Nnm+Nnm))*e1;
       
%             mi = mi+j*sum(P{i}(1, j*Nnm+1:j*Nnm+Nnm));
        
    end;
    for j = s+1 : S
        
            mi = mi+j*(P{i}(1,sNnm+(j-s-1)*Nn+1:sNnm+(j-s-1)*Nn+Nn))*e2;
             %mi = mi+j*sum((P{i}(1,sNnm+(j-s-1)*Nn+1:sNnm+(j-s-1)*Nn+Nn)));
        
    end;
end;


mr=0;
parfor i = 1 :x2
        
            mr = mr +P{i}(1,sNnm+1:sNnm+Nn)*I1D1a*e1;
        
end;
for i = 2 :x2
        
            mr = mr +P{i}(1,sNnm+1:sNnm+Nn)*t1I12*e2;
        
end;
for i=1 :x2
    for k=0:N
    mr = mr +P{i}(1, sNnm+k*n+1:sNnm+k*n+n)*(kron(k*t2*I2,a)+kron((N-k)*la2*I2,a))*e4;
    end;
end;

mo=0;
parfor i=1:x2
    mo=mo+(i-1)*P{i}*e;
end;


mp=0;
for i=1:x2
    for l=0:s
        for k=0:N
            mp=mp+(P{i}(1,(l*Nnm)+(k*nm)+1:(l*Nnm)+(k*nm)+nm))*(k*I23)*e4;
        end;
    end;
end;

for i=1:x2
    for l=s+1:S
        for k=0:N
            mp=mp+P{i}(1,sNnm+(l-s-1)*Nn+k*n+1:sNnm+(l-s-1)*Nn+k*n+n)*k*I2*e5;
        end;
    end;
end;

 ch =0.15; cs = 9; cw1 =5;cw2 =2;
%  ch = 0.105; cs = 10.5; cw =0.65;csh =3.5;
%ch=0.13,0.15,0.17;cs=9,10,11;cw1=4,5,6;cw2=0.5,2,3.5

tc = ch*mi+cs*mr+cw1*mo+cw2*mp;
%
%fprintf('\n %f   %f   %f  %f',mi,mr,mo,mp);
 fprintf('  %1.6f\t  ',tc);
AA{ar,le}(xInc,yInc) = tc;
%XX(S-S1+1,s-s1+1) = tc
 yInc = yInc + 1;

end; 
  
    xInc = xInc + 1;
fprintf('\n');
end;
    
    
% for i = i1 : i2
%     for j = j1 : j2
        Tmin = min(min(AA{ar,le}));
        [Xmin,Ymin] = find(AA{ar,le} == Tmin);
        fprintf(file,'%1.3f\t %1.2f\t %1.2f\t%1.2f\t %d\t %d\t %d \t %d\t %1.6f\r\n',ch,cs,cw1,cw2,ar,le,Xmin+S1-1,Ymin+s1-1,Tmin);
%     end;
    %fprintf('\n');
% end;
end;
end;

                %end;
%             end;
%         end;
%     end;
% end;
