%Assignment_2
%Sunanda and Abhijeet 
%-------- and 140906838

clc;
clear;
num=[0 0.01];
den=[0.005 .06 .1001 ];
sys=tf(num,den);
[a,b,c,d]=tf2ss(num,den);
A=fliplr(flipud(a));
B=flipud(b);
C=fliplr(c);

%Check Controllability 
if rank(ctrb(A,B))==rank(A)
    fprintf('The System is Controllable')
end

%Finding poles From Design Specifications
pos=5;
z=(-log(pos/100)/sqrt(pi^2+[log(pos/100)]^2));
Tsettling=2;
wn=4/(Tsettling*z);
wd=wn*sqrt(1+z^2);
poles=[-z*wn+wd*i -z*wn-wd*i];

%State Feedback Controller
aa=acker(A,B,poles);
An=A-B*aa;
Controller_System=ss(An,B,C,[0])
figure(1) 
step(An,B,C,[0])

% state feedback controller with initial conditions 
SYSTEM=ss(An,eye(2),eye(2),eye(2)); 
t=0:0.01:4; 
x=initial(SYSTEM,[1 0],t);
x1=[1 0]*x';
x2=[0 1]*x';
figure(2) 
subplot(2,1,1);
plot(t,x1),grid
subplot(2,1,2);
plot(t,x2),grid

%state feedback controller with integrator 
ahat=[A zeros(2,1);-C 0];
bhat=[B;0]; 
jhat=[poles -z*wn*5];
khat=acker(ahat,bhat,jhat) 
figure(3) 
k1=[khat(1) khat(2)]; 
k0=khat(3); 
aa1=[A-B*k1 B*k0;-C 0]; 
bb1=[0;0;1];
cc1=[c 0]; 
dd1=[0]; 
figure(3) 
step(aa1,bb1,cc1,dd1)

%Check For observability
if rank(obsv(A,C))==rank(A)
    fprintf('The System is Observable')
end

%Finding poles From Design Specifications
%pos=5;
%z=(-log(pos/100)/sqrt(pi^2+[log(pos/100)]^2));
Tsettling=2;
%five times faster
WN=4/(Tsettling*z*5);
WD=WN*sqrt(1+z^2);
Poles=[-z*WN+WD*i -z*WN-WD*i];

%Observer
L=acker(A',C',Poles')'
AN=A-L*C;
OBSE_System=ss(AN,B,C,[0])
figure(3) 
step(AN,B,C,[0])

%LQR
LQR=lqr(A,B,eye(2),1)

