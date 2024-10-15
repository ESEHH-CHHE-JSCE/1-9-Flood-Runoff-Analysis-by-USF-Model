%% USF Model Solved by Runge-Kutta-Gill Method
% Parameter Known -- 2*2 Matrix System   2022.10.22
function FunHydro(Para,rain,obsq,RIEO,qRmax,Q0,qR0)

ndata=length(RIEO);
dt=0.1;
nstep = ndata/dt;
int = 1/dt;

qsum = zeros(nstep+1,1); % qsum=Q+qR
qR = zeros(nstep+1,1);
Q = zeros(nstep+1,1);
qR(1)=qR0; % rainwater drainage out to the basin at initial stage
Q(1)=Q0;   % river discharge at initial stage

k1 = Para(1);
k2 = Para(2);
k3 = Para(3);
p1 = Para(4);
p2 = Para(5);
z  = Para(6);
alpha = Para(7);

qsum(1) = Q0+qR0; % Total discharge at time = 0 min
x1 = qsum (1)^p2;
x2 = 0;
X = [x1 ; x2]; % used during the RKG calculation
Xmat=zeros(ndata+1,2);
Xmat(1,:)=X(:)';

%% Solution using RKG Method and description of qR
for k=1:nstep
    jiten = floor(((k-1) + int)/int);
    
    u1=dt*FunUsf(X,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X1=X+u1./2;
    u2=dt*FunUsf(X1,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X2=X+((sqrt(2)-1)/2).*u1+((2-sqrt(2))/2).*u2;
    u3=dt*FunUsf(X2,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X3=X-(sqrt(2)/2).*u2+((2+sqrt(2))/2).*u3;
    u4=dt*FunUsf(X3,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X=X+(u1+(2-sqrt(2)).*u2+(2+sqrt(2)).*u3+u4)./6;  
    
    qsum(k+1)=X(1)^(1/p2);
    if k==1; X;
       end
    if mod(k,int)==0  
       Xmat(k/int+1,:)=X(:)';
       end
    if alpha*(qsum(k+1)-Q0)<qRmax;
         qR(k+1)=alpha*(qsum(k+1)-Q0);
       else
         qR(k+1)=qRmax;
       end
       if qR(k+1)<0;qR(k+1); 
          qR(k)=0;end
       if qR(k+1)>qsum(k+1); qR(k+1);
          qR(k+1)=qsum(k+1);end
       
end

Q=qsum-qR;

j=1:ndata+1;
tX = [(0:ndata)' Xmat(j,:)];

jj=int+1:int:nstep+1;
trq=[(1:ndata)' rain Q(jj)];

%% error function "Root Mean Square Error"
RMSE=sqrt(mean((obsq-trq(:,3)).^2));
%% Nash–Sutcliffe model efficiency coefficient (NSE)
NSE=1-sum((obsq-trq(:,3)).^2)/sum((obsq-mean(obsq)).^2);
fprintf('\nErro Evaluation Value\n')
fprintf('  RMSE=%.6f  NSE=%.4f\n',RMSE,NSE)

%% graphs
% Rainfall
subplot(2,1,1)
bar(0.5:ndata,rain,1)
axis ij
xlim([0 ndata])
grid on
ylabel('rainfall(mm/min)');
xlabel('time(min)');
% Discharge
subplot(2,1,2)
t=0:dt:ndata;
plot(t(jj),trq(:,3),'r-')
xlim([0 ndata])
grid on
ylabel('discharge(mm/min)');
xlabel('time(min)');
hold on
plot(1:ndata,obsq,'k','LineWidth',1)
hold off
title('Urban Storage Function model solved by R-K-G method');
legend('comq','obsq');

%% OUTPUT OF rainfall obsq calq
out=[(1:ndata)',rain,obsq,Q(jj)];
fid=fopen('hydro.out','w');
fprintf(fid,'Urban Storage Function Model Solved by R-K-G Method\n\n');
fprintf(fid,'   k       rain       obsq       calq\n');
fprintf(fid,'   0                       %8.4f\n',Q(1));
fprintf(fid,'%4.0f %10.9f %10.9f %10.9f\n',out');
fclose(fid);