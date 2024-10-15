%% USF Model by RKG Method 
% Parameter Known -- 2*2 Matrix System 
function [RMSE]=FunUsfRkg(Para,obsq,RIEO,qRmax,Q0,qR0)

n=length(RIEO);
dt = 0.1;
nstep = n/dt;
int = 1/dt;

qsum = zeros(nstep+1,1); % qsum=Q+qR
qR = zeros(nstep+1,1);
Q = zeros(nstep+1,1);
qR(1)=qR0; % rainwater drainage from the basin at initial stage
Q(1)=Q0; % river discharge at initial stage

%% parameters of USF Model
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
X = [x1 ; x2];% used during the RKG calculation

%% Solution using RKG Method and description of qR
for k=1:nstep
    jiten = floor(((k-1) + int)/int);
    
    u1=dt*FunUsf(X,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X1=X+u1./2;
    u2=dt*FunUsf(X,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X2=X+((sqrt(2)-1)/2).*u1+((2-sqrt(2))/2).*u2;
    u3=dt*FunUsf(X,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X3=X-(sqrt(2)/2).*u2+((2+sqrt(2))/2).*u3;
    u4=dt*FunUsf(X,RIEO(jiten),k1,k2,k3,p1,p2,z);
    X=X+(u1+(2-sqrt(2)).*u2+(2+sqrt(2)).*u3+u4)./6;  
 
    if X(1) < 0
       RMSE=999999;
       return
       end
       
    qsum(k+1)=X(1)^(1/p2);
    if qsum(k+1) > 999999
        RMSE=999999;
        return
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
jj=int+1:int:nstep+1;
%% error function "Root Mean Square Error"
RMSE=sqrt(mean((obsq-Q(jj)).^2));
