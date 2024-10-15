%%% SCE-UA PROGRAM FOR URBAN STORAGE FUNCTION MODEL %%%
clc;
clear;
tic;

% Parameter Input for SCE-UA %
nopt=7;        % Number of variables to be optimized
nc=20;         % Number of Complex; 20 is recommended
ncp=2*nopt+1;  % Number of Complex Population
nscp=nopt+1;   % Number of Sub-Complex Population
ntp=nc*ncp;    % Number of Total Population
maxgene=50;    % Maximum Generation
fprintf('\nParamer Value for SCE-UA\n')
formatSpec='  nopt=%d  nc=%d  ncp=%d  nscp=%d  ntp=%d  maxgene=%d\n';
fprintf(formatSpec,nopt,nc,ncp,nscp,ntp,maxgene)
OptResult=zeros(maxgene,nopt+1);

%% Data for Urban Storage Function Model %% 
load rq36.dat;     % input file
rain=rq36(:,1);    % observeed basin-average rainfall
obsq=rq36(:,2);    % observed discharge
ndata=length(rain);% the number of data
I=0.0017;          % inflow from other basins
ET=0;              % evapotranspiration
OUT=0;             % outflow from the basin
RIEO = rain+I-ET-OUT; % total input 
qRmax=0.033; % maximum rainwater (RW) drainage out to the basin
Q0=obsq (1); % river discharge before the rain starts
qR0=0; % RW drainage out to the basin from a sewer system at time = 0 
fprintf('\nInput Values for USF Model\n')
formatSpec='  ndata=%d  I=%.5f  ET=%.5f  OUT=%.5f  qRmax=%.5f  qR0=%.5f\n';
fprintf(formatSpec,ndata,I,ET,OUT,qRmax,qR0)

%%---------------------------------%%
x=zeros(ntp,nopt);
xf=zeros(ntp,1);
cx=zeros(ncp,nopt);
cf=zeros(ncp,1);
scx=zeros(nscp,nopt);
scf=zeros(nscp,1);
xbest=zeros(maxgene,nopt);
xworst=zeros(maxgene,nopt);
xfbest=zeros(maxgene,1);
xfworst=zeros(maxgene,1);
xrange=zeros(maxgene,nopt);

bl=[ 10   100   0.001   0.1    0.1     0    0.1]; % Lower Boundary
bu=[500  5000   0.05    1.0    1.0    50    1.0]; % Upper Boundary
range=bu-bl;
fprintf('\nLower and Upper Boundaries for Paremeter Searching\n')
formatSpec='  bl=[%7.1f  %8.1f   %.4f   %.4f   %.4f  %6.1f   %.3f]\n';
fprintf(formatSpec,bl)
formatSpec='  bu=[%7.1f  %8.1f   %.4f   %.4f   %.4f  %6.1f   %.3f]\n';
fprintf(formatSpec,bu)
if range > 0
   else
   '*** ERROR *** At least one range is MINUS';
end 

rand('state',0) %Initial value setting for randum number generation 

%% -- START -- %%
nrun=1;
for irun=1:nrun
    fprintf ('\n SCE-UA Run Number --> %2g\n',irun);

%%% Initial Generation igene=1 %%% 
igene=1;
fprintf ('\n GeneNo-->%3g',igene);

R=rand(ntp,nopt);
% nspace=10^6;  % search space = 10*nspace
% R=round(R*nspace)/nspace;
for ii=1:nopt
    x(:,ii)=bl(ii)+range(ii).*R(:,ii);
end
% x(1,:)=1*[50  500  0.05  0.6  0.465  5  0.4]; %initially included x
for j=1:ntp
    [xf(j)]=FunUsfRkg(x(j,:),obsq,RIEO,qRmax,Q0,qR0);
end

xfx=[xf,x];
xfx=sortrows(xfx,1);
xf=xfx(:,1);
x=xfx(:,2:nopt+1);
xbest(igene,:)=x(1,:);
xfbest(igene)=xf(1);
xworst(igene,:)=x(ntp,:);
xfworst(igene)=xf(ntp);
xsort=sort(x);
xrange(igene,:)=xsort(ntp,:)-xsort(1,:);
formatSpec='  xf=%6.7f   x=%.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f';
fprintf(formatSpec,xf(1),x(1,:));
OptResult(igene,:)=xfx(1,:);

%%%%%% Generation Loop %%%%%%
for igene=2:maxgene % Generation Loop
    fprintf ('\n GeneNo-->%3g',igene);
    for inc=1:nc  % Complex Loop
        %% Assign Points into Complexes from Total Population %%
        k1=1:ncp;
        k2=inc+(k1-1)*nc;
        cx(k1,:)=x(k2,:);
        cf(k1)=xf(k2);
        nbeta=ncp; % recommended
        for ibeta=1:nbeta
%%%% Select Individuals into Sub-Complex (Extermination Room) %%%%
%%%% using trapezoidal probability %%%%  
            for inscp=1:nscp
                ihit=1;
                while ihit > 0
                   ran1val=rand;
                   prob=0;
                   for incp=1:ncp
                        prob=prob+2*(ncp+1-incp)/(ncp*(ncp+1));
                        if ran1val <= prob
                            j=incp;
                            break;
                        end
                   end % incp
                   ihit=0;
                   if inscp == 1
                      k3(1)=j; % store in k3
                   else
                      for jcheck=1:inscp-1
                          if j == k3(jcheck)
                             ihit=ihit+1;
                          end
                      end % jcheck
                      if ihit==0, k3(inscp)=j; end
                   end
                end % ihit  
                            
                scx(inscp,:)=cx(k3(inscp),:);
                scf(inscp)=cf(k3(inscp));
             end % inscp
 
%--- Sub-Complex Manipulation (exterminate the worst one) ---%             
             nalpha=1;  % recommended
             for ialpha=1:nalpha
             scfscx=[scf,scx];
             scfscx=sortrows(scfscx,1);
             scf=scfscx(:,1);
             scx=scfscx(:,2:nopt+1);
             scxbest=scx(1,:);
             scxworst=scx(nscp,:);
             scx1=scx(1:nscp-1,:);
             scx1mean=mean(scx1);% centriod computation
             scfworst=scf(nscp);
           % Reflection Step %
             scxnew=2*scx1mean-scxworst;
             if bl <= scxnew & scxnew <= bu
             else
                ran1=rand(1,nopt);
                scxnew=bl+range.*ran1;% mutation step
             end
             [scfnew]=FunUsfRkg(scxnew,obsq,RIEO,qRmax,Q0,qR0);
             
             if scfnew >= scfworst
                % Contraction Step %
                scxnew=(scx1mean + scxworst)/2;
                [scfnew]=FunUsfRkg(scxnew,obsq,RIEO,qRmax,Q0,qR0);
                
                if scfnew >= scfworst
                   ran1=rand(1,nopt);
                   scxnew=bl+range.*ran1;% mutation step
                  [scfnew]=FunUsfRkg(scxnew,obsq,RIEO,qRmax,Q0,qR0);
                  end
             end
             scx(nscp,:)=scxnew(1,:);
             scf(nscp)=scfnew;
          end % ialpha
 %--- Get out of Extermination Room ---%
          cx(k3,:)=scx(1:nscp,:);%replace scx into cx using location in k3
          cf(k3)=scf(1:nscp);
          cfcx=[cf,cx];
          cfcx=sortrows(cfcx,1);
          cf=cfcx(:,1);
          cx=cfcx(:,2:nopt+1);
        end % ibeta   
      x(k2,:)=cx(k1,:); % replace cx into x using location stored in k2
      xf(k2)=cf(k1);
    end % inc
    
    xfx=[xf,x];
    xfx=sortrows(xfx,1);
    xf=xfx(:,1);
    x=xfx(:,2:nopt+1);
    xbest(igene,:)=x(1,:);
    xfbest(igene)=xf(1);
    xworst(igene,:)=x(ntp,:);
    xfworst(igene)=xf(ntp);
    xsort=sort(x);
    xrange(igene,:)=xsort(ntp,:)-xsort(1,:);
formatSpec='  xf=%6.7f   x=%.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f';
fprintf(formatSpec,xf(1),x(1,:));
OptResult(igene,:)=xfx(1,:);

end % igene

%Hyeto & Hydro graph 
FunHydro(x(1,:),rain,obsq,RIEO,qRmax,Q0,qR0)
      
end % irun

fprintf('\n');
toc 