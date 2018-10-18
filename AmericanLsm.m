%matlabpool local 2; % use two cores 
clear all;
S0 = 36;      %
K  = 40;      %
rf  = 0.06;    %
T  = 1;       %Time to maturity of option
sigma=0.2;     %
N  =50;       %Number of points in time grid to use (default is 50)
M  =50000;    %Number of points in asset price grid to use (default is 50)
dt = T/N;
%---------%
s=1;
rng(s);
WSimulate1=normrnd(0,1,M/2,N);
WSimulate2=-1*WSimulate1;
WSimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
S=zeros(M,N,'single');
Sn=zeros(M,1,'single');
Sn(:,1)=S0;
Ss=[Sn,S]; %
clear ('S','Sn');
TrueCashFlow=zeros(M,N,'single');
 %
for c = 1 : 1 : N
    Ss(:,c+1)=Ss(:,c).*(exp((rf-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*WSimulate(:,c)));
end
clear WSimulate;
TrueCashFlow(:,end)=max(K-Ss(:,end),0); 
for b =  (N-1) : -1 : 1
    Y1=zeros(M,1,'single');
    Xx=zeros(M,1,'single');
    X1=zeros(M,1,'single');
    R1=zeros(4,M,'single');
    a2=zeros(4,M,'single');
    ContinValue=zeros(M,1,'single');
    Idx = find(Ss(:,b+1)<K);%
    
    if b == (N-1)
        Y1=TrueCashFlow(Idx,(b+1)).*exp(-rf*dt); %Y1a
        Xx=Ss(Idx,b+1); %
        X1=Xx/S0; %
        %
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1');
        a1 = (R1')*R1;
        a2 = a1\(R1');
        clear a1;
        TempValue=R1*a2;
        clear ('a2','R1');
        ContinValue= TempValue*Y1;
        clear ('Y1','TempValue');
        Jdx = find((K-Ss(Idx,b+1) ) > ContinValue);%   
        clear('ContinValue');
        TrueCashFlow(Idx(Jdx),b) = K-Ss(Idx(Jdx),b+1); % 
        TrueCashFlow(Idx(Jdx),(b+1):N)=0; %
        clear ('Idx','Jdx');
    else % 
        Y1=zeros(size(Idx),'single');
        for d = (b+1) : 1 : N
            [ExTime,~,~]=find(TrueCashFlow(Idx,d)~=0);   %Extime
            Y1(ExTime,1)=TrueCashFlow(Idx(ExTime),d).*exp(-rf*(d-b)*dt); %
            clear ExTime;
        end
        Xx=Ss(Idx,b+1); %
        X1=Xx/S0; %
        %
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1');
        a1 = (R1')*R1;
        a2 = a1\(R1');
        clear a1;
        TempValue=R1*a2;
        clear ('a2','R1');
        ContinValue= TempValue*Y1;
        clear ('Y1','TempValue');
        Jdx = find((K-Ss(Idx,b+1) ) > ContinValue);%    
        clear('ContinValue');
        TrueCashFlow(Idx(Jdx),b) = K-Ss(Idx(Jdx),b+1); % 
        TrueCashFlow(Idx(Jdx),(b+1):N)=0; %
        clear ('Idx','Jdx');
    end
end 
clear Ss;
    TotalPrice=0;
for e = 1 : 1 : N
    [RowFinalExTime,~,~] = find( TrueCashFlow(:,e) ~= 0 ); %
    TempPrice1=sum(TrueCashFlow(RowFinalExTime,e).*exp(-rf*e*dt)); %
    clear RowFinalExTime;
    TotalPrice=TempPrice1+TotalPrice; %
end
FinalPrice=TotalPrice/M; %