%matlabpool local 2;            % 
clear all;
%-------------------------%
S0  = 36;                 %
rf  = 0.015;               %
sigma=0.2;                %
N  =52;                   %Number of points in time grid to use (minimum is 3, default is 50)
M  =10000;                %Number of points in asset price grid to use (minimum is 3, default is 50)
dt = 1/N;
%----------------------------%
L0 = 100 ;                %
g  = 0;                %
T  = 30;                  %
MaxT=30;                  %
TotalTime = N*MaxT;       %
ActualTime= N*T;        %
SimulateTimes = M ;       %
%---------%
TimeCell = ActualTime/13; %
PeriodCell = 13;        %
FirstCell = 13;         %
Celldt = 13/52;         %
%---------%
s=1;
%----------%
rng(s);
WSimulate1=normrnd(0,1,M/2,TotalTime);
WSimulate2=-1*WSimulate1;
TotalWSimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
ActualWSimulate=TotalWSimulate(1:SimulateTimes,1:ActualTime);
clear('TotalWSimulate');
%-----------%
Szero=ones(M,1);
S=zeros (M,N);
Sn=zeros(M,1);
Ss=zeros(M,N+1);
Sn(:,1)=S0.*Szero(:,1);
Ss=[Sn,S];                      %
%---------------%   
for c = 1 : 1 : ActualTime
    Ss(:,c+1)=Ss(:,c).*(exp((rf-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*ActualWSimulate(:,c)));
end
clear ('ActualWSimulate','Szeros','S','Sn');
%--------------%
TrueCashFlow=zeros(M,ActualTime);
Y1=zeros(M,1,'single');
Xx=zeros(M,1,'single');
X1=zeros(M,1,'single');
R1=zeros(4,M,'single');
a2=zeros(4,M,'single');
ContinValue=zeros(M,1,'single');
TrueCashFlow(:,end)=max(L0*(1+g)^T,L0*(Ss(:,end)./S0)); 
for b =  (ActualTime-PeriodCell) : -PeriodCell : FirstCell

    Idx = find((L0*(1+g)^T)*exp(-rf*(ActualTime-b)/PeriodCell*Celldt) < max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*(Ss(:,b+1)./S0)));%
    Idx=uint32(Idx);
    if b == (ActualTime-PeriodCell)
        Y1=TrueCashFlow(Idx,(b+PeriodCell)).*exp(-rf*Celldt); %
        Xx=Ss(Idx,b+1);          %
        X1=Xx;                % Use Price as variable
        %------------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1');
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;
        clear ('R1','a1','a2','Y1');
        %------------------%
    else                         % 
        Y1=zeros(size(Idx));
        for d = (b+PeriodCell) : PeriodCell : ActualTime
            [ExTime,~,~]=find(TrueCashFlow(Idx,d)~=0);                                    %[rowIndex,colIndex,value]=find(x)
            Y1(ExTime,1)=TrueCashFlow(Idx(ExTime),d).*exp(-rf*((d-b)/PeriodCell)*Celldt); %Extime  , celldt=13/52;
        end
        Xx=Ss(Idx,b+1);   %
        X1=Xx;         %
        %------------------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1'); 
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;       
        clear ('R1','a1','a2','Y1');
        %----------------------%
    end 
    Jdx=find(ContinValue>max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*(Ss(Idx,b+1)./S0)));   %
    clear ContinValue;
    TrueCashFlow(Idx(Jdx),b) = max(L0*(1+g)^((b/PeriodCell)*Celldt),L0.*(Ss(Idx(Jdx),b+1)./S0));   % 
    TrueCashFlow(Idx(Jdx),(b+PeriodCell):ActualTime)=0;  
    clear ('Idx','Jdx');
end 
    TotalPrice=0;
    nn=0; %check Extime Distribue
for DateCount = FirstCell : PeriodCell : ActualTime
    [RowFinalExTime,~,~] = find( TrueCashFlow(:,DateCount) ~= 0 ); %
    if DateCount~=ActualTime
        nn=size(RowFinalExTime,1)+nn;
    end
    PeriodDebtValue=max(L0*(1+g)^(DateCount/N)- L0.*(Ss(RowFinalExTime,DateCount+1)./S0),0);
    TempPrice1=sum(PeriodDebtValue.*exp(-rf*((DateCount/PeriodCell)*Celldt))); %%Celldt=13/52
    TotalPrice=TempPrice1+TotalPrice; %
end
%-----Euro-----%
EuroExPrice=max(L0*(1+g)^T-L0*(Ss(:,end)./S0),0).*exp(-rf*T);
EuroFinalPrice=sum(EuroExPrice)/M;
%------%
AmerFinalPrice=TotalPrice/M; %