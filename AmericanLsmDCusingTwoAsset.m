%matlabpool local 2;            % 
clear all;
%-------------------------%
S0A= 36;                
S0B= 36;
rf  = 0.015;             
sigmaA=0.2;           %0.13;                
sigmaB=0.2;           %0.21;
N  =52;                   %Number of points in time grid to use (minimum is 3, default is 50)
M  =10000;                %Number of points in asset price grid to use (minimum is 3, default is 50)
dt = 1/N;
%----------------------------%
BetaA =0; %0.75;
BetaB =0; %1.28;
RiskPremium = 0.0475; 
uA = rf+(BetaA*RiskPremium);  %
uB = rf+(BetaB*RiskPremium) ;  %
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
%-----------%
corrAB=0;
s=1;
v=2;   
rng(s);
%----------------------------%
WSimulate1=normrnd(0,1,M/2,TotalTime);
WSimulate2=-1*WSimulate1;
WASimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
ActualWASimulate=WASimulate(1:SimulateTimes,1:ActualTime);
%----------------------------%
rng(v);
WSimulate1=normrnd(0,1,M/2,TotalTime);
WSimulate2=-1*WSimulate1;
WBSimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
alpha21=corrAB;     	%a21*a11=1
alpha22=(1-alpha21^2)^(1/2);  %a21^2+a22^2=1
WBs=(alpha21.*WASimulate)+(alpha22.*WBSimulate);
ActualWBSimulate=WBs(1:SimulateTimes,1:ActualTime);
YhatA = zeros(SimulateTimes,ActualTime,'single');
YhatB = zeros(SimulateTimes,ActualTime,'single');
SzeroA=ones(SimulateTimes,1);
SzeroB=ones(SimulateTimes,1);
SA=zeros(SimulateTimes,ActualTime);
SB=zeros(SimulateTimes,ActualTime);
SnA=zeros(SimulateTimes,1);
SnB=zeros(SimulateTimes,1);
SAa=zeros(SimulateTimes,ActualTime+1);
SBb=zeros(SimulateTimes,ActualTime+1);
RAa=zeros(SimulateTimes,ActualTime+1);
RBb=zeros(SimulateTimes,ActualTime+1);
SnA(:,1)=S0A.*SzeroA(:,1);
SnB(:,1)=S0B.*SzeroB(:,1);
RnA(:,1)=S0A.*SzeroA(:,1);
RnB(:,1)=S0B.*SzeroB(:,1);
SAa=[SnA,SA]; 
SBb=[SnB,SB]; 
RAa=[RnA,SA]; 
RBb=[RnB,SB];   
Y1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
X1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);     
Y1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
X1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
ReturnA=zeros(SimulateTimes,TimeCell);
ReturnB=zeros(SimulateTimes,TimeCell);
PeriodReturnA=zeros(SimulateTimes,TimeCell);
PeriodReturnB=zeros(SimulateTimes,TimeCell);
TotalswitchAssetA=0;
TotalswitchAssetB=0;
TempswitchAssetA=0;
TempswitchAssetB=0;
PeriodswitchAssetA=zeros(T,1);
PeriodswitchAssetB=zeros(T,1);
for c = 1:1:ActualTime
        SAa(:,c+1)=SAa(:,c).*(exp((uA-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*ActualWASimulate(:,c)));
        SBb(:,c+1)=SBb(:,c).*(exp((uB-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*ActualWBSimulate(:,c))); 
        %--Risk Neutual Price
        RAa(:,c+1)=RAa(:,c).*(exp((rf-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*ActualWASimulate(:,c))); 
        RBb(:,c+1)=RBb(:,c).*(exp((rf-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*ActualWBSimulate(:,c)));
end
for b = FirstCell : PeriodCell : ActualTime-PeriodCell 
        Y1a(:,b/PeriodCell)=SAa(:,b+PeriodCell+1);
        Y1aa=gpuArray(Y1a(:,b/PeriodCell));
        X1a(:,b/PeriodCell)=SAa(:,b+1);
        R1a =[ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-X1a(:,b/PeriodCell))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-(2.*X1a(:,b/PeriodCell)))+(X1a(:,b/PeriodCell).^2)]; 
        R1aa=gpuArray(R1a);
        clear X1a;
        a1 = (R1a')*R1a;
        a2 = a1\(R1a');
        YhatA = R1aa*a2*Y1aa;
        clear ('R1a','a2','a1','Y1a'); 
        Y1b(:,b/PeriodCell)=SBb(:,b+PeriodCell+1);
        Y1bb=gpuArray(Y1b(:,b/PeriodCell));
        X1b(:,b/PeriodCell)=SBb(:,b+1);    
        R1b =[ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-X1b(:,b/PeriodCell))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-(2.*X1b(:,b/PeriodCell)))+(X1b(:,b/PeriodCell).^2)];
        R1bb=gpuArray(R1b);
        clear X1b;
        b1 = (R1b')*R1b;
        b2 = b1\(R1b');
        YhatB = R1bb*b2*Y1bb;
        clear('R1b','b2','b1','Y1b');
        YhatAa=gather(YhatA);
        YhatBb=gather(YhatB);
        EA(:,b/PeriodCell)=log(YhatAa./SAa(:,b+1));
        EB(:,b/PeriodCell)=log(YhatBb./SBb(:,b+1));
end
ReturnA( : ,FirstCell/PeriodCell)  =    log(RAa(:,FirstCell+1)./S0A) ;
ReturnB( : ,FirstCell/PeriodCell)  =    log(RBb(:,FirstCell+1)./S0B) ;
for d = FirstCell : PeriodCell : ActualTime-PeriodCell
        [SwitchToA,~,~]=find(EA(:,d/PeriodCell)>EB(:,d/PeriodCell));
        ReturnA(SwitchToA,(d+PeriodCell)/PeriodCell)    =   log(RAa(SwitchToA,d+PeriodCell+1)./RAa(SwitchToA,d+1));
        ReturnB(SwitchToA,(d+PeriodCell)/PeriodCell)    =   log(RAa(SwitchToA,d+PeriodCell+1)./RAa(SwitchToA,d+1));
        [ForntAssetIsB,~,~]=find(ReturnA(SwitchToA,d/PeriodCell)==log(RBb(SwitchToA,d+1)./RBb(SwitchToA,d-PeriodCell+1)));       
        TotalswitchAssetA=size(ForntAssetIsB,1)+TotalswitchAssetA;       
        TempswitchAssetA=size(ForntAssetIsB,1)+TempswitchAssetA;
        [ForntAssetIsB,~,~]=find(ReturnB(SwitchToA,d/PeriodCell)==log(RBb(SwitchToA,d+1)./RBb(SwitchToA,d-PeriodCell+1)));
        TotalswitchAssetB=size(ForntAssetIsB,1)+TotalswitchAssetB;       
        TempswitchAssetB=size(ForntAssetIsB,1)+TempswitchAssetB;
        [SwitchToB,~,~]=find(EA(:,d/PeriodCell)< EB(:,d/PeriodCell));
        ReturnA(SwitchToB,(d+PeriodCell)/PeriodCell)    =   log(RBb(SwitchToB,d+PeriodCell+1)./RBb(SwitchToB,d+1));
        ReturnB(SwitchToB,(d+PeriodCell)/PeriodCell)    =   log(RBb(SwitchToB,d+PeriodCell+1)./RBb(SwitchToB,d+1));        
        [ForntAssetIsA,~,~]=find(ReturnA(SwitchToB,d/PeriodCell)==log(RAa(SwitchToB,d+1)./RAa(SwitchToB,d-PeriodCell+1)));       
        TotalswitchAssetA=size(ForntAssetIsA,1)+TotalswitchAssetA;
        TempswitchAssetA=size(ForntAssetIsA,1)+TempswitchAssetA;
        [ForntAssetIsA,~,~]=find(ReturnB(SwitchToB,d/PeriodCell)==log(RAa(SwitchToB,d+1)./RAa(SwitchToB,d-PeriodCell+1)));       
        TotalswitchAssetB=size(ForntAssetIsB,1)+TotalswitchAssetB;
        TempswitchAssetB=size(ForntAssetIsB,1)+TempswitchAssetB;
        [NoSwitch,~,~]=find(EA(:,d/PeriodCell)==EB(:,d/PeriodCell));
        ReturnA(NoSwitch,(d+PeriodCell)/PeriodCell)    =   log(RAa(NoSwitch,d+PeriodCell+1)./RAa(NoSwitch,d+1));
        ReturnB(NoSwitch,(d+PeriodCell)/PeriodCell)    =   log(RBb(NoSwitch,d+PeriodCell+1)./RBb(NoSwitch,d+1)); 
        if  mod(d+PeriodCell,4)==0
               PeriodswitchAssetA((d+PeriodCell)/N,1)= PeriodswitchAssetA((d+PeriodCell)/N,1)+TempswitchAssetA;   
               PeriodswitchAssetB((d+PeriodCell)/N,1)= PeriodswitchAssetB((d+PeriodCell)/N,1)+TempswitchAssetB;  
               TempswitchAssetA=0;
               TempswitchAssetB=0;
        end
end
AverageswitchA=sum(PeriodswitchAssetA)/T;
AverageswitchB=sum(PeriodswitchAssetB)/T;
clear('EA','EB');
for b = FirstCell : PeriodCell : ActualTime
    PeriodReturnA(:,b/PeriodCell)=exp(sum(ReturnA(:,1:b/PeriodCell),2));
    PeriodReturnB(:,b/PeriodCell)=exp(sum(ReturnB(:,1:b/PeriodCell),2));
end
%--------------%
TrueCashFlowA=zeros(M,ActualTime);
TrueCashFlowB=zeros(M,ActualTime);
Y1A=zeros(M,1,'single');
XxA=zeros(M,1,'single');
X1A=zeros(M,1,'single');
R1A=zeros(4,M,'single');
a2A=zeros(4,M,'single');
ContinValueA=zeros(M,1,'single');
Y1B=zeros(M,1,'single');
XxB=zeros(M,1,'single');
X1B=zeros(M,1,'single');
R1B=zeros(4,M,'single');
a2B=zeros(4,M,'single');
ContinValueB=zeros(M,1,'single');
%-EndPeriod CashFlow-%
TrueCashFlowA(:,end)=max(L0*(1+g)^T,L0*(RAa(:,end)./S0A));
TrueCashFlowB(:,end)=max(L0*(1+g)^T,L0*(RBb(:,end)./S0B)); 
%--------------------------%
for b =  (ActualTime-PeriodCell) : -PeriodCell : FirstCell
            IdxA = find((L0*(1+g)^T)*exp(-rf*(ActualTime-b)/PeriodCell*Celldt) < max(L0*(1+g)^((b/PeriodCell)*Celldt),L0.*PeriodReturnA(:,b/PeriodCell)));%
            IdxA=uint32(IdxA);
            IdxB = find((L0*(1+g)^T)*exp(-rf*(ActualTime-b)/PeriodCell*Celldt) < max( L0*(1+g)^((b/PeriodCell)*Celldt),L0.*PeriodReturnB(:,b/PeriodCell)));%
            IdxB=uint32(IdxB);
            %-T=LastPeriod-%
            if b == (ActualTime-PeriodCell)
                    Y1A=TrueCashFlowA(IdxA,(b+PeriodCell)).*exp(-rf*Celldt); %
                    Y1Aa=gpuArray(Y1A);
                    X1A=PeriodReturnA(IdxA,b/PeriodCell).*S0A;          %
                    %------------%
                    R1A =[ones(size(X1A))  (exp(-X1A./2)).*ones(size(X1A))  (exp(-X1A./2)).*(1.-X1A)  (exp(-X1A./2)).*(1.-(2.*X1A)+(X1A.^2))]; 
                    R1Aa = gpuArray(R1A); 
                    a1A = (R1A')*R1A;
                    a2A = a1A\(R1A');
                    ContinValueA= R1Aa*a2A*Y1Aa;
                    Y1B=TrueCashFlowB(IdxB,(b+PeriodCell)).*exp(-rf*Celldt); %
                    Y1Bb=gpuArray(Y1B);
                    X1B=PeriodReturnB(IdxB,b/PeriodCell).*S0B;
                    %------------%
                    R1B =[ones(size(X1B))  (exp(-X1B./2)).*ones(size(X1B))  (exp(-X1B./2)).*(1.-X1B)  (exp(-X1B./2)).*(1.-(2.*X1B)+(X1B.^2))]; 
                    R1Bb = gpuArray(R1B); 
                    a1B = (R1B')*R1B;
                    a2B = a1B\(R1B');
                    ContinValueB= R1Bb*a2B*Y1Bb;
                    JdxA=find(ContinValueA>max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnA(IdxA,b/PeriodCell)));   %
                    JdxB=find(ContinValueB>max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnB(IdxB,b/PeriodCell)));   %
                    TrueCashFlowA(IdxA(JdxA),b) = max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnA(IdxA(JdxA),b/PeriodCell));   % 
                    TrueCashFlowA(IdxA(JdxA),(b+PeriodCell):ActualTime)=0;  
                    TrueCashFlowB(IdxB(JdxB),b) = max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnB(IdxB(JdxB),b/PeriodCell));   % 
                    TrueCashFlowB(IdxB(JdxB),(b+PeriodCell):ActualTime)=0;     
            else
                    Y1A=zeros(size(IdxA));
                    Y1B=zeros(size(IdxB));
                    for f = (b+PeriodCell) : PeriodCell : ActualTime
                            [ExTimeA,~,~]=find(TrueCashFlowA(IdxA,f)~=0);                                    %[rowIndex,colIndex,value]=find(x)
                            [ExTimeB,~,~]=find(TrueCashFlowB(IdxB,f)~=0);
                            Y1A(ExTimeA,1)=TrueCashFlowA(IdxA(ExTimeA),f).*exp(-rf*((f-b)/PeriodCell)*Celldt); %Extime  , celldt=13/52;
                            Y1B(ExTimeB,1)=TrueCashFlowB(IdxB(ExTimeB),f).*exp(-rf*((f-b)/PeriodCell)*Celldt);
                    end
                    Y1Aa=gpuArray(Y1A);
                    X1A=PeriodReturnA(IdxA,b/PeriodCell).*S0A;  
                    %------------%
                    R1A =[ones(size(X1A))  (exp(-X1A./2)).*ones(size(X1A))  (exp(-X1A./2)).*(1.-X1A)  (exp(-X1A./2)).*(1.-(2.*X1A)+(X1A.^2))]; 
                    R1Aa = gpuArray(R1A);
                    a1A = (R1A')*R1A;
                    a2A = a1A\(R1A');
                    ContinValueA= R1Aa*a2A*Y1Aa;
                    %----------------------%
                    Y1Bb=gpuArray(Y1B);
                    X1B=PeriodReturnB(IdxB,b/PeriodCell).*S0B;
                    %------------%
                    R1B =[ones(size(X1B))  (exp(-X1B./2)).*ones(size(X1B))  (exp(-X1B./2)).*(1.-X1B)  (exp(-X1B./2)).*(1.-(2.*X1B)+(X1B.^2))]; 
                    R1Bb=gpuArray(R1B);
                    a1B = (R1B')*R1B;
                    a2B = a1B\(R1B');
                    ContinValueB= R1Bb*a2B*Y1Bb; 
                    JdxA=find(ContinValueA>max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnA(IdxA,b/PeriodCell)));   %
                    JdxB=find(ContinValueB>max(L0*(1+g)^((b/PeriodCell)*Celldt),L0*PeriodReturnB(IdxB,b/PeriodCell)));   %
                    TrueCashFlowA(IdxA(JdxA),b) = max(L0*(1+g)^((b/PeriodCell)*Celldt),L0.*PeriodReturnA(IdxA(JdxA),b/PeriodCell));   % 
                    TrueCashFlowA(IdxA(JdxA),(b+PeriodCell):ActualTime)=0;  
                    TrueCashFlowB(IdxB(JdxB),b) = max(L0*(1+g)^((b/PeriodCell)*Celldt),L0.*PeriodReturnB(IdxB(JdxB),b/PeriodCell));   % 
                    TrueCashFlowB(IdxB(JdxB),(b+PeriodCell):ActualTime)=0;      
            end  
end 
TotalPriceA=0;
TotalPriceB=0; 
nnA=0; %check Extime Distribue
nnB=0;
for DateCount = FirstCell : PeriodCell : ActualTime
        [RowFinalExTimeA,~,~] = find( TrueCashFlowA(:,DateCount) ~= 0 ); %
        PeriodDebtValueA=max(L0*(1+g)^(DateCount/N)- L0.*PeriodReturnA(RowFinalExTimeA,DateCount/PeriodCell),0);
        TempPrice1A=sum(PeriodDebtValueA.*exp(-rf*((DateCount/PeriodCell)*Celldt))); %%Celldt=13/52
        TotalPriceA=TempPrice1A+TotalPriceA ; %
        %-------------------%
        [RowFinalExTimeB,~,~] = find( TrueCashFlowB(:,DateCount) ~= 0 ); %
        PeriodDebtValueB=max(L0*(1+g)^(DateCount/N)-L0.*PeriodReturnB(RowFinalExTimeB,DateCount/PeriodCell),0);
        TempPrice1B=sum(PeriodDebtValueB.*exp(-rf*((DateCount/PeriodCell)*Celldt))); %%Celldt=13/52
        TotalPriceB=TempPrice1B+TotalPriceB; %
end
AmerFinalPriceA=TotalPriceA/M; %
AmerFinalPriceB=TotalPriceB/M; %