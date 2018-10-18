%matlabpool local 2;            % use two cores 計算
clear all;
%-------------Paper變數資料------------%
S0  = 36;                 %履約價格
rf  = 0.06;               %無風險利率
sigma=0.2;                %資產的波動率
N  =50;                   %Number of points in time grid to use (minimum is 3, default is 50)
M  =10000;                %Number of points in asset price grid to use (minimum is 3, default is 50)
dt = 1/N;
T  = 1;
K  = 40;                  %使L0*(1+g)^T用K取代
%-------------計畫設定變數之名稱---------------%
L0 = 100 ;                %提撥金額
g  = 0.06;                %最低收益保證             

%----步驟一：模擬標準常態隨機亂數(W)seed的使用  (只要控制住seed 不館任何時間點跑出的結果都要一樣)-----%
s=1;
%-----產生標準常態分配亂數矩陣-----%
rng(s);
WSimulate1=normrnd(0,1,M/2,N,'double');
WSimulate2=-1*WSimulate1;
TotalWSimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
ActualWSimulate=TotalWSimulate(1:M,1:N);
clear('TotalWSimulate');
%-----步驟二：模擬資產實際價------%
Szero=ones(M,1,'double');
S=zeros (M,N,'double');
Sn=zeros(M,1,'double');
Ss=zeros(M,N+1,'double');
Sn(:,1)=S0.*Szero(:,1);
Ss=[Sn,S];                      %合併矩陣 形成一個 M x (N*T+1)的矩陣
%-------各期資產的價格--------%   
for c = 1 : 1 : N
    Ss(:,c+1)=Ss(:,c).*(exp((rf-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*ActualWSimulate(:,c)));
end
clear ('ActualWSimulate','Szeros','S','Sn');
%-------找出最後一期投資人實際拿得退休金的路徑-------%
TrueCashFlow=zeros(M,N,'double');
Y1=zeros(M,1,'double');
Xx=zeros(M,1,'double');
X1=zeros(M,1,'double');
R1=zeros(4,M,'double');
a2=zeros(4,M,'double');
ContinValue=zeros(M,1,'double');
TrueCashFlow(:,end)=max(K,Ss(:,end)); 
for b =  (N-1) : -1 : 1

    Idx = find( K > (Ss(:,b+1)*exp(rf*(N-b)*dt)));%找出T-1在價內的路徑,並記錄起來
    
    if b == (N-1)
        Y1=TrueCashFlow(Idx,(b+1)).*exp(-rf*dt); %Y1a為第T期折現到T-1期的現金流量
        Xx=Ss(Idx,b+1);          %Xx為第T-1期的股價 有包含第0期 故+1
        X1=Xx/S0;                %X1a為第T-1期的累積退休金
        %------跑回歸------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1');
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;
        clear ('R1','a1','a2','Y1');
        %------------------%
    else                         % 從48期開始 有兩期以上的折現路徑,需要另外的迴圈來判斷各價內路經的折現期數
        Y1=zeros(size(Idx),1);
        for d = (b+1) : 1 : N
            [ExTime,~,~]=find(TrueCashFlow(Idx,d)~=0);                                    %Extime 第d期價內路徑非為0者 
            Y1(ExTime,1)=TrueCashFlow(Idx(ExTime),d).*exp(-rf*(d-b)*dt); %將第d期Extime的路徑折現  ,接著在迴圈到下一期找該期的Extime  , celldt=13/52;
        end
        Xx=Ss(Idx,b+1);   %Xx為第T-1期的股價 有包含第0期 故+1
        X1=Xx/S0;         %X1a為第T-1期的累積退休金
        %----------跑回歸--------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1'); 
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;       
        clear ('R1','a1','a2','Y1');
        %----------------------%
    end  
    Jdx = find(ContinValue.*exp(rf*(N-b)*dt) > K );  %記錄立即履約價值大於繼續持有價值的路徑 
    clear('ContinValue');
    TrueCashFlow(Idx(Jdx),b) = Ss(Idx(Jdx),b+1);                   % 求出在T-1立即履約的價值
    TrueCashFlow(Idx(Jdx),(b+1):N)=0;                                 %選擇權只能履約一次 ,故將T-1之後期的現金流量均歸零
    clear ('Idx','Jdx');
end 
    TotalPrice=0;
for DateCount = 1 : 1 : N
    [RowFinalExTime,~,~] = find( TrueCashFlow(:,DateCount) ~= 0 ); %找尋第e期現金流量不為零的路徑
    PeriodDebtValue=max(K-TrueCashFlow(RowFinalExTime,DateCount),0);
    TempPrice1=sum(PeriodDebtValue.*exp(-rf*(DateCount*dt))); %折現 %Celldt=13/52
    TotalPrice=TempPrice1+TotalPrice; %將值存入Totalprice 
end
FinalPrice=TotalPrice/M; %取平均