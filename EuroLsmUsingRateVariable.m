%matlabpool local 2; % use two cores 計算

L0 = 100 ;  %提撥金額
g  = 0  ;  %最低收率保證
T  = 30;     % 最久幾年退休
S0A  = 10    ;  %A標的資產t=0 初始價格
S0B  = 10    ;  %B標的資產t=0 初始價格

%用paper的資料模擬
rf  = 0.015;   %0.015假設的無風險利率
%無風險利率0.130131579  Period:2014~2015/03/20 52weeks T-bill coupon rate 來源: http://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=billrates 
BetaA =0.75;
BetaB =1.28;
RiskPremium = 0.0475; 
uA = rf+(BetaA*RiskPremium);  %A資產理論預期報酬率 
uB = rf+(BetaB*RiskPremium) ;  %B資產理論預期報酬率
sigmaA = 0.13;  %A資產的波動率 (Use Paper的數據)
sigmaB = 0.21;  %B資產的波動率 
N  =52;       %分割期數(一年分成52周) 
M  =10000;      %模擬次數
dt = 1/N; %Delta
%現實情況
% Rm-Rf Period: 2014/02~2015/02  來源:http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
%RiskPremium = 0.1353 ;

number=1; %模擬迴圈次數
TotalFinaldebtA =zeros(number,1);
TotalFinaldebtB =zeros(number,1);
BSputvalueA = zeros(number,1);
BSputvalueB = zeros(number,1);
standardErrorA = zeros(number,1);
standardErrorB = zeros(number,1);
%%%%%%%%-------設定變數名稱------%%%%%%%%%%%%

MaxT=40; %最久工作40年退休
TotalTime = N*MaxT; %最長總期間(周)
ActualTime= N*T; %實際工作幾年為退休
SimulateTimes = M ; %模擬次數
%----時間間隔-----%
TimeCell = ActualTime/13; %總期間(周)
PeriodCell = 13;  %期間周期(周)
FirstCell = 13; %第一季開始點(周)
Celldt = 13/52; %每季換資產的dt (時間間隔)

corrAB=0;
for t=1:1:1
    %步驟一：模擬標準常態隨機亂數(W) 有錯誤  seed的使用  (只要控制住seed 不館任何時間點跑出的結果都要一樣)
    %固定seed
    s=1;
    v=2;
   
    %產生標準常態分配亂數矩陣
    rng(s);
    WAOrigin=normrnd(0,1,SimulateTimes,TotalTime);
    % alpha11= 1 ;
    % WAs=alpha11.*WAOrigin;
    WASimulate=WAOrigin(1:SimulateTimes,1:ActualTime);
    % %---B資產---%
    %步驟一：模擬標準常態隨機亂數(W)
    %固定seed
    rng(v);
    WBOrigin=normrnd(0,1,SimulateTimes,TotalTime);  %隨機項獨立與不獨立  , s1=X1 s2=px1-X2*(1-p^2)^1/2 (兩資產)  / 三資產的不同
    WBs=zeros(SimulateTimes,TotalTime);
    alpha21=corrAB;     	%a21*a11=1
    alpha22=(1-alpha21^2)^(1/2);  %a21^2+a22^2=1
    WBs=(alpha21.*WAOrigin)+(alpha22.*WBOrigin);
    WBSimulate=WBs(1:SimulateTimes,1:ActualTime);
    %步驟二：模擬資產實際價格
    YhatA = zeros(SimulateTimes,ActualTime);
    YhatB = zeros(SimulateTimes,ActualTime);
    TotalPensionA = zeros(SimulateTimes,1) ;
    TotalPensionB = zeros(SimulateTimes,1) ;
    EA = zeros(SimulateTimes,TimeCell); %預期報酬   每季換一次的話就要以13周為間隔
    EB = zeros(SimulateTimes,TimeCell);
    L1A=zeros(SimulateTimes,1);
    L1B=zeros(SimulateTimes,1);
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
    SAa=[SnA,SA]; %合併矩陣 形成一個 M x (N*T+1)的矩陣
    SBb=[SnB,SB]; %合併矩陣 形成一個 M x (N*T+1)的矩陣
    RAa=[RnA,SA]; %合併矩陣 形成一個 M x (N*T+1)的矩陣
    RBb=[RnB,SB]; %合併矩陣 形成一個 M x (N*T+1)的矩陣    
    Y1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    X1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);     
    Y1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    X1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    ReturnA=zeros(SimulateTimes,TimeCell);
    ReturnB=zeros(SimulateTimes,TimeCell);
    %第一期累積退休金
    %情境A 假設第0期投資A資產
    for c = 1:1:ActualTime
        SAa(:,c+1)=SAa(:,c).*(exp((rf-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*WASimulate(:,c)));
        SBb(:,c+1)=SBb(:,c).*(exp((rf-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*WBSimulate(:,c))); 
        RAa(:,c+1)=RAa(:,c).*(exp((rf-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*WASimulate(:,c)));
        RBb(:,c+1)=RBb(:,c).*(exp((rf-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*WBSimulate(:,c)));
    end
    %步驟三：計算各期的預測價格
    for b = FirstCell : PeriodCell : ActualTime-PeriodCell
        Y1a(:,b/PeriodCell)=SAa(:,b+PeriodCell+1)./SAa(:,b+1); %由於有考慮第0周 , 故這裡使用b+1 ,起點是從第1期開始
        X1a(:,b/PeriodCell)=SAa(:,b+1)./SAa(:,b-PeriodCell+1);        
        R1a =[ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-X1a(:,b/PeriodCell))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-(2.*X1a(:,b/PeriodCell)))+(X1a(:,b/PeriodCell).^2)]; 
        a1 = (R1a')*R1a;
        a2 = a1\(R1a');
        YhatA(:,b) = R1a*a2*Y1a(:,b/PeriodCell);
        Y1b(:,b/PeriodCell)=SBb(:,b+PeriodCell+1)./SBb(:,b+1); %由於有考慮第0周 , 故這裡使用b+1 ,起點是從第1期開始
        X1b(:,b/PeriodCell)=SBb(:,b+1)./SBb(:,b-PeriodCell+1);     
        R1b =[ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-X1b(:,b/PeriodCell))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-(2.*X1b(:,b/PeriodCell)))+(X1b(:,b/PeriodCell).^2)];
        b1 = (R1b')*R1b;
        b2 = b1\(R1b');
        YhatB(:,b) = R1b*b2*Y1b(:,b/PeriodCell);
    end

       %步驟五    
     for c = 1:1:SimulateTimes
        for d = FirstCell:PeriodCell:ActualTime-PeriodCell
            if YhatA(c,d)>YhatB(c,d)
                ReturnA(c,d/PeriodCell)=log(RAa(c,d+PeriodCell+1)/RAa(c,d+1));
                ReturnB(c,d/PeriodCell)=log(RAa(c,d+PeriodCell+1)/RAa(c,d+1));
            elseif YhatB(c,d)>YhatA(c,d)
                ReturnA(c,d/PeriodCell)=log(RBb(c,d+PeriodCell+1)/RBb(c,d+1));
                ReturnB(c,d/PeriodCell)=log(RBb(c,d+PeriodCell+1)/RBb(c,d+1));
            elseif YhatA(c,d)==YhatB(c,d)
                if d>FirstCell
                    if     YhatA(c,d-PeriodCell)==YhatA(c,d)
                            ReturnA(c,d/PeriodCell)=log(RAa(c,d+PeriodCell+1)/RAa(c,d+1));
                            ReturnB(c,d/PeriodCell)=log(RAa(c,d+PeriodCell+1)/RAa(c,d+1));
                    elseif YhatB(c,d-PeriodCell)==YhatB(c,d)
                            ReturnA(c,d/PeriodCell)=log(RBb(c,d+PeriodCell+1)/RBb(c,d+1));
                            ReturnB(c,d/PeriodCell)=log(RBb(c,d+PeriodCell+1)/RBb(c,d+1));
                    end
                else
                     ReturnA(c,d/PeriodCell)=log(RAa(c,d+PeriodCell+1)/RAa(c,d+1));
                     ReturnB(c,d/PeriodCell)=log(RBb(c,d+PeriodCell+1)/RBb(c,d+1));                    
                end
            end
        end
     end

    % %步驟六：計算風險中立下的累積退休金
    parfor c = 1:SimulateTimes
        L1A(c,1) = L0*(RAa(c,FirstCell+1)./S0A);
        L1B(c,1) = L0*(RBb(c,FirstCell+1)./S0B); 
    end
    %各期累積退休金(考慮轉換資產)
    parfor c = 1:SimulateTimes
         TotalPensionA(c,1)=(exp(sum(ReturnA(c,:))))*L1A(c,1) ; 
         TotalPensionB(c,1)=(exp(sum(ReturnB(c,:))))*L1B(c,1) ; 
    end
    %潛在負債現值平均(考慮轉換資產)
    debtA=(exp(-rf*T))*max((L0*((1+g)^T)-TotalPensionA),0);

    TotalFinaldebtA(t,1) = sum(debtA)/SimulateTimes;
   
    %潛在負債現值平均(考慮轉換資產)
    debtB=(exp(-rf*T))*max(((L0*(1+g)^T)-TotalPensionB),0);
    TotalFinaldebtB(t,1) = sum(debtB)/SimulateTimes;

    K=L0*((1+g)^T);
    bsS0=L0;
    %一開始投資A資產的BS Put value
    d1A =(log(bsS0/K)+(rf*T))/(sigmaA*sqrt(T))+1/2*sigmaA*sqrt(T);
    d2A = d1A - (sigmaA*sqrt(T));
    BSputvalueA(t,1) = K*exp(-rf*T)*normcdf(-d2A,0,1)-bsS0*normcdf(-d1A,0,1);  %賣權 為N(-d1) N(-d2) 買權 不用負


    %一開始投資B資產的BS Put value
    d1B =(log(bsS0/K)+(rf*T))/(sigmaB*sqrt(T))+1/2*sigmaB*sqrt(T);
    d2B = d1B - (sigmaB*sqrt(T));
    BSputvalueB(t,1) = K*exp(-rf*T)*normcdf(-d2B,0,1)-bsS0*normcdf(-d1B,0,1);

    %標準誤
    standardErrorA(t,1)=std(debtA)/sqrt(SimulateTimes);
    standardErrorB(t,1)=std(debtB)/sqrt(SimulateTimes);

%-------標準差變化變數--------%
 %     sigmaA = sigmaA+0.005;
 %     sigmaB = sigmaB+0.005;

%---------時間間隔變化變數--------%
%       TimeCell = TimeCell/2; %總期間(周)
%       PeriodCell = PeriodCell*2;  %期間周期(周)
%       FirstCell = FirstCell*2; %第一季開始點(周)
%       Celldt = Celldt*2; %每季換資產的dt (時間間隔)
%--------退休時間變化變數---------%
%      T=T+5;
%      ActualTime=N*T;
%      TimeCell = ActualTime/13; %總期間(周)
%-------初始資產價格變化變數--------%
%         S0A=S0A+5;
%         S0B=S0B+5;
%         corrAB=corrAB+0.1;
%-------Rf變化變數--------%         
%         rf=rf+0.005;
         t

end