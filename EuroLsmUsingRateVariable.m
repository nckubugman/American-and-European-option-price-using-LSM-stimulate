%matlabpool local 2; % use two cores �p��

L0 = 100 ;  %�������B
g  = 0  ;  %�̧C���v�O��
T  = 30;     % �̤[�X�~�h��
S0A  = 10    ;  %A�Ъ��겣t=0 ��l����
S0B  = 10    ;  %B�Ъ��겣t=0 ��l����

%��paper����Ƽ���
rf  = 0.015;   %0.015���]���L���I�Q�v
%�L���I�Q�v0.130131579  Period:2014~2015/03/20 52weeks T-bill coupon rate �ӷ�: http://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=billrates 
BetaA =0.75;
BetaB =1.28;
RiskPremium = 0.0475; 
uA = rf+(BetaA*RiskPremium);  %A�겣�z�׹w�����S�v 
uB = rf+(BetaB*RiskPremium) ;  %B�겣�z�׹w�����S�v
sigmaA = 0.13;  %A�겣���i�ʲv (Use Paper���ƾ�)
sigmaB = 0.21;  %B�겣���i�ʲv 
N  =52;       %���δ���(�@�~����52�P) 
M  =10000;      %��������
dt = 1/N; %Delta
%�{�걡�p
% Rm-Rf Period: 2014/02~2015/02  �ӷ�:http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
%RiskPremium = 0.1353 ;

number=1; %�����j�馸��
TotalFinaldebtA =zeros(number,1);
TotalFinaldebtB =zeros(number,1);
BSputvalueA = zeros(number,1);
BSputvalueB = zeros(number,1);
standardErrorA = zeros(number,1);
standardErrorB = zeros(number,1);
%%%%%%%%-------�]�w�ܼƦW��------%%%%%%%%%%%%

MaxT=40; %�̤[�u�@40�~�h��
TotalTime = N*MaxT; %�̪��`����(�P)
ActualTime= N*T; %��ڤu�@�X�~���h��
SimulateTimes = M ; %��������
%----�ɶ����j-----%
TimeCell = ActualTime/13; %�`����(�P)
PeriodCell = 13;  %�����P��(�P)
FirstCell = 13; %�Ĥ@�u�}�l�I(�P)
Celldt = 13/52; %�C�u���겣��dt (�ɶ����j)

corrAB=0;
for t=1:1:1
    %�B�J�@�G�����зǱ`�A�H���ü�(W) �����~  seed���ϥ�  (�u�n�����seed ���]����ɶ��I�]�X�����G���n�@��)
    %�T�wseed
    s=1;
    v=2;
   
    %���ͼзǱ`�A���t�üƯx�}
    rng(s);
    WAOrigin=normrnd(0,1,SimulateTimes,TotalTime);
    % alpha11= 1 ;
    % WAs=alpha11.*WAOrigin;
    WASimulate=WAOrigin(1:SimulateTimes,1:ActualTime);
    % %---B�겣---%
    %�B�J�@�G�����зǱ`�A�H���ü�(W)
    %�T�wseed
    rng(v);
    WBOrigin=normrnd(0,1,SimulateTimes,TotalTime);  %�H�����W�߻P���W��  , s1=X1 s2=px1-X2*(1-p^2)^1/2 (��겣)  / �T�겣�����P
    WBs=zeros(SimulateTimes,TotalTime);
    alpha21=corrAB;     	%a21*a11=1
    alpha22=(1-alpha21^2)^(1/2);  %a21^2+a22^2=1
    WBs=(alpha21.*WAOrigin)+(alpha22.*WBOrigin);
    WBSimulate=WBs(1:SimulateTimes,1:ActualTime);
    %�B�J�G�G�����겣��ڻ���
    YhatA = zeros(SimulateTimes,ActualTime);
    YhatB = zeros(SimulateTimes,ActualTime);
    TotalPensionA = zeros(SimulateTimes,1) ;
    TotalPensionB = zeros(SimulateTimes,1) ;
    EA = zeros(SimulateTimes,TimeCell); %�w�����S   �C�u���@�����ܴN�n�H13�P�����j
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
    SAa=[SnA,SA]; %�X�֯x�} �Φ��@�� M x (N*T+1)���x�}
    SBb=[SnB,SB]; %�X�֯x�} �Φ��@�� M x (N*T+1)���x�}
    RAa=[RnA,SA]; %�X�֯x�} �Φ��@�� M x (N*T+1)���x�}
    RBb=[RnB,SB]; %�X�֯x�} �Φ��@�� M x (N*T+1)���x�}    
    Y1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    X1a=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);     
    Y1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    X1b=zeros(SimulateTimes, (ActualTime-PeriodCell)/PeriodCell);
    ReturnA=zeros(SimulateTimes,TimeCell);
    ReturnB=zeros(SimulateTimes,TimeCell);
    %�Ĥ@���ֿn�h���
    %����A ���]��0�����A�겣
    for c = 1:1:ActualTime
        SAa(:,c+1)=SAa(:,c).*(exp((rf-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*WASimulate(:,c)));
        SBb(:,c+1)=SBb(:,c).*(exp((rf-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*WBSimulate(:,c))); 
        RAa(:,c+1)=RAa(:,c).*(exp((rf-(sigmaA^2)*0.5).*dt+sigmaA*sqrt(dt).*WASimulate(:,c)));
        RBb(:,c+1)=RBb(:,c).*(exp((rf-(sigmaB^2)*0.5).*dt+sigmaB*sqrt(dt).*WBSimulate(:,c)));
    end
    %�B�J�T�G�p��U�����w������
    for b = FirstCell : PeriodCell : ActualTime-PeriodCell
        Y1a(:,b/PeriodCell)=SAa(:,b+PeriodCell+1)./SAa(:,b+1); %�ѩ󦳦Ҽ{��0�P , �G�o�̨ϥ�b+1 ,�_�I�O�q��1���}�l
        X1a(:,b/PeriodCell)=SAa(:,b+1)./SAa(:,b-PeriodCell+1);        
        R1a =[ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*ones(size(X1a(:,b/PeriodCell)))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-X1a(:,b/PeriodCell))  (exp(-X1a(:,b/PeriodCell)./2)).*(1.-(2.*X1a(:,b/PeriodCell)))+(X1a(:,b/PeriodCell).^2)]; 
        a1 = (R1a')*R1a;
        a2 = a1\(R1a');
        YhatA(:,b) = R1a*a2*Y1a(:,b/PeriodCell);
        Y1b(:,b/PeriodCell)=SBb(:,b+PeriodCell+1)./SBb(:,b+1); %�ѩ󦳦Ҽ{��0�P , �G�o�̨ϥ�b+1 ,�_�I�O�q��1���}�l
        X1b(:,b/PeriodCell)=SBb(:,b+1)./SBb(:,b-PeriodCell+1);     
        R1b =[ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*ones(size(X1b(:,b/PeriodCell)))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-X1b(:,b/PeriodCell))  (exp(-X1b(:,b/PeriodCell)./2)).*(1.-(2.*X1b(:,b/PeriodCell)))+(X1b(:,b/PeriodCell).^2)];
        b1 = (R1b')*R1b;
        b2 = b1\(R1b');
        YhatB(:,b) = R1b*b2*Y1b(:,b/PeriodCell);
    end

       %�B�J��    
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

    % %�B�J���G�p�⭷�I���ߤU���ֿn�h���
    parfor c = 1:SimulateTimes
        L1A(c,1) = L0*(RAa(c,FirstCell+1)./S0A);
        L1B(c,1) = L0*(RBb(c,FirstCell+1)./S0B); 
    end
    %�U���ֿn�h���(�Ҽ{�ഫ�겣)
    parfor c = 1:SimulateTimes
         TotalPensionA(c,1)=(exp(sum(ReturnA(c,:))))*L1A(c,1) ; 
         TotalPensionB(c,1)=(exp(sum(ReturnB(c,:))))*L1B(c,1) ; 
    end
    %��b�t�Ų{�ȥ���(�Ҽ{�ഫ�겣)
    debtA=(exp(-rf*T))*max((L0*((1+g)^T)-TotalPensionA),0);

    TotalFinaldebtA(t,1) = sum(debtA)/SimulateTimes;
   
    %��b�t�Ų{�ȥ���(�Ҽ{�ഫ�겣)
    debtB=(exp(-rf*T))*max(((L0*(1+g)^T)-TotalPensionB),0);
    TotalFinaldebtB(t,1) = sum(debtB)/SimulateTimes;

    K=L0*((1+g)^T);
    bsS0=L0;
    %�@�}�l���A�겣��BS Put value
    d1A =(log(bsS0/K)+(rf*T))/(sigmaA*sqrt(T))+1/2*sigmaA*sqrt(T);
    d2A = d1A - (sigmaA*sqrt(T));
    BSputvalueA(t,1) = K*exp(-rf*T)*normcdf(-d2A,0,1)-bsS0*normcdf(-d1A,0,1);  %���v ��N(-d1) N(-d2) �R�v ���έt


    %�@�}�l���B�겣��BS Put value
    d1B =(log(bsS0/K)+(rf*T))/(sigmaB*sqrt(T))+1/2*sigmaB*sqrt(T);
    d2B = d1B - (sigmaB*sqrt(T));
    BSputvalueB(t,1) = K*exp(-rf*T)*normcdf(-d2B,0,1)-bsS0*normcdf(-d1B,0,1);

    %�зǻ~
    standardErrorA(t,1)=std(debtA)/sqrt(SimulateTimes);
    standardErrorB(t,1)=std(debtB)/sqrt(SimulateTimes);

%-------�зǮt�ܤ��ܼ�--------%
 %     sigmaA = sigmaA+0.005;
 %     sigmaB = sigmaB+0.005;

%---------�ɶ����j�ܤ��ܼ�--------%
%       TimeCell = TimeCell/2; %�`����(�P)
%       PeriodCell = PeriodCell*2;  %�����P��(�P)
%       FirstCell = FirstCell*2; %�Ĥ@�u�}�l�I(�P)
%       Celldt = Celldt*2; %�C�u���겣��dt (�ɶ����j)
%--------�h��ɶ��ܤ��ܼ�---------%
%      T=T+5;
%      ActualTime=N*T;
%      TimeCell = ActualTime/13; %�`����(�P)
%-------��l�겣�����ܤ��ܼ�--------%
%         S0A=S0A+5;
%         S0B=S0B+5;
%         corrAB=corrAB+0.1;
%-------Rf�ܤ��ܼ�--------%         
%         rf=rf+0.005;
         t

end