%matlabpool local 2;            % use two cores �p��
clear all;
%-------------Paper�ܼƸ��------------%
S0  = 36;                 %�i������
rf  = 0.06;               %�L���I�Q�v
sigma=0.2;                %�겣���i�ʲv
N  =50;                   %Number of points in time grid to use (minimum is 3, default is 50)
M  =10000;                %Number of points in asset price grid to use (minimum is 3, default is 50)
dt = 1/N;
T  = 1;
K  = 40;                  %��L0*(1+g)^T��K���N
%-------------�p�e�]�w�ܼƤ��W��---------------%
L0 = 100 ;                %�������B
g  = 0.06;                %�̧C���q�O��             

%----�B�J�@�G�����зǱ`�A�H���ü�(W)seed���ϥ�  (�u�n�����seed ���]����ɶ��I�]�X�����G���n�@��)-----%
s=1;
%-----���ͼзǱ`�A���t�üƯx�}-----%
rng(s);
WSimulate1=normrnd(0,1,M/2,N,'double');
WSimulate2=-1*WSimulate1;
TotalWSimulate=[WSimulate1;WSimulate2];
clear ('WSimulate1','WSimulate2');
ActualWSimulate=TotalWSimulate(1:M,1:N);
clear('TotalWSimulate');
%-----�B�J�G�G�����겣��ڻ�------%
Szero=ones(M,1,'double');
S=zeros (M,N,'double');
Sn=zeros(M,1,'double');
Ss=zeros(M,N+1,'double');
Sn(:,1)=S0.*Szero(:,1);
Ss=[Sn,S];                      %�X�֯x�} �Φ��@�� M x (N*T+1)���x�}
%-------�U���겣������--------%   
for c = 1 : 1 : N
    Ss(:,c+1)=Ss(:,c).*(exp((rf-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*ActualWSimulate(:,c)));
end
clear ('ActualWSimulate','Szeros','S','Sn');
%-------��X�̫�@�����H��ڮ��o�h��������|-------%
TrueCashFlow=zeros(M,N,'double');
Y1=zeros(M,1,'double');
Xx=zeros(M,1,'double');
X1=zeros(M,1,'double');
R1=zeros(4,M,'double');
a2=zeros(4,M,'double');
ContinValue=zeros(M,1,'double');
TrueCashFlow(:,end)=max(K,Ss(:,end)); 
for b =  (N-1) : -1 : 1

    Idx = find( K > (Ss(:,b+1)*exp(rf*(N-b)*dt)));%��XT-1�b���������|,�ðO���_��
    
    if b == (N-1)
        Y1=TrueCashFlow(Idx,(b+1)).*exp(-rf*dt); %Y1a����T����{��T-1�����{���y�q
        Xx=Ss(Idx,b+1);          %Xx����T-1�����ѻ� ���]�t��0�� �G+1
        X1=Xx/S0;                %X1a����T-1�����ֿn�h���
        %------�]�^�k------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1');
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;
        clear ('R1','a1','a2','Y1');
        %------------------%
    else                         % �q48���}�l ������H�W����{���|,�ݭn�t�~���j��ӧP�_�U�������g����{����
        Y1=zeros(size(Idx),1);
        for d = (b+1) : 1 : N
            [ExTime,~,~]=find(TrueCashFlow(Idx,d)~=0);                                    %Extime ��d���������|�D��0�� 
            Y1(ExTime,1)=TrueCashFlow(Idx(ExTime),d).*exp(-rf*(d-b)*dt); %�N��d��Extime�����|��{  ,���ۦb�j���U�@����Ӵ���Extime  , celldt=13/52;
        end
        Xx=Ss(Idx,b+1);   %Xx����T-1�����ѻ� ���]�t��0�� �G+1
        X1=Xx/S0;         %X1a����T-1�����ֿn�h���
        %----------�]�^�k--------%
        R1 =[ones(size(X1))  (exp(-X1./2)).*ones(size(X1))  (exp(-X1./2)).*(1.-X1)  (exp(-X1./2)).*(1.-(2.*X1)+(X1.^2))]; 
        clear ('Xx','X1'); 
        a1 = (R1')*R1;
        a2 = a1\(R1');
        ContinValue= R1*a2*Y1;       
        clear ('R1','a1','a2','Y1');
        %----------------------%
    end  
    Jdx = find(ContinValue.*exp(rf*(N-b)*dt) > K );  %�O���ߧY�i�����Ȥj���~��������Ȫ����| 
    clear('ContinValue');
    TrueCashFlow(Idx(Jdx),b) = Ss(Idx(Jdx),b+1);                   % �D�X�bT-1�ߧY�i��������
    TrueCashFlow(Idx(Jdx),(b+1):N)=0;                                 %����v�u��i���@�� ,�G�NT-1��������{���y�q���k�s
    clear ('Idx','Jdx');
end 
    TotalPrice=0;
for DateCount = 1 : 1 : N
    [RowFinalExTime,~,~] = find( TrueCashFlow(:,DateCount) ~= 0 ); %��M��e���{���y�q�����s�����|
    PeriodDebtValue=max(K-TrueCashFlow(RowFinalExTime,DateCount),0);
    TempPrice1=sum(PeriodDebtValue.*exp(-rf*(DateCount*dt))); %��{ %Celldt=13/52
    TotalPrice=TempPrice1+TotalPrice; %�N�Ȧs�JTotalprice 
end
FinalPrice=TotalPrice/M; %������