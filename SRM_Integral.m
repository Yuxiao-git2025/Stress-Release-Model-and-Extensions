function CumInt=SRM_Integral(time,mag,Mmin,param)
% SRM cumulative integrated intensity at each event time
% CumInt(k)=∫_{T1}^{time(k)}lambda(u)du,lambda(t)=exp(A+B*(t-C*St(t)))
time=time(:);
mag=mag(:);
[time,ord]=sort(time,'ascend');
mag=mag(ord);
A=param(1);
B=param(2);
C=param(3);
eta=0.75;
Si=10.^(eta.*(mag-Mmin));
T1=time(1);
N=numel(time);
Ind=time<time.';        % N x N,Ind(i,k)=1 if time(i)<time(k)
St=Si.'*Ind;            % 1 x N,St(k)=sum_{i:ti<tk}Si
prev=[T1;time(1:end-1)];
dt=time-prev;
dInt=zeros(N,1);
if abs(B)<1e-12
    dInt(:)=exp(A).*dt;
else
    % On interval (prev(k),time(k)] stress is constant =St(k)
    K=exp(A-B*C.*St);   % 1 x N
    dInt(:)=(K.*exp(B.*prev.')).*(expm1(B.*dt.')./B);
end
CumInt=cumsum(dInt);
end
