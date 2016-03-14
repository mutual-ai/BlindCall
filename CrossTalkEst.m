function [ k,CT ] = CrossTalkEst(  x,b,qual )
%CROSSTALKEST Summary of this function goes here
%   This function estimates the cross talk. bo

QUALITY=1;
QS_TH=0.5;
if(~exist('qual','var')||isempty(qual))
    QUALITY=0;
end
x=reshape(x,4,size(x,1)*size(x,2)/4);
b=reshape(b,4,size(b,1)*size(b,2)/4);

if QUALITY
    qual=reshape(qual,1,size(x,1)*size(x,2)/4);
    QS_QT=quantile(qual,QS_TH);
    Mask=qual>QS_QT;
    x=x(:,Mask);
    b=b(:,Mask);
    
end
%b(b<0)=0;
%x(x<0)=0;

%robust estimation based on quantiles
%{
k=zeros(4);
QT_val=0.5;%QT_val=0.5;
for i=1:4
    Idx=x(i,:)>0;
    for j=1:4        
        k(j,i)=quantile(b(j,Idx)./(x(i,Idx)+eps),QT_val);%%consider to use the Yaroslavsky filter to fine tune    
    end
end

for i=1:4
    k(:,i)=k(:,i)./k(i,i);
end

CT=[k(2,1),k(1,2),k(4,3),k(3,4)]; 
%}

%%%%%%%%%The below are the least squared estimates
%
xtx=zeros(4);
for i=1:4
    xtx(i,i)=sum(x(i,:).^2);
    for j=1:4
        xtb(i,j)=sum(x(i,:).*b(j,:));
    end
end

for i=1:4
    k2(:,i)=xtx\xtb(:,i);
end
k=k2';
for i=1:4
    k(:,i)=k(:,i)./k(i,i);
end
CT=[k(2,1),k(1,2),k(4,3),k(3,4)]; 

%}

end

