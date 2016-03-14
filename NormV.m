function [ V ] = NormV( I0,I1,qual,R_max )
%NORMV Summary of this function goes here
%   Detailed explanation goes here
%I0 base call
%I1 raw signal


NullCallTh=200;
norm_cycles=min(4,size(I0,2)/4);
%intensities=cell(1,4);
if(exist('qual','var')&&(~isempty(qual)))
    if norm_cycles<10
        qual_qt=quantile(qual(:),0.2);
    else
    	qual_qt=quantile(qual(:),0.2);
    end   
end


if(~exist('R_max','var'))
   R_max=10;
end

MaxLen=norm_cycles;%size(I0,2)/4;
n_reads=size(I0,1);
channel=zeros(size(I0,1),MaxLen);
tmp=I0;
cycle=0;
Intensities=channel;
for j=1:4:MaxLen*4
    cycle=cycle+1;
    channels_tmp=tmp(:,j:j+3);
    [MaxIntensity,max_idx]=max(channels_tmp,[],2);
    channel(:,cycle)=max_idx;
    channel((MaxIntensity<=NullCallTh),cycle)=0;%%null calls
    Intensities(:,cycle)=MaxIntensity;
end
%{
mean_intensities=zeros(1,4);
for c=1:4
   mean_intensities(c)=mean(Intensities(channel==c));
end
mean_intensities./max(mean_intensities)
%}

%%C(4,2)=6

A=zeros(6*size(I0,1),4);
A2=zeros(6*size(I0,1),4);
Idx1=zeros(1,6);
Idx2=zeros(1,6);
cnt=0;

for i=1:4
    for j=1:4
        if(i<j)
            cnt=cnt+1;
            Idx1(cnt)=i;
            Idx2(cnt)=j;
        end
    end
end

row_idx=[0,1,2,3;0,0,4,5;0,0,0,6;0,0,0,0];
D=zeros(1,4);

for cycle=1:norm_cycles-1
        
    cur_idx=channel(:,cycle);
    next_idx=channel(:,cycle+1);
    next_idx(cur_idx==0)=0;
    cur_idx(next_idx==0)=0;
    
    
    if(exist('qual','var')&&(~isempty(qual)))
       qual_cur=qual(:,cycle);
       qual_next=qual(:,cycle+1);
    end
    
    larger_idx=(cur_idx>next_idx);
    tmp=cur_idx(larger_idx);
    cur_idx(larger_idx)=next_idx(larger_idx);
    next_idx(larger_idx)=tmp;
    
    rows_idx1=(cur_idx)+(next_idx-1).*4;
    rows_idx1(rows_idx1<=0)=1;
    rows1=row_idx(rows_idx1);
      
    rows1(rows1==0)=7;
    %
    nn=hist(rows1,7);
    %nn(7)=[];
    w=1./(nn+1);
    w=(w./max(w(1:end-1))).^.5;
    
    %%w is reshaped
    w=w(rows1);
    w0=w;
    %}
    %%%make the range within bound
    
    
    diff_idx=(cur_idx~=next_idx);
    diff_idx=diff_idx&(rows1~=7);
    rows1(rows1==7)=1;
    
    if(exist('qual','var')&&(~isempty(qual)))
       diff_idx=(diff_idx&((qual_cur>qual_qt)&(qual_next>qual_qt)));
    end
    
    
   
    diff_idx=(diff_idx&(Intensities(:,cycle)>NullCallTh)&(Intensities(:,cycle+1)>NullCallTh)); 
    diff_idx=(diff_idx&(Intensities(:,cycle)<R_max*Intensities(:,cycle+1))&(Intensities(:,cycle)>(Intensities(:,cycle+1)/R_max)));
    
    ratio=Intensities(diff_idx,cycle)./(Intensities(diff_idx,1+cycle)+eps);
    
    larger_idx=larger_idx(diff_idx);
    ratio(larger_idx)=1./ratio(larger_idx);
    
    %max(ratio(:))
    w=w(diff_idx);
    
    entry1=[0:n_reads-1]*6+rows1'+(Idx1(rows1)-1)*(6*n_reads);
    entry2=[0:n_reads-1]*6+rows1'+(Idx2(rows1)-1)*(6*n_reads);

    entry1=entry1(diff_idx);    
    entry2=entry2(diff_idx);

    %w=1;
    A(entry1)=(A(entry1)+1).*w;
    A(entry2)=(A(entry2)-ratio').*w;
    
    A2(entry2)=(A2(entry2)+1).*w;
    A2(entry1)=(A2(entry1)-1./ratio').*w;
    
    %
    for c=1:4
        D(c)=D(c)+2*nnz((Idx1(rows1)==c).*diff_idx');
        D(c)=D(c)+2*nnz((Idx2(rows1)==c).*diff_idx');
    end
    %}
    %{
    for c=1:4
        D(c)=D(c)+2*sum(w0.*(Idx1(rows1)==c));
        D(c)=D(c)+2*sum(w0.*(Idx2(rows1)==c));
    end
    %}
    
    
    
end


%{
%the non-normalized solution
AtA=A'*A;
[V,d]=eig(AtA);
V=V(:,1)./sum(V(:,1));
V=V(:)';
%}

%%the normalized solution
%
AtA=A'*A+A2'*A2;
D=diag(D);

[V,d]=eig(AtA,D);
V=abs(V(:,1))./max(abs(V(:,1)));
V=V(:)';
%}



end

