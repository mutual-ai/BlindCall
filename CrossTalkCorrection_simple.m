function [ x,qual,x_0,QTs ] = CrossTalkCorrection_simple( K,b,QTs )
%CROSSTALKCORRECTION Summary of this function goes here
%   Detailed explanation goes here

if ~exist('kappa','var')
    kappa = 100.0;%500
end

if ~exist('lambda','var')
    lambda = 2e-2;    
end


if(exist('QTs','var')&&~isempty(QTs))
    NORMALIZE_CHANNELS=0; 
else
    NORMALIZE_CHANNELS=1;
end


max_itr=2;


size_x=size(K,1);
%%%normalize 1
K=K+0.0;
for i=1:4
    K(:,i)=K(:,i)./K(i,i);
end


%%%normalize 2
%this doesn't seem to play well
%{
for i=1:4
    K(:,i)=K(:,i)./sum(K(:,i));
end
%}


%b(b<0)=0;


betamax = 1e8;
max_it=2;

size_y=size(b,2);

 
%%%%%%%%% ISD
%%find the ratios across channels and divide the cycle
%dependent normalization constant. Then estimate the proportion and error rate. 

size_b=size(b);
x=reshape(b,4,size(b,1)/4*size(b,2));
b=x;
beta =lambda;

itr=0;
LastIt=0;



%For the first round, no reg is better because reg suppress the
%intensity in an unwanted way.
LHS=(K'*K);
RHS=(K'*b);
x=LHS\RHS;





x_0=x;

x_tmp=x_0;

while beta < betamax && itr<max_it
    itr=itr+1;    
    x0=x;
    %x_tmp=x;
    %%%%%%%%%%%%%%%%%%%%%%%%%begin of penalty subproblem
    if(NORMALIZE_CHANNELS)
        
        if(itr<=1)
        
            A_Samples=x(1:4:size_x,:);
            C_Samples=x(2:4:size_x,:);
            G_Samples=x(3:4:size_x,:);
            T_Samples=x(4:4:size_x,:);
            A_Samples=A_Samples(:);
            C_Samples=C_Samples(:);
            G_Samples=G_Samples(:);
            T_Samples=T_Samples(:);

            QT_val=0.95;

            A_QT=quantile(A_Samples,QT_val);
            C_QT=quantile(C_Samples,QT_val);
            G_QT=quantile(G_Samples,QT_val);
            T_QT=quantile(T_Samples,QT_val);

            
            QTs=[A_QT,C_QT,G_QT,T_QT];
            %QTs./max(QTs)%%print the normalization factor

     
        else
            if(exist('QTs','var')&&~isempty(QTs))
                A_QT=QTs(1);
                C_QT=QTs(2);
                G_QT=QTs(3);
                T_QT=QTs(4);
            end
        end
            
    
        %%%This is for non maximum suppression, with the suppression we can get
        %%%a better support set, and we derive statistics in this set.


        for i=1:4:size_x-3
            tmp=x(i:i+3,:);

            %%normalize 
            tmp(1,:)=tmp(1,:)./A_QT;
            tmp(2,:)=tmp(2,:)./C_QT;
            tmp(3,:)=tmp(3,:)./G_QT;
            tmp(4,:)=tmp(4,:)./T_QT;

            %%pick the maximum
            [MaxI,~]=max(tmp,[],1);
            mask=bsxfun(@eq,tmp,MaxI);
            mask=~mask;
            tmp=x(i:i+3,:);
            tmp(mask)=0;
            x(i:i+3,:)=tmp;    
        end
        
        %%% below round is for more accurate support set estimation
        %we shall use the maximums in the above round
        
        
        x_r=reshape(x,size_b(1),size_b(2));
        x_tmp_r=reshape(x_tmp,size_b(1),size_b(2));
        %qual_tmp=reshape(qual,size_b(1)/4,size_b(2));
        R_max=5;
       [ QTs ] = NormV(x_r',x_tmp_r',[],R_max);
       %
       R_max=5;
        [ QTs1 ] = NormV( bsxfun(@rdivide,x_r',repmat(QTs,1,size(x_r,1)/4)),bsxfun(@rdivide,x_tmp_r',repmat(QTs,1,size(x_r,1)/4)),[] ,R_max);
        QTs=QTs.*QTs1;    
        %}
       %QTs./max(QTs)%%print again

        A_QT=QTs(1);
        C_QT=QTs(2);
        G_QT=QTs(3);
        T_QT=QTs(4);


    else
        
       x_r=reshape(x,size_b(1),size_b(2));
       x_tmp_r=reshape(x_tmp,size_b(1),size_b(2));
       if(itr>1)
           qual_tmp=reshape(qual,size_b(1)/4,size_b(2));
           [ QTs ] = NormV(x_r',x_tmp_r',qual_tmp');
       else
           [ QTs ] = NormV(x_r',x_tmp_r');
       end
           
        A_QT=QTs(1);
        C_QT=QTs(2);
        G_QT=QTs(3);
        T_QT=QTs(4);
    end

    %redo it
    x=x0;
    

    if(itr==1)
        qual=zeros(size_x/4,size(x,2));
        
        cycle=0;

        for i=1:4:size_x-3
            %
            cycle=cycle+1;
            tmp=x_tmp(i:i+3,:);

            tmp(1,:)=tmp(1,:)./A_QT;
            tmp(2,:)=tmp(2,:)./C_QT;
            tmp(3,:)=tmp(3,:)./G_QT;
            tmp(4,:)=tmp(4,:)./T_QT;

            [MaxI,~]=max(tmp,[],1);
            %%TODO: we need the non maximal channels to measure the confidence.

            mask=bsxfun(@eq,tmp,MaxI);
            sum_m=0;
            for r=1:4
                sum_m=sum_m+mask(r,:);
                mask(r,:)=mask(r,:).*(sum_m<=1);
            end


            tmp(mask)=0;%maximum suppression
            MaxI1=MaxI;
            MaxI1(MaxI<0)=0;
            [MaxI2,~]=max(tmp,[],1);%second largest
            MaxI2(MaxI2<0)=0;
            qual(cycle,:)=MaxI1./(MaxI1+MaxI2+eps);

        end

    end
    
   
    cycle=0;
    for i=1:4:size_x-3
        %
        cycle=cycle+1;
        tmp=x(i:i+3,:);
        tmp(1,:)=tmp(1,:)./A_QT;
        tmp(2,:)=tmp(2,:)./C_QT;
        tmp(3,:)=tmp(3,:)./G_QT;
        tmp(4,:)=tmp(4,:)./T_QT;
        [MaxI,~]=max(tmp,[],1);

        mask=bsxfun(@eq,tmp,MaxI);
        sum_m=0;
        for r=1:4
            sum_m=sum_m+mask(r,:);
            mask(r,:)=mask(r,:).*(sum_m<=1);
        end


        mask=~mask;
        
        
        %{
        qt_tmp=median(tmp,1);%quantile(tmp,0.6,1);
        r_tmp=tmp;
        for j=1:4
            r_tmp(j,:)=x(i+j-1,:).*(tmp(j,:)>qt_tmp);%.*((tmp(j,:)-qt_tmp)./tmp(j,:)) %L1
        end
        
        %}
            
        tmp=x(i:i+3,:);
        tmp(mask)=0;

        x(i:i+3,:)=tmp;  
        
        if exist('r_tmp','var')
            x(i:i+3,:)=r_tmp;
        end

    end

    
    
    
    
  %%%%%%%%%end of penalty sub problem


        
       
    LHS=(K'*K+beta*eye(length(K)));
    RHS=(K'*b+beta*x);


    x=LHS\RHS;
    

   
  
    
    if LastIt
        
        cycle=0;
        for i=1:4:size_x-3
            %
            cycle=cycle+1;
            tmp=x(i:i+3,:);
            tmp(1,:)=tmp(1,:)./A_QT;
            tmp(2,:)=tmp(2,:)./C_QT;
            tmp(3,:)=tmp(3,:)./G_QT;
            tmp(4,:)=tmp(4,:)./T_QT;
            [MaxI,MaxIdx]=max(tmp,[],1);

            %%TODO: we need the non maximal channels to measure the confidence.


            mask=bsxfun(@eq,tmp,MaxI);

            mask=~mask;
            tmp=x(i:i+3,:);
            tmp(mask)=0;

            x(i:i+3,:)=tmp;    
        end
    end

    
    
    beta=beta*kappa;
    

    if beta*kappa>betamax||itr>=max_itr-1
        LastIt=1;
    end
    %beta
    

end

x=reshape(x,size_b);
qual=reshape(qual,size_b(1)/4,size_b(2));

x_0=reshape(x_0,size_b);
end

