function [ x,qual,x_tmp,QTs ] = Deconv_ISD_simple( K,b,lambda, Final,QTs, FAST )
%DECONV_ISD Summary of this function goes here
%   get deconvolved signal x from blurry observation B.

if ~exist('Final','var')
    Final = 0;
end

if ~exist('FAST','var')
    FAST = 0;%100
end

if ~exist('kappa','var')
    if Final
        kappa = 20.0;%100 %smaller kappa uses more time but gives better results
    else
        kappa = 10.0;%100 %smaller kappa uses more time but gives better results
    end
   
end

if ~exist('lambda','var')
    lambda = 2e-2;
    
end

if Final
    betamax = 1e1;
    max_itr=0;
else
     betamax = 1e6;
     max_itr=7;
end

if(exist('QTs','var')&&~isempty(QTs))
    NORMALIZE_CHANNELS=0; 
else
    NORMALIZE_CHANNELS=1;
end

size_x=size(K,1);
   
%to do ISD: 
%coarse round: find a quantile value and normalize all channels, pick the strongest
%channel using the normalized intensity.
%fine round: use the above result to accurately determine the normalization
%factor

x=b;
beta =2*lambda;%2*
itr=0;
LastIt=0;

while beta < betamax &&itr<=max_itr
    %QTs
    itr=itr+1;
    x0=x; 
    %%%%%%%%%begin of penalty sub problem
    
    if(NORMALIZE_CHANNELS)        
        if (itr<=2)
        
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

            if(itr>1)
                QTs=[A_QT,C_QT,G_QT,T_QT];
                %QTs./max(QTs)%%print the normalization factor
            end
     
        else
            
            A_QT=QTs(1);
            C_QT=QTs(2);
            G_QT=QTs(3);
            T_QT=QTs(4);
            
        end
            
    
        %%%This is for non maximum suppression, with the suppression we can get
        %%%a better support set, and we derive statistics in this set.

        if(itr>1)
            for i=1:4:size_x-3
                tmp=x(i:i+3,:);

                %%normalize 
                
                    tmp(1,:)=tmp(1,:)./A_QT;
                    tmp(2,:)=tmp(2,:)./C_QT;
                    tmp(3,:)=tmp(3,:)./G_QT;
                    tmp(4,:)=tmp(4,:)./T_QT;
                %%pick the maximum
                %%TODO: we need the non maximal to measure the confidence,

                [MaxI,~]=max(tmp,[],1);
              
                
                %{
                qt_tmp=quantile(tmp,0.6,1);
                r_tmp=tmp;
                for j=1:4
                    r_tmp(j,:)=tmp(j,:).*(tmp(j,:)>qt_tmp);
                end
                %}
                
                
                mask=bsxfun(@eq,tmp,MaxI);
                mask=~mask;
                tmp=x(i:i+3,:);
                tmp(mask)=0;
                x(i:i+3,:)=tmp;    
            end
        end
        
        %%% below round is for more accurate support set estimation
        %we shall use the maximums in the above round
        
        if(itr>1)
            %{
            R_max=5;
            [QTs ] = NormV(x',x_tmp',qual',R_max);
            %}
            %if (~FAST)
            if(itr>=1)
                R_max=2;
                [ QTs1 ] = NormV( bsxfun(@rdivide,x',repmat(QTs,1,size(x,1)/4)),bsxfun(@rdivide,x_tmp',repmat(QTs,1,size(x,1)/4)),qual' ,R_max);
                QTs=QTs.*QTs1;
            else
                 R_max=2;
                [ QTs1 ] = NormV( bsxfun(@rdivide,x',repmat(QTs,1,size(x,1)/4)),bsxfun(@rdivide,x_tmp',repmat(QTs,1,size(x,1)/4)),[] ,R_max);
                QTs=QTs.*QTs1;
                
            end
            %end
           %}
            %QTs./max(QTs)%%print again

            A_QT=QTs(1);
            C_QT=QTs(2);
            G_QT=QTs(3);
            T_QT=QTs(4);

        end
   
    else
        if (itr>1 && ~FAST)
            %{
            R_max=5;
            [QTs ] = NormV(x',x_tmp',qual',R_max);
            %}
            %if (~FAST)
                R_max=2;
                for Norm_it=1:1
                [ QTs1 ] = NormV( bsxfun(@rdivide,x',repmat(QTs,1,size(x,1)/4)),bsxfun(@rdivide,x_tmp',repmat(QTs,1,size(x,1)/4)),qual' ,R_max);
                QTs=QTs.*QTs1;
                end
            %end
            %}
            
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
    end
    cycle=0;
    if(itr>1)


        for i=1:4:size_x-3
            %
            cycle=cycle+1;
            tmp=x(i:i+3,:);

            tmp(1,:)=tmp(1,:)./A_QT;
            tmp(2,:)=tmp(2,:)./C_QT;
            tmp(3,:)=tmp(3,:)./G_QT;
            tmp(4,:)=tmp(4,:)./T_QT;

            if(0)%old strategy
               [MaxI,~]=max(tmp,[],1);
                mask=bsxfun(@eq,tmp,MaxI);
                sum_m=0;
                for r=1:4
                    sum_m=sum_m+mask(r,:);
                    mask(r,:)=mask(r,:).*(sum_m<=1);
                end

                mask=~mask;
                
                tmp=x(i:i+3,:);
                tmp(mask)=0;
                x(i:i+3,:)=tmp; 
    
            end
        
            if(1)%new strategy
                qt_tmp=median(tmp,1);%quantile(tmp,0.6,1);
                r_tmp=tmp;
                for j=1:4
                    r_tmp(j,:)=x(i+j-1,:).*(tmp(j,:)>qt_tmp);
                end
                x(i:i+3,:)=r_tmp;
            end
            %}
            
            
        end
    end
    
    %%%%%%%%%end of penalty sub problem
    
    %figure;bar3(reshape(x(:,1),4,size(x,1)/4),'detached');
    %figure;bar3(reshape(b(:,1),4,size(b,1)/4),'detached');
    if(itr<=1)
        %For the first round, no reg is better because reg suppress the
        %intensity in an unwanted way.
        LHS=(K'*K);
        %LHS=(K'*K+.01*eye(length(K)));
        %LHS=(K'*K+.01*diag(repmat(QTs,1,length(K)/4)));
        RHS=(K'*b);
   
    else
       
        LHS=(K'*K+beta*eye(length(K)));
        RHS=(K'*b+beta*x);

    end
    x=LHS\RHS;
    if(itr<=1) 
        x_tmp=x;
    end
    

    
    if itr==1
    %%get quality score
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
    
    

    
    if (LastIt || (max_itr==0) )
        
        %{
        tmp=bsxfun(@rdivide,x,repmat(QTs(:),size_x/4,1));
        cycle=0;
        for i=1:4:size_x-3
            %
            cycle=cycle+1;
            [MaxI]=max(tmp(i:i+3,:),[],1);
            mask=bsxfun(@eq,tmp(i:i+3,:),MaxI);
            mask=~mask;
            tmp2=x(i:i+3,:);
            tmp2(mask)=0;
            x(i:i+3,:)=tmp2;    
        end
        %}
        %
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
            
            mask=bsxfun(@eq,tmp,MaxI);
            mask=~mask;
            tmp=x(i:i+3,:);
            tmp(mask)=0;

            x(i:i+3,:)=tmp;    
        end
        %}
    end

    
    
    beta=beta*kappa;
    

    if beta*kappa>betamax||itr>=max_itr-1
        LastIt=1;
    end
    %beta
    
end





