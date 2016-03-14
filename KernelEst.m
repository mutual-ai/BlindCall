function [ K,CT,K_CT] = KernelEst( A,b,qual,CT,PrePhase,Phase,thresh )
%BLUREST Summary of this function goes here
%   Kernel estimation algorithm, with a penalty term.

if(~exist('thresh','var')||isempty(thresh))
    thresh=0.1; %thresh=0.05;%%%penalty threshold, This current implementation is equivalent to L1 norm thresholding. 
  %hiseq2000:0.1, old data:0.0;
end

%%select the high/low quality portion to estimate the blur.
QS_QT=0.5;
mean_qual=min(qual,[],1);
%mean_qual=min(qual,1);
QS_TH=quantile(mean_qual,QS_QT);
Mask=mean_qual<QS_TH;
A=A(:,Mask);
b=b(:,Mask);

thresh_supp=1e-1;

size_x=size(b,1);
size_y=size(b,2);
%b(b<0)=0;
%A(A<0)=0;
size_k=(PrePhase+Phase+1);

A_Corr=zeros(size_k,size_k);
C_Corr=zeros(size_k,size_k);
G_Corr=zeros(size_k,size_k);
T_Corr=zeros(size_k,size_k);


A_Rows=A(1:4:size_x,:);
C_Rows=A(2:4:size_x,:);
G_Rows=A(3:4:size_x,:);
T_Rows=A(4:4:size_x,:);
b_A_Rows=b(1:4:size_x,:);
b_C_Rows=b(2:4:size_x,:);
b_G_Rows=b(3:4:size_x,:);
b_T_Rows=b(4:4:size_x,:);

%%The below needs to be changed to work for any intensity
%%the following strategy is taken since deconv_ISD keeps only one channel.

Kmask=ones(PrePhase+PrePhase+1,1);

A_Idx=(A_Rows>thresh_supp);
C_Idx=(C_Rows>thresh_supp);
G_Idx=(G_Rows>thresh_supp);
T_Idx=(T_Rows>thresh_supp);



CT_AC=CT(1);
CT_CA=CT(2);
CT_GT=CT(3);
CT_TG=CT(4);

%
%%unmasked
A_b_Rows=b_A_Rows;
C_b_Rows=b_C_Rows;
G_b_Rows=b_G_Rows;
T_b_Rows=b_T_Rows;

%}


 
%%get more reliable estimates
%make selections of the reference set and the observation set
%This part is suggested by Aleksey
%{
tmp1=(A_Idx(1:end-1,:)&A_Idx(2:end,:));

A_Idx(1:end-1,:)=A_Idx(1:end-1,:)&(~tmp1);
A_Idx(2:end,:)=A_Idx(2:end,:)&(~tmp1);

tmp1=(C_Idx(1:end-1,:)&C_Idx(2:end,:));
C_Idx(1:end-1,:)=C_Idx(1:end-1,:)&(~tmp1);
C_Idx(2:end,:)=C_Idx(2:end,:)&(~tmp1);


tmp1=(G_Idx(1:end-1,:)&G_Idx(2:end,:));
G_Idx(1:end-1,:)=G_Idx(1:end-1,:)&(~tmp1);
G_Idx(2:end,:)=G_Idx(2:end,:)&(~tmp1);


tmp1=(T_Idx(1:end-1,:)&T_Idx(2:end,:));
T_Idx(1:end-1,:)=T_Idx(1:end-1,:)&(~tmp1);
T_Idx(2:end,:)=T_Idx(2:end,:)&(~tmp1);


tmp2=(A_Idx(1:end-2,:)&A_Idx(3:end,:));
A_Idx(1:end-2,:)=A_Idx(1:end-2,:)&(~tmp2);
A_Idx(3:end,:)=A_Idx(3:end,:)&(~tmp2);


tmp2=(C_Idx(1:end-2,:)&C_Idx(3:end,:));
C_Idx(1:end-2,:)=C_Idx(1:end-2,:)&(~tmp2);
C_Idx(3:end,:)=C_Idx(3:end,:)&(~tmp2);

tmp2=(G_Idx(1:end-2,:)&G_Idx(3:end,:));
G_Idx(1:end-2,:)=G_Idx(1:end-2,:)&(~tmp2);
G_Idx(3:end,:)=G_Idx(3:end,:)&(~tmp2);

tmp2=(T_Idx(1:end-2,:)&T_Idx(3:end,:));
T_Idx(1:end-2,:)=T_Idx(1:end-2,:)&(~tmp2);
T_Idx(3:end,:)=T_Idx(3:end,:)&(~tmp2);



A_mask=conv2(double(A_Idx),Kmask,'same');
A_mask=A_mask>0.001;
C_mask=conv2(double(C_Idx),Kmask,'same');
C_mask=C_mask>0.001;
G_mask=conv2(double(G_Idx),Kmask,'same');
G_mask=G_mask>0.001;
T_mask=conv2(double(T_Idx),Kmask,'same');
T_mask=T_mask>0.001;
A_b_Rows=b_A_Rows.*A_mask;
C_b_Rows=b_C_Rows.*C_mask;
G_b_Rows=b_G_Rows.*G_mask;
T_b_Rows=b_T_Rows.*T_mask;
A_Rows=A(1:4:size_x,:).*A_Idx;
C_Rows=A(2:4:size_x,:).*C_Idx;
G_Rows=A(3:4:size_x,:).*G_Idx;
T_Rows=A(4:4:size_x,:).*T_Idx;
%}


Atb=zeros(size_k,0);
Ctb=zeros(size_k,0);
Gtb=zeros(size_k,0);
Ttb=zeros(size_k,0);
A_Corr=zeros(size_k,size_k);
C_Corr=zeros(size_k,size_k);
G_Corr=zeros(size_k,size_k);
T_Corr=zeros(size_k,size_k);


%min_k : ||Ak-b||^2
% (A^t)(Ak-b)=0
% AtA k=Atb


%%in the following A_corr is AtA as above, in the calculations we take into
%%consideration the cross talk effects/blur across channels. Better
%%double checked, I spend *a lot of time* to write the following code.


A_Rows=[zeros(PrePhase,size_y);A_Rows;zeros(Phase,size_y)];
C_Rows=[zeros(PrePhase,size_y);C_Rows;zeros(Phase,size_y)];
G_Rows=[zeros(PrePhase,size_y);G_Rows;zeros(Phase,size_y)];
T_Rows=[zeros(PrePhase,size_y);T_Rows;zeros(Phase,size_y)];


for i=1:size_k    
    Atb(i,1)=0;
    Ctb(i,1)=0;
    Gtb(i,1)=0;
    Ttb(i,1)=0;

    %write down the matrices to understand the calculation

    %%since we generate the blur we can hope to get more accurate
    %%calculations here by tapering the blurry observations to deal
    %%with the boundary effects, like the commentted part

    %{

    if(exist('A_general_blurred','var'))
        A_Rows=[zeros(PrePhase,size_y);A_Rows;zeros(Phase,size_y)];
        C_Rows=[zeros(PrePhase,size_y);C_Rows;zeros(Phase,size_y)];
        G_Rows=[zeros(PrePhase,size_y);G_Rows;zeros(Phase,size_y)];
        T_Rows=[zeros(PrePhase,size_y);T_Rows;zeros(Phase,size_y)];
        A_b_Rows=[A_general_blurred(1:4:PrePhase*4,:);b_A_Rows;A_general_blurred(end-4*Phase+1:4:end,:)];
        C_b_Rows=[A_general_blurred(2:4:PrePhase*4,:);b_C_Rows;A_general_blurred(end-4*Phase+2:4:end,:)];
        G_b_Rows=[A_general_blurred(3:4:PrePhase*4,:);b_G_Rows;A_general_blurred(end-4*Phase+3:4:end,:)];
        T_b_Rows=[A_general_blurred(4:4:PrePhase*4,:);b_T_Rows;A_general_blurred(end-4*Phase+4:4:end,:)];

    end
    %}
    for k=1:size(A_b_Rows,1)
        Atb(i,1)=Atb(i,1)+sum(A_Rows(k+i-1,:).*A_b_Rows(k,:));
        Atb(i,1)=Atb(i,1)+CT_CA*sum(C_Rows(k+i-1,:).*A_b_Rows(k,:));
        Ctb(i,1)=Ctb(i,1)+sum(C_Rows(k+i-1,:).*C_b_Rows(k,:));
        Ctb(i,1)=Ctb(i,1)+CT_AC*sum(A_Rows(k+i-1,:).*C_b_Rows(k,:));
        Gtb(i,1)=Gtb(i,1)+sum(G_Rows(k+i-1,:).*G_b_Rows(k,:));
        Gtb(i,1)=Gtb(i,1)+CT_TG*sum(T_Rows(k+i-1,:).*G_b_Rows(k,:));
        Ttb(i,1)=Ttb(i,1)+sum(T_Rows(k+i-1,:).*T_b_Rows(k,:));
        Ttb(i,1)=Ttb(i,1)+CT_GT*sum(G_Rows(k+i-1,:).*T_b_Rows(k,:));
    end

    for j=1:size_k
        if j<i
            A_Corr(i,j)=A_Corr(j,i);
            C_Corr(i,j)=C_Corr(j,i);
            G_Corr(i,j)=G_Corr(j,i);
            T_Corr(i,j)=T_Corr(j,i);

        else
            A_Corr(i,j)=0;
            C_Corr(i,j)=0;
            G_Corr(i,j)=0;
            T_Corr(i,j)=0;

            %
            for k=j:size(A_Rows,1)
                A_Corr(i,j)=A_Corr(i,j)+sum(A_Rows(k-j+i,:).*A_Rows(k,:));
                A_Corr(i,j)=A_Corr(i,j)+CT_CA.*CT_CA.*sum(C_Rows(k-j+i,:).*C_Rows(k,:));
                C_Corr(i,j)=C_Corr(i,j)+sum(C_Rows(k-j+i,:).*C_Rows(k,:));
                C_Corr(i,j)=C_Corr(i,j)+CT_AC.*CT_AC.*sum(A_Rows(k-j+i,:).*A_Rows(k,:));
                G_Corr(i,j)=G_Corr(i,j)+sum(G_Rows(k-j+i,:).*G_Rows(k,:));
                G_Corr(i,j)=G_Corr(i,j)+CT_TG.*CT_TG.*sum(T_Rows(k-j+i,:).*T_Rows(k,:));
                T_Corr(i,j)=T_Corr(i,j)+sum(T_Rows(k-j+i,:).*T_Rows(k,:));
                T_Corr(i,j)=T_Corr(i,j)+CT_GT.*CT_GT.*sum(G_Rows(k-j+i,:).*G_Rows(k,:));
            end
            %}
            %{
            for k=i:size(A_Rows,1)-size_k+i
                A_Corr(i,j)=A_Corr(i,j)+sum(A_Rows(k,:).*A_Rows(k+j-i,:));
                A_Corr(i,j)=A_Corr(i,j)+CT_CA.*CT_CA.*sum(C_Rows(k,:).*C_Rows(k+j-i,:));
                C_Corr(i,j)=C_Corr(i,j)+sum(C_Rows(k,:).*C_Rows(k+j-i,:));
                C_Corr(i,j)=C_Corr(i,j)+CT_AC.*CT_AC.*sum(A_Rows(k,:).*A_Rows(k+j-i,:));
                G_Corr(i,j)=G_Corr(i,j)+sum(G_Rows(k,:).*G_Rows(k+j-i,:));
                G_Corr(i,j)=G_Corr(i,j)+CT_TG.*CT_TG.*sum(T_Rows(k,:).*T_Rows(k+j-i,:));
                T_Corr(i,j)=T_Corr(i,j)+sum(T_Rows(k,:).*T_Rows(k+j-i,:));
                T_Corr(i,j)=T_Corr(i,j)+CT_GT.*CT_GT.*sum(G_Rows(k,:).*G_Rows(k+j-i,:));
                
            end

            %}
        end
    end
end

%Now we calculate the phasing blur of AC channel.

AC_LHS=A_Corr+C_Corr;
AC_RHS=Atb+Ctb;
AC_K=AC_LHS\AC_RHS;
AC_K(AC_K<0)=0;
AC_K=AC_K./sum(AC_K);

%%GT channel

GT_LHS=G_Corr+T_Corr;
GT_RHS=Gtb+Ttb;
GT_K=GT_LHS\GT_RHS;
GT_K(GT_K<0)=0;
GT_K=GT_K./sum(GT_K);

K=[AC_K,AC_K,GT_K,GT_K];



%%consider some penalty here? L1 L0,ISD,L2
    %L1 is taken here
%%It is beneficial to add this term
K=K-thresh*max(K(:));

K(K<0)=0;
for c=1:size(K,2)       
    K(:,c)=K(:,c)./sum(K(:,c));
end






 CT=eye(4);
[K_m]=CreateK(size_x/4,K,CT,PrePhase,Phase );
A_blurred=K_m*A;

A_Rows=A(1:4:size_x,:);
C_Rows=A(2:4:size_x,:);
G_Rows=A(3:4:size_x,:);
T_Rows=A(4:4:size_x,:);

A_Rows_blurred=A_blurred(1:4:size_x,:);
C_Rows_blurred=A_blurred(2:4:size_x,:);
G_Rows_blurred=A_blurred(3:4:size_x,:);
T_Rows_blurred=A_blurred(4:4:size_x,:);

b_A_Rows=b(1:4:size_x,:);
b_C_Rows=b(2:4:size_x,:);
b_G_Rows=b(3:4:size_x,:);
b_T_Rows=b(4:4:size_x,:);

A_Idx=(A_Rows>1e2);
C_Idx=(C_Rows>1e2);
G_Idx=(G_Rows>1e2);
T_Idx=(T_Rows>1e2);


CT_AC=median(b_C_Rows(A_Idx)./(A_Rows_blurred(A_Idx)+eps));
CT_CA=median(b_A_Rows(C_Idx)./(C_Rows_blurred(C_Idx)+eps));
CT_GT=median(b_T_Rows(G_Idx)./(G_Rows_blurred(G_Idx)+eps));
CT_TG=median(b_G_Rows(T_Idx)./(T_Rows_blurred(T_Idx)+eps));

%}
CT=[CT_AC,CT_CA,CT_GT,CT_TG];

K_CT=eye(4);

K_CT(2,1)=CT(1);
K_CT(1,2)=CT(2);
K_CT(4,3)=CT(3);
K_CT(3,4)=CT(4);





end

