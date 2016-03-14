
%function [ output,output_tmp,Flows_all,qual_all ] = BaseCaller_simple( Flows_all,MaxLen,BatchID )
%%avoid copying of the parameter

%BASECALLER Summary of this function goes here
% This base calling algorithm fixes the cross talk estimation for all cycles.
TrainingSize=5000;    
ReadsPerBatch=10000;
BatchBeg=1;%1
DEBUG=1;

if(~exist('UPDATE_CROSSTALK','var')||isempty(UPDATE_CROSSTALK))
    UPDATE_CROSSTALK=0;
end


if(~exist('deconv_time','var')||isempty(deconv_time))
    deconv_time=0.0;
end
if(~exist('norm_time','var')||isempty(norm_time))
    norm_time=0.0;
end
if(~exist('est_time','var')||isempty(est_time))
    est_time=0.0;
end
if(~exist('final_deconv_time','var')||isempty(final_deconv_time))
    final_deconv_time=0.0;
end
if(~exist('writing_time','var')||isempty(writing_time))
    writing_time=0.0;
end
if(~exist('kernel_th','var')||isempty(kernel_th))
    kernel_th=0.1;
  %hiseq2000:0.1, old data:0.0;
end

if(~exist('CutOffTh','var')||isempty(CutOffTh))
    CutOffTh=5;
end

if(~exist('RAW_QS','var')||isempty(RAW_QS))
    RAW_QS=0;
end

if(~exist('Fastq','var')||isempty(Fastq))
    Fastq=0;
end
if(~exist('N_th','var')||isempty(N_th))
    N_th=0.501;
end

if(~exist('FAST','var')||isempty(FAST))
    FAST=1;
end

%Flows_all(Flows_all<0)=0;
nReads=size(Flows_all,1);

TrainAll=0;
WindowSize0=min(ceil(MaxLen/4/3),21);

%%This is the flow order:
FlowSequence(1:MaxLen)='A';
FlowSequence(2:4:MaxLen)='C';
FlowSequence(3:4:MaxLen)='G';
FlowSequence(4:4:MaxLen)='T';
CrossTalkWinSize=4;%%window size used for estimate the cross talk and the estimates will be fixed for all future cycles.




if Fastq
    fid = fopen(strcat('Batch_',num2str(BatchID),'.fastq'), 'w');
else
    fid = fopen(strcat('Batch_',num2str(BatchID),'.fasta'), 'w');
end

if RAW_QS
    fid2 = fopen(strcat('Batch_',num2str(BatchID),'_rqs.txt'), 'w');
end

%tic


%%%%%%training



kernels={};

%TrainIdx = randsample(nReads,TrainingSize);
TrainIdx=round([(nReads/TrainingSize):(nReads/TrainingSize):nReads-1]);
TrainingSize=length(TrainIdx);
x_crop=double(Flows_all(TrainIdx,:));


WindowSize=WindowSize0;
 
x_crop=x_crop(:,1:MaxLen);

base_order='ACGT';

output=zeros(size(x_crop));

if(DEBUG)
    output_tmp=zeros(size(x_crop));
    qual_all=zeros(size(x_crop,1),size(x_crop,2)/4);
end



%%%%%%%initial values


tTraining=tic;

CT_AC=0.8;
CT_CA=0.25;
CT_GT=0.5;
CT_TG=0.05;

%kernel size
kernel_sz=5;%must be an odd number
%%Phasing:blur in the -> direction; Prephasing: <- direction
PrePhase=floor(kernel_sz/2);
Phase=floor(kernel_sz/2);
deconv_rounds=2;
pos=1;

   
[xx]=[1:kernel_sz];
xx=xx-ceil(kernel_sz/2);
sigma2=0.1;
K0=exp(-(xx.^2).^2/2/sigma2);
K0=K0(:)./sum(K0(:));
K0=repmat(K0,[1,4]);

%CT: cross talk, the blur between channels
CT=[CT_AC,CT_CA,CT_GT,CT_TG];

x0=x_crop;

K_CT=eye(4);
K_CT(1,2)=CT_CA;%how C is spread to A
K_CT(2,1)=CT_AC;%how A is spread to C
K_CT(3,4)=CT_TG;%how T is spread to G
K_CT(4,3)=CT_GT;%how G is spread to T

CT_Flows=x_crop(:,1:CrossTalkWinSize*4)';
K_CT2=K_CT;
for ct_it=1:2

    if ct_it<=1
        [ Flows_CT_Corr,qual,Flows_CT_Corr_0,QT_V] = CrossTalkCorrection_simple( K_CT2,CT_Flows ); 
        [ K_CT2,CT ] = CrossTalkEst(Flows_CT_Corr, CT_Flows);
    else                
       [ Flows_CT_Corr,qual,Flows_CT_Corr_0,QT_V] = CrossTalkCorrection_simple( K_CT2,CT_Flows,QT_V); 
        [ K_CT2,CT ] = CrossTalkEst(Flows_CT_Corr, CT_Flows,qual);
    end
end
%


Flows_CT_Corr=reshape(Flows_CT_Corr,size(CT_Flows));
Flows_CT_Corr_0=reshape(Flows_CT_Corr_0,size(CT_Flows));
qual=reshape(qual,size(CT_Flows,1)/4,size(CT_Flows,2));
%R_max=5;
%[ QT_V ] = NormV( Flows_CT_Corr',Flows_CT_Corr_0',qual' ,R_max);
%
R_max=2;
for Norm_it=1:1
[ QT_V1 ] = NormV( bsxfun(@rdivide,Flows_CT_Corr',repmat(QT_V,1,size(Flows_CT_Corr,1)/4)),bsxfun(@rdivide,Flows_CT_Corr_0',repmat(QT_V,1,size(Flows_CT_Corr,1)/4)),qual' ,R_max);
QT_V=QT_V.*QT_V1;
end
%}
K_CT=K_CT2;  
QTs=[];
QTs=QT_V;

n_kernels=0;
while pos<MaxLen-1%-WindowSize*4
    if pos > (MaxLen-(WindowSize)*4)
        pos=MaxLen-(WindowSize+kernel_sz)*4+1;
        WindowSize=(WindowSize+kernel_sz);
    else
        if(pos>(kernel_sz*4))
            pos=pos-(kernel_sz)*4;
        end
    end
    %fprintf('Cycle %d.\r\n',ceil(pos/4));
    FlowValues=x_crop(:,pos:pos-1+WindowSize*4)';

    n_kernels=n_kernels+1;

    [ K ] = CreateK( WindowSize,K0,K_CT,PrePhase,Phase );

    kernels{n_kernels}.K=K0;
    kernels{n_kernels}.pos=pos;
    kernels{n_kernels}.WindowSize=WindowSize;
    kernels{n_kernels}.CT=CT;
    kernels{n_kernels}.QT_V=QTs;
    kernels{n_kernels}.PrePhase=PrePhase;
    kernels{n_kernels}.Phase=Phase;

    lambda=0.01;

    [FlowValuesC1,qual,FlowValuesC1_tmp,QTs]= Deconv_ISD_simple(K,FlowValues,lambda,0,QTs);
    R_max=2;
    for Norm_it=1:1
    [ QTs1 ] = NormV( bsxfun(@rdivide,FlowValuesC1',repmat(QTs,1,size(FlowValuesC1,1)/4)),bsxfun(@rdivide,FlowValuesC1_tmp',repmat(QTs,1,size(FlowValuesC1,1)/4)),qual' ,R_max);
    QTs=QTs.*QTs1;
    end 

    for r=1:deconv_rounds-1
        if(~UPDATE_CROSSTALK)
        [ K0 ] =KernelEst(FlowValuesC1,FlowValues,qual,CT,PrePhase,Phase,kernel_th);
        else
        [ K0,CT,K_CT ] =KernelEst(FlowValuesC1,FlowValues,qual,CT,PrePhase,Phase,kernel_th);

        CT
        K_CT

        end

        [ K ] = CreateK( WindowSize,K0,K_CT,PrePhase,Phase );

        kernels{n_kernels}.K=K;
        if(exist('K0_CT','var'))
            kernels{n_kernels}.K0_CT=K0_CT;
        end
        kernels{n_kernels}.K0=K0;
        kernels{n_kernels}.CT=CT;
        kernels{n_kernels}.pos=pos;
        kernels{n_kernels}.WindowSize=WindowSize;
        kernels{n_kernels}.PrePhase=PrePhase;
        kernels{n_kernels}.Phase=Phase;

        if(r==deconv_rounds-1)
            Final=0;
        else
            Final=0;
        end
        [FlowValuesC1 ,qual,FlowValuesC1_tmp,QTs]= Deconv_ISD_simple(K,FlowValues,lambda,Final,QTs);
        R_max=2;
        for Norm_it=1:1
        [ QTs1 ] = NormV( bsxfun(@rdivide,FlowValuesC1',repmat(QTs,1,size(FlowValuesC1,1)/4)),bsxfun(@rdivide,FlowValuesC1_tmp',repmat(QTs,1,size(FlowValuesC1,1)/4)),qual' ,R_max);
        QTs=QTs.*QTs1;
        end
        %}
        kernels{n_kernels}.QT_V=QTs;


    end


     if(pos==1)
        output(:,pos:pos+WindowSize*4-1)=FlowValuesC1';
        if DEBUG 
            output_tmp(:,pos:pos+WindowSize*4-1)=FlowValuesC1_tmp';
            qual_all(:,pos:pos+WindowSize-1)=qual';
        end
    else
        offset=floor(kernel_sz/2);

        output(:,pos+offset*4:pos+WindowSize*4-1)=FlowValuesC1((1+offset*4):end,:)'; 
        if DEBUG            
            output_tmp(:,pos+offset*4:pos+WindowSize*4-1)=FlowValuesC1_tmp((1+offset*4):end,:)';
            qual_all(:,(pos+3)/4+offset:(pos+3)/4+WindowSize-1)=qual((1+offset):end,:)';
        end

    end
    pos=pos+(WindowSize)*4;


end            


Training_Elapsed=toc(tTraining);

fprintf(strcat('Training time: ',num2str(Training_Elapsed),' sec.\r\n'));



%%%%%%end of training




batch=BatchBeg-1;
out_r=0;




for i=1+ReadsPerBatch*(BatchBeg-1):ReadsPerBatch:nReads
    batch=batch+1;
    if(i+ReadsPerBatch-1<=nReads-ReadsPerBatch-1)
        x_crop=double(Flows_all(i:i+ReadsPerBatch-1,:));
    else
        if(i>nReads)
            break;
        end
        x_crop=double(Flows_all(i:end,:)); 
        nReads=nReads-ReadsPerBatch-1;
    end
    if(size(x_crop,1)<ReadsPerBatch)
        break;
    end
    if (mod(batch,10)==1)
        fprintf(strcat('Batch',num2str(batch),'\n'));
    end
    WindowSize=WindowSize0;
 
    x_crop=x_crop(:,1:MaxLen);

    %%for debugging
    %FlowStemPlot(FlowSequence,x_crop(100,:))
    %figure;bar3(reshape(x_crop(100,:),4,size(x_crop,2)/4))

    output=zeros(size(x_crop));
    
    if(DEBUG)
        output_tmp=zeros(size(x_crop));
        qual_all=zeros(size(x_crop,1),size(x_crop,2)/4);
    end
    %%%%%%%initial values

    pos=1;

   
    if (TrainAll==1)
        [xx]=[1:kernel_sz];
        xx=xx-ceil(kernel_sz/2);
        sigma2=0.1;
        K0=exp(-(xx.^2).^2/2/sigma2);
        K0=K0(:)./sum(K0(:));
        K0=repmat(K0,[1,4]);

        %CT: cross talk, the blur between channels
        CT=[CT_AC,CT_CA,CT_GT,CT_TG];

        x0=x_crop;
        
        K_CT=eye(4);
        K_CT(1,2)=CT_CA;%how C is spread to A
        K_CT(2,1)=CT_AC;%how A is spread to C
        K_CT(3,4)=CT_TG;%how T is spread to G
        K_CT(4,3)=CT_GT;%how G is spread to T

        CT_Flows=x_crop(:,1:CrossTalkWinSize*4)';
        K_CT2=K_CT;
        for ct_it=1:2
            
            if ct_it<=1
                [ Flows_CT_Corr,qual,Flows_CT_Corr_0,QT_V] = CrossTalkCorrection_simple( K_CT2,CT_Flows ); 
                [ K_CT2,CT ] = CrossTalkEst(Flows_CT_Corr, CT_Flows);
            else                
               [ Flows_CT_Corr,qual,Flows_CT_Corr_0,QT_V] = CrossTalkCorrection_simple( K_CT2,CT_Flows,QT_V); 
                [ K_CT2,CT ] = CrossTalkEst(Flows_CT_Corr, CT_Flows,qual);
            end
        end
        %
        
        
        Flows_CT_Corr=reshape(Flows_CT_Corr,size(CT_Flows));
        Flows_CT_Corr_0=reshape(Flows_CT_Corr_0,size(CT_Flows));
        qual=reshape(qual,size(CT_Flows,1)/4,size(CT_Flows,2));
        %R_max=5;
        %[ QT_V ] = NormV( Flows_CT_Corr',Flows_CT_Corr_0',qual' ,R_max);
        %
        R_max=2;
        for Norm_it=1:1
        [ QT_V1 ] = NormV( bsxfun(@rdivide,Flows_CT_Corr',repmat(QT_V,1,size(Flows_CT_Corr,1)/4)),bsxfun(@rdivide,Flows_CT_Corr_0',repmat(QT_V,1,size(Flows_CT_Corr,1)/4)),qual' ,R_max);
        QT_V=QT_V.*QT_V1;
        end
        %}
        K_CT=K_CT2;
    end
    
    QTs=[];
    QTs=QT_V;
    
    n_kernels=0;
    while pos<MaxLen-1%-WindowSize*4
        if pos > (MaxLen-(WindowSize)*4)
            pos=MaxLen-(WindowSize+kernel_sz)*4+1;
            WindowSize=(WindowSize+kernel_sz);
        else
            if(pos>(kernel_sz*4))
                pos=pos-(kernel_sz)*4;
            end
        end
        %fprintf('Cycle %d.\r\n',ceil(pos/4));
        FlowValues=x_crop(:,pos:pos-1+WindowSize*4)';
        
        n_kernels=n_kernels+1;
        if (TrainAll==1)
            
            [ K ] = CreateK( WindowSize,K0,K_CT,PrePhase,Phase );

            kernels{n_kernels}.K=K0;
            kernels{n_kernels}.pos=pos;
            kernels{n_kernels}.WindowSize=WindowSize;
            kernels{n_kernels}.CT=CT;
            kernels{n_kernels}.QT_V=QTs;
            kernels{n_kernels}.PrePhase=PrePhase;
            kernels{n_kernels}.Phase=Phase;

            lambda=0.01;

            t_d=tic;
            [FlowValuesC1,qual,FlowValuesC1_tmp,QTs]= Deconv_ISD_simple(K,FlowValues,lambda,0,QTs);
            %{
            R_max=5;
            [ QTs ] = NormV( FlowValuesC1',FlowValuesC1_tmp',qual',R_max );
            %}
            R_max=2;
            for Norm_it=1:1
            [ QTs1 ] = NormV( bsxfun(@rdivide,FlowValuesC1',repmat(QTs,1,size(FlowValuesC1,1)/4)),bsxfun(@rdivide,FlowValuesC1_tmp',repmat(QTs,1,size(FlowValuesC1,1)/4)),qual' ,R_max);
            QTs=QTs.*QTs1;
            end 
            
            %}

            td2=toc(t_d);
            deconv_time=deconv_time+td2;
    %        change these parameters if necessary
    %{
            if(pos>100)
                PrePhase=2;
                Phase=2;
                deconv_rounds=2;
            end
            if(pos>200)

                PrePhase=2;
                Phase=2;
                deconv_rounds=2;

            end
            if(pos>300)
                PrePhase=2;
                Phase=2;
                deconv_rounds=2;
            end
%}
            for r=1:deconv_rounds-1
                t_e=tic;
                if(~UPDATE_CROSSTALK)
                [ K0 ] =KernelEst(FlowValuesC1,FlowValues,qual,CT,PrePhase,Phase,kernel_th);
                else
                [ K0,CT,K_CT ] =KernelEst(FlowValuesC1,FlowValues,qual,CT,PrePhase,Phase,kernel_th);
                    
                end
                CT
                K_CT
                te2=toc(t_e);
                est_time=est_time+te2;
                
                [ K ] = CreateK( WindowSize,K0,K_CT,PrePhase,Phase );
                
                kernels{n_kernels}.K=K;
                if(exist('K0_CT','var'))
                    kernels{n_kernels}.K0_CT=K0_CT;
                end
                kernels{n_kernels}.K0=K0;
                kernels{n_kernels}.CT=CT;
                kernels{n_kernels}.pos=pos;
                kernels{n_kernels}.WindowSize=WindowSize;
                kernels{n_kernels}.PrePhase=PrePhase;
                kernels{n_kernels}.Phase=Phase;
                
                if(r==deconv_rounds-1)
                    Final=0;
                else
                    Final=0;
                end
                t_d=tic;
                %[FlowValuesC1 ,qual,FlowValuesC1_tmp,NormFactor]= Deconv_ISD2(K,FlowValues,lambda,Final);
                [FlowValuesC1 ,qual,FlowValuesC1_tmp,QTs]= Deconv_ISD_simple(K,FlowValues,lambda,Final,QTs);
                %{ 
                R_max=5;
                [ QTs ] = NormV( FlowValuesC1',FlowValuesC1_tmp',qual',R_max );
                %}
                R_max=2;
                for Norm_it=1:1
                [ QTs1 ] = NormV( bsxfun(@rdivide,FlowValuesC1',repmat(QTs,1,size(FlowValuesC1,1)/4)),bsxfun(@rdivide,FlowValuesC1_tmp',repmat(QTs,1,size(FlowValuesC1,1)/4)),qual' ,R_max);
                QTs=QTs.*QTs1;
                end
                %}
                kernels{n_kernels}.QT_V=QTs;
                
                td2=toc(t_d);
                deconv_time=deconv_time+td2;
                if Final
                final_deconv_time=final_deconv_time+td2;
                end
            end
        
        else
            
            
            K=kernels{n_kernels}.K;
            K0=kernels{n_kernels}.K0;
            QTs=kernels{n_kernels}.QT_V;
%            K0_CT=kernels{n_kernels}.K0_CT;
            CT=kernels{n_kernels}.CT;
            PrePhase=kernels{n_kernels}.PrePhase;
            Phase=kernels{n_kernels}.Phase;
            Final=1;

            t_d=tic;
            [FlowValuesC1 ,qual,FlowValuesC1_tmp,QTs]= Deconv_ISD_simple(K,FlowValues,lambda,Final,QTs,FAST);

            td2=toc(t_d);
            deconv_time=deconv_time+td2;
            if Final
            final_deconv_time=final_deconv_time+td2;
            end
            
        end
        if(pos==1)
            output(:,pos:pos+WindowSize*4-1)=FlowValuesC1';
            if DEBUG 
                output_tmp(:,pos:pos+WindowSize*4-1)=FlowValuesC1_tmp';
                qual_all(:,pos:pos+WindowSize-1)=qual';
            end
        else
            offset=floor(kernel_sz/2);
            
            output(:,pos+offset*4:pos+WindowSize*4-1)=FlowValuesC1((1+offset*4):end,:)'; 
            if DEBUG            
                output_tmp(:,pos+offset*4:pos+WindowSize*4-1)=FlowValuesC1_tmp((1+offset*4):end,:)';
                qual_all(:,(pos+3)/4+offset:(pos+3)/4+WindowSize-1)=qual((1+offset):end,:)';
            end
            
        end
        pos=pos+(WindowSize)*4;

    end


    save('kernels.mat','kernels');
    %save('output.mat','output','output_tmp','qual_all');

    
    %fprintf('Output sequences.\r\n');
    if RAW_QS
        qual_all(qual_all<0.5)=0.5;
        qual_all(qual_all>1.0)=1.0;
    end
    qs_format=repmat('%f ',[1,MaxLen/4]);
    qs_format(end)=[];

    seqs=repmat('N',[size(output,1),MaxLen/4]);
    tmp=output(:,1:MaxLen);
    tmp2=x_crop(:,1:MaxLen);
    tmp(tmp<CutOffTh)=0;
    tmp(tmp2<CutOffTh)=0;
    cycle=0;
    for j=1:4:MaxLen
        cycle=cycle+1;
        %channels=tmp(:,j:j+3);
        %[MaxIntensity,max_idx]=max(channels,[],2);
        [MaxIntensity,max_idx]=max(tmp(:,j:j+3),[],2);
        seqs(:,cycle)=base_order(max_idx);
        seqs(MaxIntensity<CutOffTh,cycle)='N';
    end

    qual_all(qual_all<0.5)=0.5;
    qual_all(qual_all>1.0)=1.0;

    qs_mask=(qual_all<N_th);
    seqs(qs_mask)='N';
    
    tWrite=tic();
    if Fastq
        max_qs=max(qual_all(:));
        min_qs=min(qual_all(:));
        qual_all=(qual_all-min_qs)/(max_qs-min_qs);
        qual_all=floor(qual_all*length(QC))+1;
        qual_all(qual_all>=length(QC))=length(QC);
        qs_c=QC(qual_all);
        
        for r=1:size(output,1)
            out_r=out_r+1;
            %{
            fprintf(fid,'@seq_%d\r\n',out_r);
            fprintf(fid,seqs(r,:));
            fprintf(fid,'\r\n');
            fprintf(fid,'+\r\n');
            fprintf(fid,qs_c(r,:));
            fprintf(fid,'\r\n');
            %}
            fprintf(fid,'@seq_%d\r\n%s\r\n+\r\n%s\r\n',out_r,seqs(r,:),qs_c(r,:));
        end
    end
        
    if ~Fastq

        if RAW_QS
            for r=1:size(output,1)
                out_r=out_r+1;
                %{
                fprintf(fid,'>seq_%d\r\n',out_r);
                fprintf(fid,seqs(r,:));
                fprintf(fid,'\r\n');
                %}
                fprintf(fid,'>seq_%d\r\n%s\r\n',out_r,seqs(r,:));
                
                fprintf(fid2,'%d\r\n',MaxLen/4);
                fprintf(fid2,qs_format,qual_all(r,:));
                fprintf(fid2,'\r\n');
                
                %}
                
            end
        else
            
            for r=1:size(output,1)
                out_r=out_r+1;
                 %{
                fprintf(fid,'>seq_%d\r\n',out_r);
                fprintf(fid,seqs(r,:));
                fprintf(fid,'\r\n');
                %}
                
                fprintf(fid,'>seq_%d\r\n%s\r\n',out_r,seqs(r,:));
                
            end
        end

    end

    tWrite2=toc(tWrite);
    writing_time=writing_time+tWrite2;


end

out_r
if(out_r~=totClust)
   [out_r,totClust]  
end

fclose(fid);

fprintf(strcat('Deconv time: ',num2str(deconv_time),' sec.\r\n'));
fprintf(strcat('Final deconv time: ',num2str(final_deconv_time),' sec.\r\n'));
fprintf(strcat('Blur est time: ',num2str(est_time),' sec.\r\n'));
fprintf(strcat('Writinging time: ',num2str(writing_time),' sec.\r\n'));
