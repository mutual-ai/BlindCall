clear all;
%LASTN = maxNumCompThreads(1)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%modify the directory below
CIF_dir='D:\TestPrograms\BaseCallingData\PhiX174_UMD_HiSeq_201305\Data\Intensities\L004';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%POS_dir='D:\TestPrograms\BaseCallingData\PhiX174_100nt_pe_UMD_HiSeq\Data\Intensities\';
N_th=0.501;
kernel_th=0.13;
Tiles_beg=1;
Tiles_end=1;
totCycles=101;
%totClust=100000;
tmp_dircs = dir(CIF_dir);

cif_dircs={};
cnt=0;
for i=1:length(tmp_dircs)
   if tmp_dircs(i).isdir==1&&(~strcmp(tmp_dircs(i).name,'.'))&&(~strcmp(tmp_dircs(i).name,'..'))
      cnt=cnt+1;
      cif_dircs{cnt} =tmp_dircs(i).name; 
   end
end


cif_files=dir(strcat(CIF_dir,'/',cif_dircs{1},'/*.cif'));

if ~exist('Tiles_beg','var')||isempty(Tiles_beg)
   Tiles_beg=1; 
end
Tiles_beg=min(Tiles_beg,length(cif_files));

if ~exist('Tiles_end','var')||isempty(Tiles_end)
   Tiles_end=length(cif_files); 
end
Tiles_end=min(Tiles_end,length(cif_files));


%pos_files=dir(strcat(POS_dir,'/*.txt'));


if ~exist('totCycles','var')||isempty(totCycles)
   totCycles=length(cif_dircs); 
end
totCycles=min(totCycles,length(cif_dircs));

tRead=tic;

for i=Tiles_beg:Tiles_end
    cif_dirc=strcat(CIF_dir,'/',cif_dircs{1});
    cif_name=cif_files(i).name;
    fid=fopen(strcat(cif_dirc,'/',cif_name));
    CIF = fread(fid, 3,'uint8');
    Version=fread(fid, 1,'uint8');
    Precision=fread(fid, 1,'uint8');
    Cycle=fread(fid, 1,'uint16');
    nCycle=fread(fid, 1,'uint16');
    nClust=fread(fid, 1,'uint32');
    if(exist('totClust','var')&&~isempty(totClust))
        totClust=min(totClust,nClust);
    else
        totClust=nClust;
    end
    fclose(fid);
    Flows_all=zeros(totClust,4*totCycles,'uint16');
    %pos_name=pos_files(i).name;
    %fid2=fopen(strcat(POS_dir,'/',pos_name));
    %cod= fscanf(fid2, '%f');
    %cod=reshape(cod,2,length(cod)/2);
    %fclose(fid2);
    for j=1:totCycles
        
        fprintf(strcat('Tile ',num2str(i),', ','cycle ',num2str(j),'\n'));
        cif_dirc=strcat(CIF_dir,'/C',num2str(j),'.1');
        cif_name=cif_files(i).name;
        fid=fopen(strcat(cif_dirc,'/',cif_name));
        CIF = fread(fid, 3,'uint8');
        Version=fread(fid, 1,'uint8');
        Precision=fread(fid, 1,'uint8');
        Cycle=fread(fid, 1,'uint16');
        nCycles=fread(fid, 1,'uint16');
        nClust=fread(fid, 1,'uint32');

       %
        if(Cycle~=j)
            fprintf('Cycle error.');
            return;
        end
        
        if Precision==2
            flows=fread(fid,nClust*4*nCycles,'int16');
            
        end
         fclose(fid);
        flows=reshape(flows,nClust,4);
        Flows_all(1:totClust,1+4*(j-1):4*j)=flows(1:totClust,:);

    end
    
    
   out_str=strcat('Tile',num2str(i),'.mat');
   % save(out_str,'Flows_all','-v7.3');
   CallLen=min(totCycles,size(Flows_all,2)/4);
   LoadingTime=toc(tRead)

     
   BatchID=i;
   MaxLen=CallLen*4;
   BaseCaller_simple;
end
