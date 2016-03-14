function [ K ] = CreateK( WindowSize,K0,CT_m,PrePhaseCycles,PhaseCycles )
%CREATEK Summary of this function goes here
%   B=K*x; x is the ground truth signal;B is the observation
%This process is a convolution, we create the blur pattern first, and we
%flip the matrix finally.

%PrePhaseCycles:cycles for pre-phasing, 
%PhaseCycles:cycles for phasing, 
%CT_m: cross talk matrix.

    K=eye(WindowSize*4);

    for i=1:4
        CT_m(:,i)=CT_m(:,i)./CT_m(i,i);
    end

    CT_AC=CT_m(2,1);
    CT_CA=CT_m(1,2);
    CT_GT=CT_m(4,3);
    CT_TG=CT_m(3,4);
    %%%A channel
    for i=1:4:WindowSize*4
        K(i,i+1)=CT_AC;
        for j=1:PrePhaseCycles
            if i-4*j>0
                PrePhasing_Int=K0(PrePhaseCycles+1-j,1);
                K(i,i-4*j)=PrePhasing_Int;
                K(i,i-4*j+1)=PrePhasing_Int*CT_AC;
            else
                break;
            end
        end

        for j=1:PhaseCycles
            if i+4*j<=(WindowSize)*4
                Phasing_Int=K0(PhaseCycles+1+j,1);
                K(i,i+4*j)=Phasing_Int;
                K(i,i+4*j+1)=Phasing_Int*CT_AC;
            else
                break;
            end
        end        
    end

    %%%C channel

    for i=2:4:WindowSize*4
        K(i,i-1)=CT_CA;

        for j=1:PrePhaseCycles
            if i-4*j>0
                PrePhasing_Int=K0(PrePhaseCycles+1-j,2);
                K(i,i-4*j)=PrePhasing_Int;
                K(i,i-4*j-1)=PrePhasing_Int*CT_CA;
            else
                break;
            end
        end

        for j=1:PhaseCycles
            if i+4*j<=(WindowSize)*4
                Phasing_Int=K0(PhaseCycles+1+j,2);
                K(i,i+4*j)=Phasing_Int;
                K(i,i+4*j-1)=Phasing_Int*CT_CA;
            else
                break;
            end
        end

    end

    %%%G channel

    for i=3:4:WindowSize*4

        K(i,i+1)=CT_GT;

        for j=1:PrePhaseCycles
            if i-4*j>0
                PrePhasing_Int=K0(PrePhaseCycles+1-j,3);
                K(i,i-4*j)=PrePhasing_Int;
                K(i,i-4*j+1)=PrePhasing_Int*CT_GT;
            else
                break;
            end
        end

        for j=1:PhaseCycles
            if i+4*j<=(WindowSize)*4
                Phasing_Int=K0(PhaseCycles+1+j,3);
                K(i,i+4*j)=Phasing_Int;
                K(i,i+4*j+1)=Phasing_Int*CT_GT;
            else
                break;
            end
        end


    end

    %%%T channel

    for i=4:4:WindowSize*4
        K(i,i-1)=CT_TG;


        for j=1:PrePhaseCycles
            if i-4*j>0
                PrePhasing_Int=K0(PrePhaseCycles+1-j,4);
                K(i,i-4*j)=PrePhasing_Int;
                K(i,i-4*j-1)=PrePhasing_Int*CT_TG;
            else
                break;
            end
        end

        for j=1:PhaseCycles
            if i+4*j<=(WindowSize)*4
                Phasing_Int=K0(PhaseCycles+1+j,4);
                K(i,i+4*j)=Phasing_Int;
                K(i,i+4*j-1)=Phasing_Int*CT_TG;
            else
                break;
            end
        end


    end


    %

    for i=1:size(K,1)
        beg_col=mod(i,4);
        if beg_col==0
            beg_col=4;
        end
        sum_k=sum(K(i,beg_col:4:end));
        K(i,:)=K(i,:)./sum_k;
       
    end


    K=K';
    
   

end

