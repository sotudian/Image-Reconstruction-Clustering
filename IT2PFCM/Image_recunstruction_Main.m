
clc;
clear all;
close all;

% ----------------------------------------------------------------------
I1= imread('image.jpg');
s=size(I1);
subplot(121)
imhist(I1(:,:,1))

I2=reshape(I1(:,:,1),[s(1)*s(2),1]);
I3=reshape(I1(:,:,2),[s(1)*s(2),1]);
I4=reshape(I1(:,:,3),[s(1)*s(2),1]);
I5=cat(2,I2,I3,I4);
I6 = im2double(I5);

% scatter3(I6(:,1),I6(:,2),I6(:,3),'.')
% I7=1-I6;  
I7=I6;
ncolors=8;
nC = ncolors;
Xin=I7;
% ----------------------------------------------------------------------
% Options
m1 = 2.0;
m2 = 4.0;
eta1 = 3.0;
eta2 = 5.0;
Cf=0.5;
Cp=0.5;


% ----------------------------------------------------------------------
% Initial weights for PFCM
[V,U] = fcm(Xin,nC);
K=1;
w = FindWeights (Xin, U, V, mean(m1,m2), K);
% w =zeros(size(w));
% ----------------------------------------------------------------------
% PFCM
for i = 1:1100
  try
    [V1,V2,U1,U2, E1,E2] = Interval_T2_PFCM_clustering (Xin,nC,m1,m2,eta1,eta2,Cf,Cp,w);

    break;  % Break out of the i-loop on success
  catch ME
    disp(ME);
    fprintf('Retrying...\n');
  end
end
    
 



%% recunstruct
nC=ncolors;
Xin=I7;
for i=1:size(U1,2)
TEMP=0;
 for j=1:nC
     
  Lower_X=U1(j,i)*V1(j,:);
  Upper_X=U2(j,i)*V2(j,:);
    TEMP=TEMP+(Lower_X+Upper_X)/2;
 end
  Reconstruct_data(i,:)=TEMP/((sum(U1(:,i))+sum(U2(:,i)))/2) ;  
    
end

ERROR=0;
for i=1:size(U1,2)
    
A1=Xin(i,:)-Reconstruct_data(i,:);
A2=norm(A1)^2;
ERROR=ERROR+A2;
end
VR=sqrt(ERROR)/size(U1,2)






R1=Reconstruct_data(:,1);
R2=Reconstruct_data(:,2);
R3=Reconstruct_data(:,3);

Reverse_reshape_1=reshape(R1,[s(1) s(2) 1]);
Reverse_reshape_2=reshape(R2,[s(1) s(2) 1]);
Reverse_reshape_3=reshape(R3,[s(1) s(2) 1]);

Recunstructed_image=cat(3,Reverse_reshape_1,Reverse_reshape_2,Reverse_reshape_3);


subplot(121)
imshow(I1)
subplot(122)
imshow(Recunstructed_image)



