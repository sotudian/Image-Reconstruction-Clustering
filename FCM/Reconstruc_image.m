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

ncolors=5;


[V,U]=fcm(I7,ncolors);

U1=U;
U2=U;
V1=V;
V2=V;
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








% 
% i1=1;
% i2=1;
% 
% k1=[];
% k2=[];
% 
% 
% 
% 
% 
% for i=1:(s(1)*s(2))
%   a=max(U(:,i))  ;
%    b=find(U(:,i)==a) ;
%    if b==1
%        k1(i,:)=I7(i,:);
%       k2(i,:)=[0.5,0.5,0.5];
%    else
%        k2(i,:)=I7(i,:);
%       k1(i,:)=[0.5,0.5,0.5];
%    
%    end   
%    
% end
% k11=1-k1;
% k22=1-k2;
% 
% cluster1R=zeros(s(1),s(2));
% cluster1G=zeros(s(1),s(2));
% cluster1B=zeros(s(1),s(2));
% cluster2R=zeros(s(1),s(2));
% cluster2G=zeros(s(1),s(2));
% cluster2B=zeros(s(1),s(2));
% for i=1:(s(1)*s(2))
%  cluster1R(i)=k11(i,1);  
%  cluster1G(i)=k11(i,2); 
%  cluster1B(i)=k11(i,3); 
%  
%  cluster2R(i)=k22(i,1);  
%  cluster2G(i)=k22(i,2); 
%  cluster2B(i)=k22(i,3); 
%     
% end
% C1=cat(3,cluster1R,cluster1G,cluster1B);
% C2=cat(3,cluster2R,cluster2G,cluster2B);
