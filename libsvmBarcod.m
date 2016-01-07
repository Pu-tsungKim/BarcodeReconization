%%%%%%%%%%%%%%%%%%%%%
%barcode located use connected component
%recognition use SVM  
% 2015/12/23 edited 
%%%%%%%%%%%%%%%%%%%%%
clc;
clf;
clear all;
close all;
%% set parameters
gthres=1.5;
neggthres=-gthres;
fstThres=0.3;
sndThres=0.3;
fstThres2=0.05;
sndThres2=0.05;
cropthres=0.15;
numlong=35;
Toprate=1;
numberheight=0.55;
side=16;
centroidPosition=20;
pjthres=0.2;
%% input image
I0=imread('barcode2.jpg');
I0=imresize(I0,[488 648]);
figure;
imshow(I0);
I0=adaHSV_Saturation(I0);%% 呼叫副程式做飽和度增強
figure;
imshow(I0);
%% gradient
I=double(I0);
G=gradient(I);
figure;
imshow(G);
Gmin=min(G(:));
Gmax=max(G(:));
N=(20*(G-Gmin)/(Gmax-Gmin))-10;
figure;
imshow(N);
[row,col,dim]=size(N);
sumth=0;
for i=1:row
     
    for j=1:col
        if N(i,j)> gthres
                N(i,j,1)=0;
                 N(i,j,2)=0;
                 N(i,j,3)=1; %blue is peak
        elseif N(i,j)<neggthres
                 N(i,j,1)=0;
                 N(i,j,2)=1; %green is valley
                 N(i,j,3)=0;
        else
            N(i,j,1)=0;
                 N(i,j,2)=0; 
                 N(i,j,3)=0;
        end
    end
end
figure;
imshow(N);
W=im2bw(N);
J = medfilt2(W,[3 3]);
figure;
imshow(J);
%% mophorlogy
SE = strel('rectangle',[2 20]); 
dia1 = imdilate(J,SE);
figure;
imshow(dia1);
dia1 =bwareaopen(dia1,8000);  %%小於8000不要
figure;
imshow(dia1);
%% regionprop
[LB,NUM]=bwlabel(dia1,8); %% connected component 
stats=regionprops(LB,'Basic');
for i=1:length(stats);
    if (stats(i).BoundingBox(3)/stats(i).BoundingBox(4)<3) && (stats(i).BoundingBox(3)/stats(i).BoundingBox(4)>Toprate)
        subimage1=imcrop(I0,stats(i).BoundingBox);
        level = graythresh(subimage1);  %%ostu
        bw=im2bw(subimage1,level);
        %         bw=imprBernsen(subimage1);
        
        [rotI,Theta]=barcode_rotate4(bw);  %%%%%%%%%rotateimage
        rotI=1-rotI;
        figure;
        imshow(rotI);
        [rota,rotb]=size(rotI);
        % bw=im2bw(rotI);
        y_mask= [-1 0 1;-2 0 2;-1 0 1];
        sobl = abs(filter2(y_mask,rotI));
        figure;
        imshow(sobl);
        SE = strel('rectangle',[1 10]); 
        dia2 = imdilate(sobl,SE);
        figure;
        imshow(dia2);
       t=vprojection(dia2);
       %% critical point
       start=1;End=2;
       label=start;index=1;
        critical(1,1)=1; critical(2,2)=1;
       for i=1:length(t)
           if label==start
               if t(i)>rota*fstThres
                   critical(index,label)=i;
                   label=End;
               end
               
           end
           
           if label==End
               if t(i)<rota*sndThres
                   critical(index,label)=i;
                   label=start;
                   index=index+1;
               end
           end
       end
       
      for i=1:length(critical)
          width=critical(i,2)-critical(i,1);
          if width<15
              continue;
          end    
         barregion=imcrop(rotI,[1,critical(i,1),rotb,width]);
          figure;
          imshow(barregion);
          t2=vprojection(barregion);
         %critical point2
        start=1;End=2;
       label=start;index=1;
%         critical2(1,1)=1; critical2(2,2)=1;
       for i=1:length(t2)
           if label==start
               if t2(i)>rota*fstThres2
                   critical2(index,label)=i;
                   label=End;
               end
               
           end
           
           if label==End
               if t2(i)<rota*sndThres2
                   critical2(index,label)=i;
                   label=start;
                   index=index+1;
               end
           end
       end

  
        
          minpj=find(t2==min(t2));
          Numregion=imcrop(barregion,[1,minpj(1),rotb,rota-minpj(1)]);
          figure;
          imshow(Numregion);
          
      end
      
      
%       Numregion=imcrop(rotI,[1,critical(2,1)-3,rotb,rota-critical(2,1)]);
%       figure;
%       imshow(Numregion);
end   %loop for roi
end   %loop for status


