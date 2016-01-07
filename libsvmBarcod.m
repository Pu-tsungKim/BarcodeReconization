%%%%%%%%%%%%%%%%%%%%%
%barcode located use connected component
%recognition use SVM  
% 2015/12/24 edited 
%%%%%%%%%%%%%%%%%%%%%
clc;
clf;
clear all;
close all;
%% set parameters
gthres=1.5;
neggthres=-gthres;
erodeRate=0.04;
cropthres=0.4;
numlong=0.25;
Toprate=1;
numberheight=0.55;
side=16;
centroidPosition=20;
pjthres=0.2;
%% input image
I0=imread('barcode16.jpg');
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
SE = strel('rectangle',[1 20]); 
dia1 = imdilate(J,SE);
figure;
imshow(dia1);
dia1 =bwareaopen(dia1,8000);  %%小於8000不要
figure;
imshow(dia1);
%% regionprop
[LB,NUM]=bwlabel(dia1,8); %% connected component 
stats=regionprops(LB,'Basic');
for irgnProps=1:length(stats);
    if (stats(irgnProps).BoundingBox(3)/stats(irgnProps).BoundingBox(4)<3) && (stats(irgnProps).BoundingBox(3)/stats(irgnProps).BoundingBox(4)>Toprate)
         subimage1=imcrop(I0,stats(irgnProps).BoundingBox);
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
        for i=1:5:300
            k=i*0.01;
        sobSE=strel('rectangle',[round(k*stats(irgnProps).BoundingBox(4)) 2]);%%%for erode
        soberd =imerode(sobl,sobSE);
              if   sum(sum(soberd==1))<erodeRate*stats(irgnProps).Area
                     figure;
                     imshow(soberd);
                  break;
              end    
        end
        figure;
        imshow(soberd);
         %% horizontal projection
        t=zeros(1,255);
        [Accm,Accn]=size(soberd);
   for i=1:Accm %%for loop vertical projection
     t(i)=0;
     for j=1:Accn
     if soberd(i,j)>=1
         t(i)=t(i)+1;
     end
     end 
     end%%end for projection
  t=uint8(t);
   figure;
   bar(t);          
       if sum(sum(soberd==1)) >50 %% 被侵蝕後像素值
        sobSE=strel('rectangle',[2 round(0.0625*rotb)]);%%% for dialate
        sobdia=imdilate(soberd,sobSE);
         figure;
        imshow(sobdia);
        title('after mophorlogy');
      %%%projection for cropping number
        [m1,n1]=size(sobdia);
        rowEnd=0;
        for i=1:m1
           k=m1-i+1;
            prjsum=0;
            for j=1:n1
                if  sobdia(k,j)>=1
                    prjsum=prjsum+1;
                end
            end
            if prjsum/n1>=cropthres
                rowEnd=k;
                break;
            end
        end 
    
        code=imcrop(rotI,[1,rowEnd+0.03*stats(irgnProps).BoundingBox(4),n1,numlong*stats(irgnProps).BoundingBox(4)]);
           figure;
        imshow(code);
        AccI=Accurateplate(code); 
        figure;
        imshow(AccI);
        %% horizontal projection
        t=zeros(1,255);
        [Accm,Accn]=size(AccI);
   for i=1:Accm %%for loop vertical projection
     t(i)=0;
     for j=1:Accn
     if AccI(i,j)>0
         t(i)=t(i)+1;
     end
     end 
     end%%end for projection
  t=uint8(t);
   figure;
   bar(t);
   %% critical point
start=1;End=2;
label=start;index=1;
% critical(1,1)=1; critical(2,2)=1;
    for i=1:length(t)
        if label==start
           if t(i)>Accn*pjthres
             critical(index,label)=i;
             label=End;
           end

        end 
 
        if label==End
             if t(i)<Accn*pjthres
              critical(index,label)=i;
              label=start;
              index=index+1;
             end
        end
    end
    if critical(1,2)-critical(1,1)<Accm*0.4
   AccI2=imcrop(AccI,[1,critical(2,1),Accn,Accm-critical(2,1)]);
    else
  AccI2=AccI;
    end
        %% segmentation
        [m2,n2]=size(AccI2);
%         AccI2=bwareaopen(AccI2,round(0.01*m2*n2));  
        [LB2,NUM2]=bwlabel(AccI2,4); %%4 connected component 
         substats=regionprops(LB2,'Basic');
         numcount=0;
         for i=1:length(substats)
              if  (substats(i).BoundingBox(3)>5) && (substats(i).BoundingBox(4)>m2*numberheight)...
                      && (substats(i).Centroid(1)>n2/centroidPosition)  && (substats(i).Centroid(1)< (n2*(centroidPosition-1)/centroidPosition))
                  numimage=imcrop(AccI2,substats(i).BoundingBox);
                  [a,b]=size(numimage);
                  if  b>a             
                  numimage=imcrop(AccI2,[substats(i).BoundingBox(1),substats(i).BoundingBox(2)...
                      substats(i).BoundingBox(3)/2,substats(i).BoundingBox(4)]);
                  end  
                  numcount=numcount+1;
                  figure;
                  imshow(numimage);
                   %Rcgn=numrecognition(numimage);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Support Vector Mechine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 numimage=imresize(numimage,[side side]);
                 numimage =double( reshape(numimage,[1,side*side])); 
                 TestSample(numcount,:)=numimage;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
                  end               
         end
        
    end  %loop for sum(sum(soberd)) >50
   end   %loop for roi
end   %loop for status


%% SVM result
%            models=load('trainingModel.mat','models');
%            models=models.models;
  testlabel=zeros(numcount,1);
            models=load('libtrainingModel5.mat','models');
    models=models.models;
%     for j=1:size(TestSample,1)
%     for k=1:numclass
%         if(svmclassify(models(k),TestSample(j,:))) 
%             break;
%         end
%     end
%     result(j) = k-1;
% end
  [predicted_label, accuracy, decision_values] = libsvmpredict(testlabel, TestSample, models);
     for i=1:length(predicted_label)
        fprintf('%d',predicted_label(i));
        

     end
