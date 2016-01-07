%%%%%%%%%%%%%%%%%%%%%
%barcode located use connected component
%recognition use SVM  
%sharpen image
% 2015/1/7 edited 
%%%%%%%%%%%%%%%%%%%%%
clc;
clf;
clear all;
close all;
%% set parameters
gthres=1.5;
neggthres=-gthres;
Toprate=1;
numRegionHeighTop=0.3;
numRegionHeighBottom=5;
numRegionwidthTop=0.3;
numRegionwidthBottom=3;
whRate=0.2;
yPosition=0.5;
hdiff=4;
pjWHrate=0.3;
yCentroid=0.35;
pjWidth=3;

side1=16;
side2=16;

%% input image
I0=imread('barcode19.jpg');
I0=imresize(I0,[488 648]);
figure;
imshow(I0);
% I0 = imsharpen(I0,'Radius',5,'Amount',2);
% figure, imshow(I0), title('Sharpened Image');
% I0=adaHSV_Saturation(I0);%% 呼叫副程式做飽和度增強
% figure;
% imshow(I0);
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
for rgnProps=1:length(stats);
    if (stats(rgnProps).BoundingBox(3)/stats(rgnProps).BoundingBox(4)<3) && (stats(rgnProps).BoundingBox(3)/stats(rgnProps).BoundingBox(4)>Toprate)
         subimage1=imcrop(I0,stats(rgnProps).BoundingBox);
         subimage1 = imsharpen(subimage1,'Radius',3,'Amount',2.5);   %%% 3    2
         figure;
         imshow(subimage1);
         level = graythresh(subimage1);  %%ostu
         bw=im2bw(subimage1,level);
%         bw=imprBernsen(subimage1);

        [rotI,Theta]=barcode_rotate4(bw);  %%%%%%%%%rotateimage
        rotI=1-rotI;           
        figure;
        imshow(rotI);
        [rota,rotb]=size(rotI);
  %% find number region
         [subLB,NUM2]=bwlabel(rotI,8); %% connected component 
         substats=regionprops(subLB,'Basic');
         NBoundingBox=zeros(length(substats),4);
         count=0;
         for i=1:length(substats)
             if (substats(i).BoundingBox(4)<rota*numRegionHeighTop) && (substats(i).BoundingBox(4)>numRegionHeighBottom)...
                     &&(substats(i).BoundingBox(3)<stats(rgnProps).BoundingBox(3)*numRegionwidthTop)&&(substats(i).BoundingBox(3)>numRegionwidthBottom)...
                 && (substats(i).Centroid(2)>rota*yPosition) && (substats(i).BoundingBox(3)/substats(i).BoundingBox(4)>whRate) 
                   count=count+1
                   NBoundingBox(count,1:4)=substats(i).BoundingBox;
%                    numRegion=imcrop(rotI, NBoundingBox(count,1:4));      %%test for number count 
%                    figure;                                               %%test for number count 
%                    imshow(numRegion)                                     %%test for number count 
             end
         end %% end of number region statistics
         yStart=min( NBoundingBox(1:count,2)); 
         yEnd=max( NBoundingBox(1:count,2)); 
         hmax=max( NBoundingBox(1:count,4)); 
         hmedian=median(NBoundingBox(1:count,4));
         cropI=imcrop(rotI,[1 yStart rotb yEnd-yStart+hmax]);
         figure;
         imshow(cropI);
         
         diffcount=0;
         N2BoundingBox=zeros(count,4);
         for i=1:count
             if  abs(hmedian-NBoundingBox(i,4))>hdiff
                 count=count-1
             else    
                 diffcount=diffcount+1;
                 N2BoundingBox(diffcount,1:4)=NBoundingBox(i,1:4);
             end
         end
         
        if count <=9
            continue;
        end
     %%    segmentation
       pjcount=0;
         if count>=13
              for i=1:count
                  pjcount=pjcount+1
                  ccaNumber=imcrop(rotI,N2BoundingBox(i,1:4));
                  figure;
                  imshow(ccaNumber); 
                  lv=graythresh(double(ccaNumber));
                  segI=im2bw(ccaNumber,lv);
                  segI=double(imresize(segI,[side1 side2]));
                  vseg=reshape(segI,[1,side1*side2]);
                  TestSample(pjcount,:)=vseg;
              end
         else   
               
                t2 =vprojection(cropI);
                [cropa,cropb]=size(cropI);
                [initial]=criticalPoint(t2,0.05,cropa) ;
                start=1;End=2;
                for i=1:length(initial)
                    width=abs(initial(i,End)-initial(i,start));
                    pjNumber=imcrop(cropI,[ initial(i,start) 1 width cropa]);
                      if width>0.05*cropb %0.05
                       t3 =vprojection(pjNumber);   
                       [initial2]=criticalPoint2(t3,0.2,cropa) ;      
                                  for i=1:length(initial2)
                                         width2=abs(initial2(i,End)-initial2(i,start));
                                         divNumber=imcrop(pjNumber,[ initial2(i,start) 1 width2 cropa]);
                                         
                                         [segIm,segIn]=size(divNumber);pixel=sum(sum(divNumber));
                                         xi = ones(segIm,1)*[1:segIn]; yi = [1:segIm]'*ones(1,segIn);
                                         meany = sum(sum(divNumber.*yi))/pixel;
                        
                                         if pixel>segIm*segIn*pjWHrate && meany>segIm*yCentroid
                                         figure;
                                         imshow(divNumber);
                                         pjcount=pjcount+1
                                         lv=graythresh(double(divNumber));
                                         segI=im2bw(divNumber,lv);
                                         segI=double(imresize(segI,[side1 side2]));
                                         vseg=reshape(segI,[1,side1*side2]);
                                         TestSample(pjcount,:)=vseg;                                      
                                         end
                                  end    
                      else
                          
                            if width<=0;
                                continue;
                            end  
                              [segIm,segIn]=size(pjNumber);pixel=sum(sum(pjNumber));
                              xi = ones(segIm,1)*[1:segIn]; yi = [1:segIm]'*ones(1,segIn); 
                              meany = sum(sum(pjNumber.*yi))/pixel;
%                         
                            if  (pixel>segIm*segIn*pjWHrate) && (meany>segIm*yCentroid) && (segIn>pjWidth)
                                 figure;
                                 imshow(pjNumber);
                                 pjcount=pjcount+1
                                 lv=graythresh(double(pjNumber));
                                 segI=im2bw(pjNumber,lv);
                                 segI=double(imresize(segI,[side1 side2]));
                                 vseg=reshape(segI,[1,side1*side2]);
                                 TestSample(pjcount,:)=vseg;
                           end
                      end %% end of width threshold
                 end
         end %% end of 13 number count 
    
    
    end  
end


%% SVM result
testlabel=zeros(pjcount,1);
models=load('libtrainingModel5.mat','models');
models=models.models;
[predicted_label, accuracy, decision_values] = libsvmpredict(testlabel, TestSample, models);
for i=1:length(predicted_label)
    fprintf('%d',predicted_label(i));
    
    
end
