% Problem 1a
% 1.Read the images
% 2. Calculate mean of each images read.
% 3. Subtract the mean from each image.

% filename=('/Users/sonalisreedhar/Documents/DIP/HW4/Prob_1a/train')
clear all;
clc;
Imagedata1=readraw1('brick1.raw',128,128);
Imagedata2=readraw1('brick2.raw',128,128);
Imagedata3=readraw1('brick3.raw',128,128);
Imagedata4=readraw1('brick4.raw',128,128);
Imagedata5=readraw1('brick5.raw',128,128);
Imagedata6=readraw1('brick6.raw',128,128);
Imagedata7=readraw1('brick7.raw',128,128);
Imagedata8=readraw1('brick8.raw',128,128);
Imagedata9=readraw1('brick9.raw',128,128);
Imagedata10=readraw1('blanket1.raw',128,128);
Imagedata11=readraw1('blanket2.raw',128,128);
Imagedata12=readraw1('blanket3.raw',128,128);
Imagedata13=readraw1('blanket4.raw',128,128);
Imagedata14=readraw1('blanket5.raw',128,128);
Imagedata15=readraw1('blanket6.raw',128,128);
Imagedata16=readraw1('blanket7.raw',128,128);
Imagedata17=readraw1('blanket8.raw',128,128);
Imagedata18=readraw1('blanket9.raw',128,128);
Imagedata19=readraw1('grass1.raw',128,128);
Imagedata20=readraw1('grass2.raw',128,128);
Imagedata21=readraw1('grass3.raw',128,128);
Imagedata22=readraw1('grass4.raw',128,128);
Imagedata23=readraw1('grass5.raw',128,128);
Imagedata24=readraw1('grass6.raw',128,128);
Imagedata25=readraw1('grass7.raw',128,128);
Imagedata26=readraw1('grass8.raw',128,128);
Imagedata27=readraw1('grass9.raw',128,128);
Imagedata28=readraw1('rice1.raw',128,128);
Imagedata29=readraw1('rice2.raw',128,128);
Imagedata30=readraw1('rice3.raw',128,128);
Imagedata31=readraw1('rice4.raw',128,128);
Imagedata32=readraw1('rice5.raw',128,128);
Imagedata33=readraw1('rice6.raw',128,128);
Imagedata34=readraw1('rice7.raw',128,128);
Imagedata35=readraw1('rice8.raw',128,128);
Imagedata36=readraw1('rice9.raw',128,128);

    
ImagedataNew1=getMeanRemoved(Imagedata1);
Imagedata1NeigborWin=getTextureFeature(ImagedataNew1);
featureVect(1,:)=getTextureFilteredImg(Imagedata1NeigborWin);

ImagedataNew2=getMeanRemoved(Imagedata2);
Imagedata2NeigborWin=getTextureFeature(ImagedataNew2);
featureVect(2,:)=getTextureFilteredImg(Imagedata2NeigborWin);

ImagedataNew3=getMeanRemoved(Imagedata3);
Imagedata3NeigborWin=getTextureFeature(ImagedataNew3);
featureVect(3,:)=getTextureFilteredImg(Imagedata3NeigborWin);

ImagedataNew4=getMeanRemoved(Imagedata4);
Imagedata4NeigborWin=getTextureFeature(ImagedataNew4);
featureVect(4,:)=getTextureFilteredImg(Imagedata4NeigborWin);

ImagedataNew5=getMeanRemoved(Imagedata5);
Imagedata5NeigborWin=getTextureFeature(ImagedataNew5);
featureVect(5,:)=getTextureFilteredImg(Imagedata5NeigborWin);

    
ImagedataNew6=getMeanRemoved(Imagedata6);
Imagedata6NeigborWin=getTextureFeature(ImagedataNew6);
featureVect(6,:)=getTextureFilteredImg(Imagedata6NeigborWin);


ImagedataNew7=getMeanRemoved(Imagedata7);
Imagedata7NeigborWin=getTextureFeature(ImagedataNew7);
featureVect(7,:)=getTextureFilteredImg(Imagedata7NeigborWin);

ImagedataNew8=getMeanRemoved(Imagedata8);
Imagedata8NeigborWin=getTextureFeature(ImagedataNew8);
featureVect(8,:)=getTextureFilteredImg(Imagedata8NeigborWin);

ImagedataNew9=getMeanRemoved(Imagedata9);
Imagedata9NeigborWin=getTextureFeature(ImagedataNew9);
featureVect(9,:)=getTextureFilteredImg(Imagedata9NeigborWin);

    
ImagedataNew10=getMeanRemoved(Imagedata10);
Imagedata10NeigborWin=getTextureFeature(ImagedataNew10);
featureVect(10,:)=getTextureFilteredImg(Imagedata10NeigborWin);

ImagedataNew11=getMeanRemoved(Imagedata11);
Imagedata11NeigborWin=getTextureFeature(ImagedataNew11);
featureVect(11,:)=getTextureFilteredImg(Imagedata11NeigborWin);

ImagedataNew12=getMeanRemoved(Imagedata12);
Imagedata12NeigborWin=getTextureFeature(ImagedataNew12);
featureVect(12,:)=getTextureFilteredImg(Imagedata12NeigborWin);

ImagedataNew13=getMeanRemoved(Imagedata13);
Imagedata13NeigborWin=getTextureFeature(ImagedataNew13);
featureVect(13,:)=getTextureFilteredImg(Imagedata13NeigborWin);

ImagedataNew14=getMeanRemoved(Imagedata14);
Imagedata14NeigborWin=getTextureFeature(ImagedataNew14);
featureVect(14,:)=getTextureFilteredImg(Imagedata14NeigborWin);

    
ImagedataNew15=getMeanRemoved(Imagedata15);
Imagedata15NeigborWin=getTextureFeature(ImagedataNew15);
featureVect(15,:)=getTextureFilteredImg(Imagedata15NeigborWin);

ImagedataNew16=getMeanRemoved(Imagedata16);
Imagedata16NeigborWin=getTextureFeature(ImagedataNew16);
featureVect(16,:)=getTextureFilteredImg(Imagedata16NeigborWin);

ImagedataNew17=getMeanRemoved(Imagedata17);
Imagedata17NeigborWin=getTextureFeature(ImagedataNew17);
featureVect(17,:)=getTextureFilteredImg(Imagedata17NeigborWin);

ImagedataNew18=getMeanRemoved(Imagedata18);
Imagedata18NeigborWin=getTextureFeature(ImagedataNew18);
featureVect(18,:)=getTextureFilteredImg(Imagedata18NeigborWin);
    
ImagedataNew19=getMeanRemoved(Imagedata19);
Imagedata19NeigborWin=getTextureFeature(ImagedataNew19);
featureVect(19,:)=getTextureFilteredImg(Imagedata19NeigborWin);

ImagedataNew20=getMeanRemoved(Imagedata20);
Imagedata20NeigborWin=getTextureFeature(ImagedataNew20);
featureVect(20,:)=getTextureFilteredImg(Imagedata20NeigborWin);

ImagedataNew21=getMeanRemoved(Imagedata21);
Imagedata21NeigborWin=getTextureFeature(ImagedataNew21);
featureVect(21,:)=getTextureFilteredImg(Imagedata21NeigborWin);

ImagedataNew22=getMeanRemoved(Imagedata22);
Imagedata22NeigborWin=getTextureFeature(ImagedataNew22);
featureVect(22,:)=getTextureFilteredImg(Imagedata22NeigborWin);

ImagedataNew23=getMeanRemoved(Imagedata23);
Imagedata23NeigborWin=getTextureFeature(ImagedataNew23);
featureVect(23,:)=getTextureFilteredImg(Imagedata23NeigborWin);


ImagedataNew24=getMeanRemoved(Imagedata24);
Imagedata24NeigborWin=getTextureFeature(ImagedataNew24);
featureVect(24,:)=getTextureFilteredImg(Imagedata24NeigborWin);

ImagedataNew25=getMeanRemoved(Imagedata25);
Imagedata25NeigborWin=getTextureFeature(ImagedataNew25);
featureVect(25,:)=getTextureFilteredImg(Imagedata25NeigborWin);

ImagedataNew26=getMeanRemoved(Imagedata26);
Imagedata26NeigborWin=getTextureFeature(ImagedataNew26);
featureVect(26,:)=getTextureFilteredImg(Imagedata26NeigborWin);

ImagedataNew27=getMeanRemoved(Imagedata27);
Imagedata27NeigborWin=getTextureFeature(ImagedataNew27);
featureVect(27,:)=getTextureFilteredImg(Imagedata27NeigborWin);

ImagedataNew28=getMeanRemoved(Imagedata28);
Imagedata28NeigborWin=getTextureFeature(ImagedataNew28);
featureVect(28,:)=getTextureFilteredImg(Imagedata28NeigborWin);
    
ImagedataNew29=getMeanRemoved(Imagedata29);
Imagedata29NeigborWin=getTextureFeature(ImagedataNew29);
featureVect(29,:)=getTextureFilteredImg(Imagedata29NeigborWin);

ImagedataNew30=getMeanRemoved(Imagedata30);
Imagedata30NeigborWin=getTextureFeature(ImagedataNew30);
featureVect(30,:)=getTextureFilteredImg(Imagedata30NeigborWin);

ImagedataNew31=getMeanRemoved(Imagedata31);
Imagedata31NeigborWin=getTextureFeature(ImagedataNew31);
featureVect(31,:)=getTextureFilteredImg(Imagedata31NeigborWin);

ImagedataNew32=getMeanRemoved(Imagedata32);
Imagedata32NeigborWin=getTextureFeature(ImagedataNew32);
featureVect(32,:)=getTextureFilteredImg(Imagedata32NeigborWin);

ImagedataNew33=getMeanRemoved(Imagedata33);
Imagedata33NeigborWin=getTextureFeature(ImagedataNew33);
featureVect(33,:)=getTextureFilteredImg(Imagedata33NeigborWin);


ImagedataNew34=getMeanRemoved(Imagedata34);
Imagedata34NeigborWin=getTextureFeature(ImagedataNew34);
featureVect(34,:)=getTextureFilteredImg(Imagedata34NeigborWin);

ImagedataNew35=getMeanRemoved(Imagedata35);
Imagedata35NeigborWin=getTextureFeature(ImagedataNew35);
featureVect(35,:)=getTextureFilteredImg(Imagedata35NeigborWin);

ImagedataNew36=getMeanRemoved(Imagedata36);
Imagedata36NeigborWin=getTextureFeature(ImagedataNew36);
featureVect(36,:)=getTextureFilteredImg(Imagedata36NeigborWin);

% for i= 1:25
%     figure()
%     hist(featureVect(i,:),30)
% end
for i=1:36
    feature_vector_std(i,:)=(featureVect(i,:)-mean(featureVect(i,:)))./std(featureVect(i,:),1);
end

function[E]=getAvgEnergy(texture_filtered)
texture_filtered=double(texture_filtered);
sum=0;
for i=1:128
    for j=1:128
        sum=sum+abs(texture_filtered(i,j));
    end
end
E=1/(128*128).*sum;
end

function[texture_filtered]= convFilter(neighbour_window,filter)
N=5;
for i=1+(N-1)/2:128+(N-1)/2
    for j=1+(N-1)/2:128+(N-1)/2
        texture_filtered(i-(N-1)/2,j-(N-1)/2)=sum(sum(neighbour_window.*filter));
    end
end
end
               
         
                
function [Energy]=getTextureFilteredImg(neighbour_window)
%to get all 25 features from an image
N=5;

%Filter declaration
L5=[1 4 6 4 1];
E5=[-1 -2 0 2 1];
S5=[-1 0 2 0 -1];
W5=[-1 2 0 -2 1];
R5=[1 -4 6 -4 1];

fil1=L5'*L5;
fil2=((L5'*E5)+(E5'*L5))/2;
fil3=((L5'*S5)+(S5'*L5))/2;
fil4=((L5'*W5)+(W5'*L5))/2;
fil5=((L5'*R5)+(R5'*L5))/2;


% E5L5=E5'*L5;
fil6=E5'*E5;
fil7=((E5'*S5)+(S5'*E5))/2;
fil8=((E5'*W5)+(W5'*E5))/2;
fil9=((E5'*R5)+(R5'*E5))/2;

% W5L5=W5'*L5;
% W5E5=W5'*E5;
fil10=((W5'*S5)+(W5'*E5))/2;
fil11=W5'*W5;
fil12=((W5'*R5)+(R5'*W5))/2;

% S5L5=S5'*L5;
% S5E5=S5'*E5;
fil13=S5'*S5;
% S5W5=S5'*W5;
fil14=((S5'*R5)+(R5'*S5))/2;

% R5L5=R5'*L5;
% R5E5=R5'*E5;
% R5S5=R5'*S5;
% R5W5=R5'*W5;
fil15=R5'*R5;
for i=1:5
    for j=1:5
        texture_filtered1=convFilter(neighbour_window,fil1);
        texture_filtered2=convFilter(neighbour_window,fil2);
        texture_filtered3=convFilter(neighbour_window,fil3);
        texture_filtered4=convFilter(neighbour_window,fil4);
        texture_filtered5=convFilter(neighbour_window,fil5);
        texture_filtered6=convFilter(neighbour_window,fil6);
        texture_filtered7=convFilter(neighbour_window,fil7);
        texture_filtered8=convFilter(neighbour_window,fil8);
        texture_filtered9=convFilter(neighbour_window,fil9);
        texture_filtered10=convFilter(neighbour_window,fil10);
        texture_filtered11=convFilter(neighbour_window,fil11);
        texture_filtered12=convFilter(neighbour_window,fil12);
        texture_filtered13=convFilter(neighbour_window,fil13);
        texture_filtered14=convFilter(neighbour_window,fil14);
        texture_filtered15=convFilter(neighbour_window,fil15);
        
    end
end
      

E1=getAvgEnergy(texture_filtered1);
E2=getAvgEnergy(texture_filtered2);
E3=getAvgEnergy(texture_filtered3);
E4=getAvgEnergy(texture_filtered4);
E5=getAvgEnergy(texture_filtered5);

E6=getAvgEnergy(texture_filtered6);
E7=getAvgEnergy(texture_filtered7);
E8=getAvgEnergy(texture_filtered8);
E9=getAvgEnergy(texture_filtered9);
E10=getAvgEnergy(texture_filtered10);

E11=getAvgEnergy(texture_filtered11);
E12=getAvgEnergy(texture_filtered12);
E13=getAvgEnergy(texture_filtered13);
E14=getAvgEnergy(texture_filtered14);
E15=getAvgEnergy(texture_filtered15);


Energy=[E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15];
end
   




function[texture_extended]= getMeanRemoved(texture)
texture_sum=0;
% texture=readraw1('brick1.raw',128,128);
% figure(1)
% imshow(texture)
texture1=double(texture);
N=5;
for i=1:128
    for j=1:128
%         
        texture_sum=texture_sum+(texture1(i,j));
%        
    end
end
texture_mean=(1/((128)*(128))).*texture_sum;
for i=1:128
    for j=1:128

        texture_mean_removed(i,j)=(texture(i,j)-texture_mean);
    end
end
% figure(2)
% imshow(texture_mean_removed)
% 
texture_extended=zeros(128+(N-1),128+(N-1));

for i=1:128          %data centered
    for j=1:128
        texture_extended(i+(N-1)/2,j+(N-1)/2)=texture_mean_removed(i,j);
    end
end

%boundary extension
for i=0:((N-1)/2)-1
    texture_extended(:,((N-1)/2)-i)=texture_extended(:,((N-1)/2)+i+2);
    texture_extended(:,128+((N-1)/2)+i+1)=texture_extended(:,128+((N-1)/2)-i-1);
    texture_extended(((N-1)/2)-i,:)=texture_extended(((N-1)/2)+i+2,:);
    texture_extended((128+(N-1)/2)+i+1,:)=texture_extended((128+(N-1)/2)-i-1,:);
end
end



%to get convolution of image and 25 5x5 filters
function[neighbour_window]=getTextureFeature(texture_mean_removed)
N=5;
for i=1+(N-1)/2:128+(N-1)/2
    for j=1+(N-1)/2:128+(N-1)/2
        neighbour_window=texture_mean_removed((i-(N-1)/2:i+(N-1)/2),(j-(N-1)/2:j+(N-1)/2));
    end
end
end

