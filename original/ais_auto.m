
function [out] = ais_auto(cell,threshold,smooth,chosen)
close all

pixconv=0.161; %%%% microns per pixel
f = 0.33; %%%% is fraction of max fluo intensity at which AIS start+end parameters taken (0.33 as default; Grubb & Burrone 2010)

%%%% opening processed tif exported from ImageJ

disp(cell)              % Prints the file name
file = [cell '.tif - Processed method 2.5.tif'];
pic = imread(file);

    n=1;
    Ch{n} = pic(:,:,n);
    figure(n)
    imagesc(Ch{n})
    colormap(gray)
    axis square
    hold on

 %screen_size = get(0, 'ScreenSize'); %Use these two lines to make AIS %image full screen
 %set(figure(n), 'Position', [0 0 screen_size(3) screen_size(4) ] );

    hold on

%pause on
%pause

%%


Chais=mat2gray(Ch{n});                                        %Converts the ais channel to a 0-1 grayscale image.

% Chais_smooth=edge(Chais,smooth);

 Chais_smooth=imfilter(Chais,fspecial(smooth,[20 20],2));    %If smooth="gaussian", 2D gaussian smoothing (20x20, sd=2) is performed.

figure(n+1)
imagesc(Chais_smooth)                                           %Plots the smoothed image.
colormap(gray)
axis square


%%

if threshold==0
    threshold=graythresh(Chais_smooth);                         %If the input threshold=0, the graythresh function is used to find the optimal threshold (Otsu method)
end

%%
but=0;
loop=0;

 while but==0                                                   %This while-loop makes it possible to rethreshold the image when the max pixel is not in the selected region.
    ais_Ch=im2bw(Chais_smooth,threshold);                       %Thresholds the smoothed image.


    figure(n+2)
    imagesc(ais_Ch)                                             %Plots the thresholded image.
    colormap(gray)
    axis square
    title('threshold')

%     pause off

    %pause
%%

    SEo=strel('square',3);                                      %Structure elements are defined, for use with imopen en imclose.
    SEc=strel('square',3);

    ais_opened=imopen(ais_Ch,SEo);                              %Performs morphological opening on the thresholded image. Removes small objects.
    ais_closed=imclose(ais_opened,SEc);                         %Performs morphological closing on the opened image. Closes holes and gaps.

    figure(n+3)
    imagesc(ais_closed)                                         %Plots the opened/closed image.
    colormap(gray)
    axis square

    %pause

%%

    ais_CC=bwconncomp(ais_closed);                              %Finds all 8-connected components in the image.
    numPixels = cellfun(@numel,ais_CC.PixelIdxList);            %Finds the amount of pixels in each connected component.
    idx=find(numPixels<70);                                     %Finds which connected components are smaller than 70 pixels.
    ais_CC.PixelIdxList(idx)=[];                                %Removes those small connected components.
    ais_CC.NumObjects=length(ais_CC.PixelIdxList);              %Corrects the NumObjects parameter of ais_CC.


    %%

    figure(n)
    disp(' ')
    disp('zoom, then click near AIS start point')
    hold on
    zoom ON;pause;zoom OFF;
    [Xcell,Ycell]=ginput(1);                                    %Script requests single mouse click input for determination of which CC to use. Click just outside the start of the ais.
    Xcell=round(Xcell);                                         %ginput gives position on x-axis and y-axis in decimals, so needs to be rounded to find pixel coordinate.
    Ycell=round(Ycell);


    %%
    L=labelmatrix(ais_CC);                                      %Creates a matrix where each pixel of a CC has a identical integer, but different CC have different integers(labels).

    I=L(Ycell,Xcell);                                           %X and Y are in reverse order because y-axis corresponds with second pixel coordinate and vice-versa.
    if I==0                                                     %If the single mouse click was on a background pixel, we want to find the nearest nonzero pixel.
        L_bw=imbinarize (L,0.0000001);                                %Weird way of getting a binary image again.
        [dist,ind_dist]=bwdist(L_bw);                           %Creates distance matrix, and a matrix where each pixel gives the index of the closest pixel in L_bw.
        closest_pixel=ind_dist(Ycell,Xcell);                    %Finds closes pixel
        [Iclosest,Jclosest]=ind2sub(size(ind_dist),closest_pixel);
        I=L(closest_pixel);                                     %I is now the CC-label of the pixel closest to your mouse click.
    end
    ais_select=(L==I);                                          %Produces a binary image, containing only pixels with the label I -> one connected component.

    figure(n+4)
    imagesc(ais_select)                                         %Plots selected ais figure.
    colormap(gray)
    axis square


%%
    if loop==0                                                  %Makes sure these steps are only performed during the first iteration.
        D=im2double(Ch{n});
        max_all=max(max(D));                              %Calculates the maximum value in the original, unsmoothed and unthresholded image.
        max_ais=max(max(ais_select.*D));                  %Calculates the maximum value in the original image, but only in the region that is considered to be the ais.

            if max_all==max_ais
            but=1;                                              %If these two maxima are equal, rethresholding is not necessary and the script can continue.
            else
            disp('max pixel is not in ais region. Image is thresholded again with max from ais region')
            threshold=(max_ais/max_all)*threshold;              %Resets the input threshold, proportional to the max/max ratio.
            loop=1;                                             %While-loop will enter second iteration. Image will be rethresholded.
            end
    else
        but=1;
    end
end
%%

ais_skeleton=bwmorph(ais_select,'thin',Inf);                    %Morphological thinning procedure is performed, resulting in an 8-neighbour, single pixel wide 'skeleton'.
% regionprops(ais_skeleton,'Area','Perimeter')

figure(n+5)
imagesc(ais_skeleton)                                           %Plots the skeleton figure.
colormap(gray)
axis square

%%

figure(n+6)
imagesc(Ch{n})                                                %Plots the original figure again, so that we can plot other stuff on top of it.
colormap(gray)
axis square


hold on

[row,column]=ind2sub(size(ais_skeleton),find(ais_skeleton));    %On top of the raw figure, we are going to plot which pixels were selected in ais_skeleton.

plot(column,row,'.r')

%%

[dist_skel,ind_dist_skel]=bwdist(ais_skeleton);
closest_pixel_skel=ind_dist_skel(Ycell,Xcell);                  %again finds closest pixel to initial single mouse click.
[I_skel,J_skel]=ind2sub(size(ind_dist_skel),closest_pixel_skel);

Xais=I_skel;                                                    %These two variables will contain all coordinates of ais skeleton pixels, in the correct order.
Yais=J_skel;

ais_trace=ais_skeleton;                                         %Initializing for overly complicated double while-loop.
ais_traceq=ais_skeleton;
ais_trace(Xais(end),Yais(end))=0;
ais_traceq(Xais(end),Yais(end))=0;

loop2=0;
but2=1;
                                                                %The
                                                                %double while-loop will place the pixels in the correct order (beginning-end ais) and makes sure there are no branches.

while but2==1;                                                  %This loop makes sure the final result will have no branches, and will produce the longest route as a result.
    loop2=loop2+1;
                                                                %This loop will trace the ais skeleton by finding the neighbour of each pixel and then deleting the previous pixel.
    while ~isempty(find(ais_trace,1))                           %As long as the matrix ais_trace is not empty, the tracing is not done.
        Xend=Xais(end);                                         %Finds the pixel coordinates of the last found pixel.
        Yend=Yais(end);
        temp=zeros(size(ais_trace));                            %Creates a temporary matrix of the same size as ais_trace.
        temp(Xend-1:Xend+1,Yend-1:Yend+1)=ais_trace(Xend-1:Xend+1,Yend-1:Yend+1);   %temp now contains all zeros, apart from the direct neighbours of the center pixel.
        n=find(temp);                                           %Finds the index of neighbouring pixels in temp. Temp is same size as ais_trace, so indices wil correspond.
            if length(Xais)>300                                 %This is just to prevent repeated loops that take a long time.
                ais_trace=zeros(size(ais_trace));               %If the ais that is found is longer than 300 pixels, the tracing is terminated.
            elseif numel(n)>1                                   %If the center pixel has more than one neighbour, it means there are multiple branches.
                q=n(1);                                         %Just use one of the neighbouring pixel indices (takes first element of n and names it q).
                [Xnextq,Ynextq]=ind2sub(size(ais_trace),q);     %The next pixel coordinates are determined.
                Xais(end+1)=Xnextq;                             %Writes away pixel coordinates of next pixel in Xais and Yais.
                Yais(end+1)=Ynextq;
                ais_trace(Xnextq,Ynextq)=0;                     %Deletes this pixel from ais_trace, so it won't be found again in next loop.
            elseif numel(n)==1                                  %If the center pixel has only one neighbour, this pixel will be the next pixel.
                [Xnext,Ynext]=ind2sub(size(ais_trace),n);       %The next pixel coordinates are determined.
                Xais(end+1)=Xnext;                              %Writes away pixel coordinates of next pixel in Xais and Yais.
                Yais(end+1)=Ynext;
                ais_trace(Xnext,Ynext)=0;                       %Deletes this pixel from ais_trace, so it won't be found again in next loop.
            elseif numel(n)==0;                                 %If the center pixel has no neighbours, the while-loop will be terminated.
                ais_trace=zeros(size(ais_trace));               %Empties ais_trace, to make sure the while-loop will terminate itself.
            else
                error('something is wrong')
            end

        clear temp n
    end

    XaisL{loop2}=Xais;                                          %Writes away Xais ans Yais in XaisL and YaisL with the loop number as a subscript, in case multiple loops are necessary (branching).
    YaisL{loop2}=Yais;
    but2=0;

    if exist('q','var')                                         %If q exists, it means a branching was encountered. If multiple branches were encountered, only the last was saved in q.
       XYq=(Xais==Xnextq).*(Yais==Ynextq);                      %Xnextq and Ynextq were the coordinates of the neighbour of the branching point.
       qdel=find(XYq);                                          %Finds after which index in Xais and Yais, the saved coordinates belong to the branch.
       Xdel=Xais(qdel:end);                                     %Finds X-coordinates of all the pixels in the branch.
       Ydel=Yais(qdel:end);                                     %Finds Y-coordinates of all the pixels in the branch.
       ais_traceq(Xdel,Ydel)=0;                                 %All the pixels in the branch are now removed from ais_traceq.
       ais_trace=ais_traceq;                                    %All the pixels in the branch are now removed from ais_trace. Ais_traceq can be used in next loop.
       clear Xais Yais XYq qdel Xdel Ydel q
       Xais=I_skel;                                             %Re-initializing for the next loop.
       Yais=J_skel;
       ais_trace(Xais(end),Yais(end))=0;
       but2=1;                                                  %Makes sure the while-loop is not terminated yet. While-loop will terminate when no branches were found.
    end
    clear q

    if loop2>10                                                 %This is just to prevent repeated loops that take a long time.
        but2=0;
    end
end

%%

for l=1:loop2;
    long(l)=length(XaisL{l});                                   %Calculates the amount of pixels that are in each possible path, stored in the different elements of XaisL (and YaisL).
end

longline=find(long==max(long));                                 %Finds which of the possible paths is the longest.
Xais=XaisL{longline};                                           %The longest path is selected as the ais.
Yais=YaisL{longline};

%%
n=1
figure(n+7)
imagesc(Ch{n})                                                %Plots the original figure again, so that we can plot other stuff on top of it.
colormap(gray)
axis square

hold on

%%

%This part of the script is copied from ais_z3_pharma_lsm and is nearly
%identical. The main difference is that the csaps function is used to fit a
%2D spline through a 2D point cloud. The parameter now set to 0.3 can be
%used to decide how closely the spline will follow the points. At the end
%of the script, the length of the AIS is calculated in microns.


xy=[Yais;Xais];
xy=double(xy);

%%%% interpolate points with a spline curve and finer spacing.
t = 1:length(xy);
ts = [[1: 0.1: length(xy)],t(end)];    %%%% fine sampling, with double deletion below, ensures pixel-by-pixel line section
xysm = csaps(t,xy,0.3,ts);
for u = 2:(length(xysm))
    ax(u) = sqrt ( (xysm(1,u) - xysm(1,u-1))^2 + (xysm(2,u)-xysm(2,u-1))^2 );
end
ax = cumsum(ax);    %%%% so distance (in pixels) along spline from start of axon
ax_um = ax .* pixconv;   %%% distance along spline in um
xys = round(xysm); %%% rounding spline to whole-number pixel co-ordinates only
plot(xys(1,:),xys(2,:),'-b')
plot(xysm(1,:),xysm(2,:),'-y')

%%


[~,~,S]=ginput(1);                                              %ginput will show crosshairs, but that data is not used. the only thing that is saved in S is the button you press.
S=char(S);



if strcmp(S,'n')                                                %If you press the letter N on the keyboard, the output is 0. You basically discard this ais based on observation by eye. Handy for batch processing.
    out=0;
else
    double_id = [];
    for i = 1:length(xys(1,:))
        for j = (i+1):length(xys(1,:))
            if xys(1,i) == xys(1,j)
                if xys(2,i) == xys(2,j)
                    double_id = [double_id; j]; %%%% so finding double co-ordinates
                end
            end
        end
    end
    double_id = unique(double_id);
    alln = 1:length(xys(1,:));
    m = setdiff(alln,double_id);    %%% so m is index of all unique xy points in line section
    xs = xys(1,:); ys = xys(2,:);
    x_pix = xs(m);  %%%% so x_pix is unique array of axon x co-ordinates
    y_pix = ys(m);  %%%% y_pix is unique array of axon y co-ordinates

    for g = 1:length(x_pix) %%%% for each pixel, finding nearest location on spline axon
        d{g} = [];
        for h = 1:length(xysm)
            d{g} = [d{g}; sqrt((x_pix(g)-xysm(1,h))^2+(y_pix(g)-xysm(2,h))^2)];
        end
        [mind(g),mindi(g)] = min(d{g});
        x_ax(g) = xysm(1,mindi(g));
        y_ax(g) = xysm(2,mindi(g));
    end
    saxon_um = ax_um(mindi);   %%% so is array of distances of each pixel along axon - more accurate version of axon_um


end

savexy = [cell '_xy.txt'];
fid = fopen(savexy,'wt');
fprintf(fid,'%f\n',length(x_pix)); %%%% start of txt file gives number of points in arrays
for n = 1:length(x_pix)
    fprintf(fid,'%f\t %f\n',x_pix(n),y_pix(n)); %%%% writing x & y co-ordinates
end


savexy = [cell '_xy.txt'];


disp(cell)
close all

%%%%  Setting optionsd


nCh = 1;    %%% sets number of colour channels (max = 2)

byeye = 0;  %%% Choose '1' for by-eye measures, '0' for not (e.g. when running with 'chosen')


draw = 1;   %%% picks channel for drawing profile along axon: '1' for Channel1,'2' for Channel2
prof = 1;   %%% picks channel for profiling fluo intensity: '1' for Channel1,'2' for Channel2, or '3' for both
    if nCh == 1 & sum([draw prof])>2
        disp(' ')
        disp('nCh and draw/prof options are incompatible')
        stop
    end


%%%% finding appropriate folder to load axon co-ordinates, if using 'chosen'

q = length(cell);
while strcmp(cell(q),'/')==0
    q = q-1;
end
folder = cell(q-1);


%%%% opening Ch1 zprojection

Ch1_file = [cell '.tif'];
Ch1_pic = imread([Ch1_file]);
figure(1)
imagesc(Ch1_pic)
colormap(gray)
axis square
title('Ch1')
hold on

fid = fopen(savexy,'wt');
fprintf(fid,'%f\n',length(x_pix)); %%%% start of txt file gives number of points in arrays
for n = 1:length(x_pix)
    fprintf(fid,'%f\t %f\n',x_pix(n),y_pix(n)); %%%% writing x & y co-ordinates
end
fclose(fid);


%%%% obtaining fluorescence intensities for axonal profile - Ch1

for i = 1:length(x_pix) %%%% so working along the axon
    lv_c(i) = Ch1_pic(y_pix(i),x_pix(i));  %%% lv_c is fluorescence intensity at current point of axonal profile. Picture co-ordinates actually rows, then columns.
    lv_1(i) = Ch1_pic(y_pix(i)+1,x_pix(i));    %%%% this and next 7 rows get values for 3x3 roi centred on lv_c(i)
    lv_2(i) = Ch1_pic(y_pix(i)-1,x_pix(i));
    lv_3(i) = Ch1_pic(y_pix(i),x_pix(i)+1);
    lv_4(i) = Ch1_pic(y_pix(i),x_pix(i)-1);
    lv_5(i) = Ch1_pic(y_pix(i)+1,x_pix(i)+1);
    lv_6(i) = Ch1_pic(y_pix(i)-1,x_pix(i)-1);
    lv_7(i) = Ch1_pic(y_pix(i)+1,x_pix(i)-1);
    lv_8(i) = Ch1_pic(y_pix(i)-1,x_pix(i)+1);
    lv_smooth(i) = mean([lv_c(i) lv_1(i) lv_2(i) lv_3(i) lv_4(i) lv_5(i) lv_6(i) lv_7(i) lv_8(i)]); %%% lv_smooth is mean fluorescence intensity over 3x3 roi
end

v = num2str(lv_c); lv_c = str2num(v);   %%% re-formatting lv_c to avoid bugs

    %%% sliding mean to smooth intensity profile
for i = 1:length(x_pix)
    width = (1.5/(pixconv)); %%%sets no of pixels each side, i.e. for d = 20, width of sliding window is 41 (want width around 3um)
    d=((round(width))+1)
    if i<(d+1)  %%% at very start of axon profile , not allowing full window
        lv_slide(i) = mean([lv_c(1:i) lv_c(i:i+d)]);    %%% sliding mean with raw lv_c values
        lv_ss(i) = mean([lv_smooth(1:i) lv_smooth(i:i+d)]); %% sliding mean with lv_smooth values
    elseif i>(length(x_pix)-(d+1))  %%% at very end of axon profile, not allowing full window
        lv_slide(i) = mean([lv_c(i-d:i) lv_c(i:length(x_pix))]);
        lv_ss(i) = mean([lv_smooth(i-d:i) lv_smooth(i:length(x_pix))]);
    else   %%% in middle of axon profile, allowing full window
        lv_slide(i) = mean([lv_c(i-d:i) lv_c(i:i+d)]);
        lv_ss(i) = mean([lv_smooth(i-d:i) lv_smooth(i:i+d)]);
    end

end

axon_um = [1:length(x_pix)]*pixconv;    %%%% full length of axonal profile
xstep = axon_um(2)-axon_um(1);  %%%% sampling frequency along axonal profile, in um (should be equal to pixconv)

norm_lvslide = (lv_slide - min(lv_slide)) ./ (max(lv_slide)-min(lv_slide)); %%% is normalised, smoothed profile for 1x1 pixel sampling
norm_lv = (lv_ss - min(lv_ss)) ./ (max(lv_ss)-min(lv_ss));  %%%% is normalised, smoothed profile for 3x3 pixel sampling

figure(3)
subplot(2,2,1)
imagesc(Ch1_pic)
colormap(gray)
axis square
title('Ch1','color','g')
hold on
plot(x_pix,y_pix,'-b')
hold off
pix_narray = 1:length(x_pix); %%%% useful index array up to end of axonal profile

%%%%% measures of AIS location & length - Ch1

max_i = (find(norm_lv==1)); max_i = max_i(1);  %%% all measures use normalised, smoothed 3x3 pixel sampling profiles
max_x = pix_narray(max_i);  %%%% point index along axon where fluorescence intensity is highest
ais_end = find( (pix_narray>max_i) & (norm_lv>f),1,'last');
if length(ais_end)>0
    ais_end = ais_end(1);
else
    ais_end = x_pix(length(x_pix));  %%%% point index along axon past max where fluorescence intensity falls to f of its peak
end
ais_start = find( (pix_narray<max_i) & (norm_lv<f));
if length(ais_start)>0
    ais_start = ais_start(length(ais_start));
else
    ais_start = 0;%%%% point index along axon pre max where fluorescence intensity falls to f of its peak
end

debut = ais_start*pixconv;  %%%% AIS start position in Ch1, in um
fin = ais_end*pixconv;  %%%% AIS end position in Ch1, in um
lngth = fin-debut;      %%%% AIS length in Ch1, in um
mid = mean([debut fin]); %%%% AIS mid position in Ch1, in um
maxi = max_x*pixconv;   %%%% AIS max position in Ch1, in um

%%%% measures of AIS location & length - Ch2

%%%%% plotting

figure(3)
subplot(2,2,3)
plot(axon_um,lv_c,'g-') %%%% plotting raw Ch1 fluorescence intensity values in green
axis square
title('Raw')
if nCh>1 & prof>1
    hold on
    plot(axon_um,lv_cb,'r-') %%%% plotting raw Ch2 fluorescence intensity values in red
    hold off
end
subplot(2,2,4)
plot(axon_um,norm_lv,'g-') %%%% plotting normalised, smoothed Ch1 fluo in green
yline(0.33,'-.','Threshold')
hold on
if nCh>1 & prof>1
    plot(axon_um,norm_lvb,'r-') %%%% plotting normalised, smoothed Ch1 fluo in red
end
if prof==2
    plot([ais_end2*pixconv ais_end2*pixconv],[0 1],'b-')
    plot([ais_start2*pixconv ais_start2*pixconv],[0 1],'b-')
    plot([max_x2*pixconv max_x2*pixconv],[0 1],'b-')
    text((max(axon_um)-25),0.9,'Ch2 prof','color','b')
else %%% so if prof==1, or prof==3
    plot([ais_end*pixconv ais_end*pixconv],[0 1],'b-')
    plot([ais_start*pixconv ais_start*pixconv],[0 1],'b-')
    plot([max_x*pixconv max_x*pixconv],[0 1],'b-')  %%%% so plotting AIS start, max and end positions
    text((max(axon_um)-25),0.9,'Ch1 prof','color','b')
end

axis square
title('Smoothed & normalised')

hold off

%%%% results output

disp(' ')
disp('AIS Start'), display([debut]')
disp('AIS End'), display([fin]')
disp('AIS Mid'), display([mid]')
disp('AIS Max'), display([maxi]')
disp('AIS Length'), display([lngth]')
disp(' ')

clipboard('copy',num2str(lngth,6));

end




