function matlab_rethreshold_reference(basepath, threshold, Xcell, Ycell, outfile)
% MATLAB_RETHRESHOLD_REFERENCE  Ground truth for the *rethreshold loop* specifically.
%
% ``matlab_reference.m`` runs the original's ``while but==0`` loop single-pass, which is
% the fixed-threshold mode. This one runs the loop **as written**: threshold, click, and
% if the image's brightest pixel is not inside the clicked component, rescale by
% max_ais/max_all, segment the whole image again and take the second click.
%
% Both clicks use the supplied (Xcell, Ycell), which is what a person does in practice:
% the prompt is "click near AIS start point" both times, and they aim at the same AIS.
%
% Everything the loop touches is dumped so the Python port can be compared step by step,
% not just on the final number. Coordinates in the output are MATLAB 1-based.

pixconv = 0.161;
f = 0.33;

pic = imread([basepath '.tif - Processed method 2.5.tif']);
n = 1;
Ch{n} = pic(:,:,n);

Chais = mat2gray(Ch{n});
Chais_smooth = imfilter(Chais, fspecial('gaussian',[20 20],2));

if threshold == 0
    threshold = graythresh(Chais_smooth);
end

thresholds = [];       % the level used on each pass of the loop
labels_sel = [];
cc_counts  = [];
cc_pixels  = [];

but  = 0;
loop = 0;

while but == 0
    thresholds(end+1) = threshold;

    ais_Ch = im2bw(Chais_smooth, threshold);

    SEo = strel('square',3);
    SEc = strel('square',3);
    ais_opened = imopen(ais_Ch, SEo);
    ais_closed = imclose(ais_opened, SEc);

    ais_CC = bwconncomp(ais_closed);
    numPixels = cellfun(@numel, ais_CC.PixelIdxList);
    idx = find(numPixels < 70);
    ais_CC.PixelIdxList(idx) = [];
    ais_CC.NumObjects = length(ais_CC.PixelIdxList);

    L = labelmatrix(ais_CC);

    I = L(Ycell, Xcell);
    if I == 0
        L_bw = imbinarize(L, 0.0000001);
        [dist, ind_dist] = bwdist(L_bw);
        closest_pixel = ind_dist(Ycell, Xcell);
        I = L(closest_pixel);
    end
    ais_select = (L == I);

    cc_counts(end+1) = ais_CC.NumObjects;
    labels_sel(end+1) = double(I);
    cc_pixels(end+1) = sum(ais_select(:));

    if loop == 0
        D = im2double(Ch{n});
        max_all = max(max(D));
        max_ais = max(max(ais_select .* D));
        out_max_all = max_all;
        out_max_ais = max_ais;
        if max_all == max_ais
            but = 1;
        else
            threshold = (max_ais/max_all)*threshold;
            loop = 1;
        end
    else
        but = 1;
    end
end

ais_skeleton = bwmorph(ais_select,'thin',Inf);

[dist_skel, ind_dist_skel] = bwdist(ais_skeleton);
closest_pixel_skel = ind_dist_skel(Ycell, Xcell);
[I_skel, J_skel] = ind2sub(size(ind_dist_skel), closest_pixel_skel);

Xais = I_skel;
Yais = J_skel;

ais_trace  = ais_skeleton;
ais_traceq = ais_skeleton;
ais_trace(Xais(end), Yais(end))  = 0;
ais_traceq(Xais(end), Yais(end)) = 0;

loop2 = 0;
but2  = 1;

while but2 == 1
    loop2 = loop2 + 1;
    while ~isempty(find(ais_trace,1))
        Xend = Xais(end);
        Yend = Yais(end);
        temp = zeros(size(ais_trace));
        temp(Xend-1:Xend+1, Yend-1:Yend+1) = ais_trace(Xend-1:Xend+1, Yend-1:Yend+1);
        nn = find(temp);
        if length(Xais) > 300
            ais_trace = zeros(size(ais_trace));
        elseif numel(nn) > 1
            q = nn(1);
            [Xnextq, Ynextq] = ind2sub(size(ais_trace), q);
            Xais(end+1) = Xnextq;
            Yais(end+1) = Ynextq;
            ais_trace(Xnextq, Ynextq) = 0;
        elseif numel(nn) == 1
            [Xnext, Ynext] = ind2sub(size(ais_trace), nn);
            Xais(end+1) = Xnext;
            Yais(end+1) = Ynext;
            ais_trace(Xnext, Ynext) = 0;
        elseif numel(nn) == 0
            ais_trace = zeros(size(ais_trace));
        else
            error('something is wrong')
        end
        clear temp nn
    end

    XaisL{loop2} = Xais;
    YaisL{loop2} = Yais;
    but2 = 0;

    if exist('q','var')
       XYq  = (Xais == Xnextq) .* (Yais == Ynextq);
       qdel = find(XYq);
       Xdel = Xais(qdel:end);
       Ydel = Yais(qdel:end);
       ais_traceq(Xdel, Ydel) = 0;
       ais_trace = ais_traceq;
       clear Xais Yais XYq qdel Xdel Ydel q
       Xais = I_skel;
       Yais = J_skel;
       ais_trace(Xais(end), Yais(end)) = 0;
       but2 = 1;
    end
    clear q

    if loop2 > 10
        but2 = 0;
    end
end

for l = 1:loop2
    long(l) = length(XaisL{l});
end
longline = find(long == max(long));
Xais = XaisL{longline};
Yais = YaisL{longline};

xy = [Yais; Xais];
xy = double(xy);

t  = 1:length(xy);
ts = [[1: 0.1: length(xy)], t(end)];
xysm = csaps(t, xy, 0.3, ts);
for u = 2:(length(xysm))
    ax(u) = sqrt((xysm(1,u) - xysm(1,u-1))^2 + (xysm(2,u) - xysm(2,u-1))^2);
end
ax = cumsum(ax);
ax_um = ax .* pixconv;
xys = round(xysm);

double_id = [];
for i = 1:length(xys(1,:))
    for j = (i+1):length(xys(1,:))
        if xys(1,i) == xys(1,j)
            if xys(2,i) == xys(2,j)
                double_id = [double_id; j];
            end
        end
    end
end
double_id = unique(double_id);
alln = 1:length(xys(1,:));
m = setdiff(alln, double_id);
xs = xys(1,:); ys = xys(2,:);
x_pix = xs(m);
y_pix = ys(m);

% The original always has the raw Ch1 z-projection beside the processed file. Some exported
% folders only ship the processed one, and aiscounter falls back to reading the profile from
% it (imaging.load_image); the fallback is mirrored here so those folders can be compared too.
raw_file = [basepath '.tif'];
if exist(raw_file, 'file')
    Ch1_pic = imread(raw_file);
else
    Ch1_pic = pic(:,:,1);
end

for i = 1:length(x_pix)
    lv_c(i) = Ch1_pic(y_pix(i), x_pix(i));
    lv_1(i) = Ch1_pic(y_pix(i)+1, x_pix(i));
    lv_2(i) = Ch1_pic(y_pix(i)-1, x_pix(i));
    lv_3(i) = Ch1_pic(y_pix(i),   x_pix(i)+1);
    lv_4(i) = Ch1_pic(y_pix(i),   x_pix(i)-1);
    lv_5(i) = Ch1_pic(y_pix(i)+1, x_pix(i)+1);
    lv_6(i) = Ch1_pic(y_pix(i)-1, x_pix(i)-1);
    lv_7(i) = Ch1_pic(y_pix(i)+1, x_pix(i)-1);
    lv_8(i) = Ch1_pic(y_pix(i)-1, x_pix(i)+1);
    lv_smooth(i) = mean([lv_c(i) lv_1(i) lv_2(i) lv_3(i) lv_4(i) lv_5(i) lv_6(i) lv_7(i) lv_8(i)]);
end

v = num2str(lv_c); lv_c = str2num(v);

for i = 1:length(x_pix)
    width = (1.5/(pixconv));
    d = ((round(width))+1);
    if i < (d+1)
        lv_slide(i) = mean([lv_c(1:i) lv_c(i:i+d)]);
        lv_ss(i)    = mean([lv_smooth(1:i) lv_smooth(i:i+d)]);
    elseif i > (length(x_pix)-(d+1))
        lv_slide(i) = mean([lv_c(i-d:i) lv_c(i:length(x_pix))]);
        lv_ss(i)    = mean([lv_smooth(i-d:i) lv_smooth(i:length(x_pix))]);
    else
        lv_slide(i) = mean([lv_c(i-d:i) lv_c(i:i+d)]);
        lv_ss(i)    = mean([lv_smooth(i-d:i) lv_smooth(i:i+d)]);
    end
end

norm_lv = (lv_ss - min(lv_ss)) ./ (max(lv_ss)-min(lv_ss));
pix_narray = 1:length(x_pix);

max_i = (find(norm_lv==1)); max_i = max_i(1);
max_x = pix_narray(max_i);
ais_end = find( (pix_narray>max_i) & (norm_lv>f),1,'last');
if length(ais_end)>0
    ais_end = ais_end(1);
else
    ais_end = x_pix(length(x_pix));
end
ais_start = find( (pix_narray<max_i) & (norm_lv<f));
if length(ais_start)>0
    ais_start = ais_start(length(ais_start));
else
    ais_start = 0;
end

debut = ais_start*pixconv;
fin   = ais_end*pixconv;
lngth = fin-debut;

out = struct();
out.click_col      = Xcell;
out.click_row      = Ycell;
out.n_passes       = length(thresholds);   % 2 whenever the rescale fired
out.thresholds     = thresholds;
out.threshold      = threshold;            % the level finally used
out.max_all        = out_max_all;
out.max_ais        = out_max_ais;
out.ratio          = out_max_ais / out_max_all;
out.cc_counts      = cc_counts;
out.labels_sel     = labels_sel;
out.cc_pixels      = cc_pixels;
out.skeleton_count = sum(ais_skeleton(:));
out.seed_row       = I_skel;
out.seed_col       = J_skel;
out.n_candidates   = loop2;
out.cand_lengths   = long;
out.trace_rows     = Xais;
out.trace_cols     = Yais;
out.n_trace        = length(Xais);
out.n_pix          = length(x_pix);
out.x_pix          = x_pix;
out.y_pix          = y_pix;
out.norm_lv        = norm_lv;
out.arclength_um   = ax_um(end);
out.max_i          = max_i;
out.ais_start      = ais_start;
out.ais_end        = ais_end;
out.lngth          = lngth;

fid = fopen(outfile,'wt');
fprintf(fid, '%s', jsonencode(out));
fclose(fid);
fprintf('WROTE %s  passes=%d  thr=%.10f  lngth=%.10f\n', outfile, out.n_passes, threshold, lngth);
end
