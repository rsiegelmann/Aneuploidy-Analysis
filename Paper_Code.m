%% Figure 1a
vhealthy = []
vtrisomy = []
for i = 1:22
r1 = R.s100kb.tpmMeanTrim{i}(:, 1);
r2 = R.s100kb.tpmMeanTrim{i}(:, 2);
for q = 1:length(r1)
    if r1(q)>11750
       gap = r1(q)-r2(q);
       r1(q) = 0;
       r2(q) = 0;
       count = count+1;
    end
end
kh = sum(r1)/length(r1)
vhealthy = [vhealthy; kh]
kv = sum(r2)/length(r2)
vtrisomy = [vtrisomy; kv]
end 
vtrisomy(7) = 2/3*vtrisomy(7)


figure
bar([vhealthy, vtrisomy])
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("RNASeq/100kb", 'FontWeight', 'Bold')
ylim([0 80])
legend("Healthy", "Trisomy")

%% Figure 1b
v = []
count = 0
for i = 1:22
r1 = R.s100kb.tpmMeanTrim{i}(:, 1);
r2 = R.s100kb.tpmMeanTrim{i}(:, 2);
for q = 1:length(r1)
    if r1(q)>11750
       gap = r1(q)-r2(q);   
       r1(q) = 0;
       r2(q) = 0;
       count = count+1;
    end
end
k = sum(r2-r1)/sum(r1);
if i == 7
k = sum(2/3*r2-r1)/sum(r1)
q = k
end
v = [v; k];
end
v = v*100;
m = mean(v);


figure
b = bar(v)
hold on
br = bar(7,v(7),'r');
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("Percent Change in Gene Expression", 'FontWeight', 'Bold')
hold on
yline(m, 'Color', 'r', 'LineWidth', 3)
ylim([-25 25])
hold off

v = []
for i = 1:22
r1 = R.s100kb.tpmMeanTrim{i}(:, 1);
r2 = R.s100kb.tpmMeanTrim{i}(:, 2);
k = sum(r2-r1)/sum(r1);
v = [v; k];
end
v = v.*100;
m = sum(abs(v))/length(v);

%% Figure 1c
ratiomatrix = [];
for num = 1:22
ratiovec = [];
r1 = R.s100kb.tpmMeanTrim{num}(:,1);
r2 = R.s100kb.tpmMeanTrim{num}(:,2);
for i = 1:length(r1)
    if num ~= 7 
    if r1(i) < eps && r2(i) < eps
    elseif r2(i) > 6*r1(i)
        ratiovec = [ratiovec; 6];
    else
        k = r2(i)/r1(i);
        ratiovec = [ratiovec; k];
    end
    end
    if num == 7
        if r1(i) < eps && r2(i) < eps
    elseif 2/3*r2(i) > 6*r1(i)
       ratiovec = [ratiovec; 6];
    else
        k = 2/3*r2(i)/r1(i);
        ratiovec = [ratiovec; k];
    end
    end
end
ratiomatrix{num} = ratiovec;
end

rvec = [];
cvec = [];
groups = [];
for num = 1:22
    for i = 1:length(ratiomatrix{1, num})
        rvec = [rvec; ratiomatrix{1, num}(i)];
        cvec = [cvec; strcat("", num2str(num))];
    end
    groups = [groups, strcat("", num2str(num))]
end
groups = num2cell(groups)
for num = 1:22
    groups{num} = char(groups{num})
end
hold on
violinplot(rvec,cvec, 'GroupOrder', groups);
yline(1, 'k', 'LineWidth', 2)
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("Expression Ratio", 'FontWeight', 'Bold')
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
set(gca, 'YTickLabels', ["0", "1", "2", "3", "4", "5", "6+"])


%% Figure 2a, 2b
inc = zeros(22, 1);
dec = zeros(22, 1);
for num = 1:22
r1 = R.s100kb.tpmMeanTrim{num}(:,1);
r2 = R.s100kb.tpmMeanTrim{num}(:,2);
if num == 7
    r2 = 2/3*r2
end
ratios = zeros(length(r2), 1);
for k = 1:length(ratios)
ratios(k) = (r2(k)/r1(k));
if r1(k)*1.2 < r2(k)
    inc(num) = inc(num)+1;
end
if r1(k) == r2(k)
    ratios(k) = 1;
end
if r1(k) > r2(k)*1.2
    dec(num) = dec(num)+1;
end
end
dec(num) = 100*dec(num)/length(r1);
inc(num) = 100*inc(num)/length(r1);
end

figure
h = bar([inc dec])
color1raw = '#FFA500';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#008000';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
for n = 1:22
    h(2).CData(n,:) = color1;
    h(1).CData(n,:) = color2; 
end
y1 = yline(mean(inc), 'LineWidth', 3);
y2 = yline(mean(dec), 'LineWidth', 3);
y2.Color = color1;  
y1.Color = color2; 
legend(["Increasing", "Decreasing"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("Expression: % Length of Chromosome", 'FontWeight', 'Bold')
ylim([0 40])


inc = zeros(22, 1);
dec = zeros(22, 1);
for num = 1:22
r1 = R.s100kb.tpmMeanTrim{num}(:,1);
r2 = R.s100kb.tpmMeanTrim{num}(:,2);
if num == 7
    r2 = 2/3*r2
end
ratios = zeros(length(r2), 1);
for k = 1:length(ratios)
ratios(k) = (r2(k)/r1(k));
if r1(k)*2 < r2(k)
    inc(num) = inc(num)+1;
end
if r1(k) == r2(k)
    ratios(k) = 1;
end
if r1(k) > r2(k)*2
    dec(num) = dec(num)+1;
end
end
dec(num) = 100*dec(num)/length(r1);
inc(num) = 100*inc(num)/length(r1);
end

figure
h = bar([inc dec])
color1raw = '#FFA500';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#008000';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
for n = 1:22
    h(2).CData(n,:) = color1;
    h(1).CData(n,:) = color2; 
end
y1 = yline(mean(inc), 'LineWidth', 3);
y2 = yline(mean(dec), 'LineWidth', 3);
y2.Color = color1;  
y1.Color = color2; 
legend(["Increasing", "Decreasing"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("Expression: % Length of Chromosome", 'FontWeight', 'Bold')
ylim([0 40])

%% Figure 2c, 2d
vv = [];
vvneg = [];
vvpos = [];
vvboth = [];
vvboth2 = [];
vvdiff = [];
vvdevpos = [];
vvdevneg = [];
vvdev = [];
for k = 0:100
percentpos = [];
percentneg = [];
for num = 1:22
vpos = [];
vneg = [];
r1 = R.s100kb.tpmMeanTrim{num, 1}(:,1);
r2 = R.s100kb.tpmMeanTrim{num, 1}(:,2);
if num == 7
    r2 = 2/3*r2;
end
for i = 1:length(r1)
    if r2(i)/r1(i)>1+k/10
        vpos = [vpos; i];
     elseif r1(i)/r2(i)>1+k/10
        vneg = [vneg; i];
    end
end
percentpos = [percentpos; 100*length(vpos)/length(r1)];
percentneg = [percentneg; 100*length(vneg)/length(r1)];
end

vv = [vv; mean(percentpos+percentneg)];
vvneg = [vvneg; mean(percentneg)];
vvpos = [vvpos; mean(percentpos)];
vvboth = [vvboth; mean(percentneg) mean(percentpos)];
vvboth2 = [vvboth2; mean(percentpos) mean(percentneg)];
vvdiff = [vvdiff;mean(percentpos)-mean(percentneg)];
vvdevpos = [vvdevpos; std(percentpos)];
vvdevneg = [vvdevneg; std(percentneg)];
vvdev = [vvdev; std(percentneg+percentpos)];
end
x = linspace(1, 11, 101);
figure
b = bar(x, vvboth2, 'stacked')
hold on
color1raw = '#FFA500';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#008000';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
b(1).FaceColor = color2;
b(2).FaceColor = color1;
legend(['Increase'; 'Decrease'])
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Minimum Gene Expression Ratio Threshold", 'FontWeight', 'Bold')
ylabel("Average Percent of Chromosome", 'FontWeight', 'Bold')


vv = [];
vvneg = [];
vvpos = [];
vvboth = [];
vvboth2 = [];
vvdiff = [];
vvdevpos = [];
vvdevneg = [];
vvdev = [];
for k = 0:100
percentpos = [];
percentneg = [];
for num = 1:22
vpos = [];
vneg = [];
r1 = R.s100kb.tpmMeanTrim{num, 1}(:,1);
r2 = R.s100kb.tpmMeanTrim{num, 1}(:,2);
if num == 7
    r2 = 2/3*r2
end
for i = 1:length(r1)
    if r2(i)/r1(i)>1+k/10 && r2(i)/r1(i)<1.1+k/10
        vpos = [vpos; i];
     elseif r1(i)/r2(i)>1+k/10 && r1(i)/r2(i)<1.1+k/10
        vneg = [vneg; i];
    end
    if k == 100
        if r2(i)/r1(i)>1+k/10
        vpos = [vpos; i];
     elseif r1(i)/r2(i)>1+k/10
        vneg = [vneg; i];
    end
    end
end
percentpos = [percentpos; 100*length(vpos)/length(r1)];
percentneg = [percentneg; 100*length(vneg)/length(r1)];
end

vv = [vv; mean(percentpos+percentneg)];
vvneg = [vvneg; mean(percentneg)];
vvpos = [vvpos; mean(percentpos)];
vvboth = [vvboth; mean(percentneg) mean(percentpos)];
vvboth2 = [vvboth2; mean(percentpos) mean(percentneg)];
vvdiff = [vvdiff;mean(percentpos)-mean(percentneg)];
vvdevpos = [vvdevpos; std(percentpos)];
vvdevneg = [vvdevneg; std(percentneg)];
vvdev = [vvdev; std(percentneg+percentpos)];
end
x = linspace(1, 11, 101);
figure
b = bar(x, vvboth2, 'stacked')
hold on
color1raw = '#FFA500';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#008000';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
b(1).FaceColor = color2;
b(2).FaceColor = color1;
legend(['Increase'; 'Decrease'])
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Gene Expression Ratio Threshold", 'FontWeight', 'Bold')
ylabel("Average Percent of Chromosome", 'FontWeight', 'Bold')

%% Figure 3a

count = zeros(22, 2)
for num = 1:22
    pc = 0;
    nc = 0;
    A1 = H.s100kb.oeTrim{num}(:,:,1);
    A2 = H.s100kb.oeTrim{num}(:,:,2);
    pc1 = pca(A1);
    pc1 = pc1(:,1);
    pc2 = pca(A2);
    pc2 = pc2(:,1);
    if num == 4 || num == 5
        pc1 = pca(A1);
        pc1 = pc1(:,2);
        pc2 = pca(A2);
        pc2 = pc2(:,2);
    end
    if corr(pc1, pc2) < 0
    pc1 = pc1*(-1);
    end
    for i=1:size(pc1,1)-1
        if pc2(i)>0 && pc1(i)<0
            pc = pc+1;
        end
        if pc2(i)<0 && pc1(i)>0
            nc = nc+1;
        end
    end
    count(num, :) = [pc/length(A1), nc/length(A1)];
end
hold on
color1raw = '#FFA500';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#008000';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
figure
h = bar([100*count(:, 1) 100*count(:, 2)])
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
for n = 1:22
    h(1).CData(n,:) = color2;
    h(2).CData(n,:) = color1; 
end
y1 = yline(100*mean(count(:, 1)), 'LineWidth', 3);
y2 = yline(100*mean(count(:, 2)), 'LineWidth', 3);
y1.Color = color2;  
y2.Color = color1; 
legend(["B to A", "A to B"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome Number", 'FontWeight', 'Bold')
ylabel("Percent Compartment Change", 'FontWeight', 'Bold')

%% Figure 3b

MC = []
for num = 1:22
    A1 = H.s100kb.krTrim{num}(:,:,1);
    A2 = H.s100kb.krTrim{num}(:,:,2);
    [FN1,Fvec1] = Fiedlervec(A1);
    [FN2,Fvec2] = Fiedlervec(A2);
    MC = [MC;100*(FN2-FN1)/FN1];
end
bar(MC)
hold on
bar(7,MC(7),'r');
yline(mean(MC),'r', 'LineWidth', 3)
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0])
xlabel("Chromosome", 'FontWeight', 'Bold')
ylabel("Percent Change Algebraic Connectivity", 'FontWeight', 'Bold')

%% Figure 4a, Supplementary Figure 2

for num = 1:22
    t1 = graph(H.s100kb.oeTrim{num}(:,:,1));
    t2 = graph(H.s100kb.oeTrim{num}(:,:,2));
    bet1 = centrality(t1, 'betweenness');
    bet2 = centrality(t2, 'betweenness');
    close1 = centrality(t1, 'closeness');
    close2 = centrality(t2, 'closeness');   
    deg1 = centrality(t1, 'degree');
    deg2 = centrality(t2, 'degree');
    eigen1 = centrality(t1, 'eigenvector');
    eigen2 = centrality(t2, 'eigenvector');
    figure
    legend('Healthy','Trisomy-7', 'Location', 'nw');
    s1 = subplot(4,3,1);
    bar(bet1, 'b')
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Betweenness Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    ylabel("Betweenness", 'FontWeight', 'Bold')
    title("Healthy", 'FontWeight', 'Bold')
    set(gca,'TickLength',[0 0])
    s2 = subplot(4,3,4);
    bar(close1, 'b');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Closeness Centrality', 'FontSize', 16,'FontName', 'Times New Roman'); 
    ylabel("Closeness", 'FontWeight', 'Bold')
    set(gca,'TickLength',[0 0])
    s3 = subplot(4,3,7);
    bar(deg1, 'b');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Degree Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    ylabel("Degree", 'FontWeight', 'Bold')
    set(gca,'TickLength',[0 0])
    s4 = subplot(4,3,10);
    bar(eigen1, 'b');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Eigenvector Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    ylabel("Eigenvector", 'FontWeight', 'Bold')
    set(gca,'TickLength',[0 0])
    s5 = subplot(4,3,2);
    bar(bet2, 'r');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Betweenness Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    set(gca,'TickLength',[0 0])
    title("Trisomy-7", 'FontWeight', 'Bold')
    s6 = subplot(4,3,5);
    bar(close2, 'r');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Closeness Centrality', 'FontSize', 16,'FontName', 'Times New Roman'); 
    set(gca,'TickLength',[0 0])
    s7 = subplot(4,3,8);
    bar(deg2, 'r');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Degree Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    set(gca,'TickLength',[0 0])
    s8 = subplot(4,3,11);
    bar(eigen2, 'r');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', [])
    %title('Eigenvector Centrality', 'FontSize', 16,'FontName', 'Times New Roman') ;
    set(gca,'TickLength',[0 0])
    s9 = subplot(4,3,3);
    bar(100*(bet2./bet1-1), 'm');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[])
    %title('Betweenness Centrality Change', 'FontSize', 16,'FontName', 'Times New Roman');
    set(gca,'TickLength',[0 0])
    title("Percent Difference", 'FontWeight', 'Bold')
    ylim([-50 100])
    s10 = subplot(4,3,6);
    bar(100*(close2./close1-1), 'm');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[])
    %title('Closeness Centrality Change', 'FontSize', 16,'FontName', 'Times New Roman');
    set(gca,'TickLength',[0 0])
    ylim([-50 100])
    s11 = subplot(4,3,9);
    bar(100*(deg2./deg1-1), 'm');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[])
    %title('Degree Centrality Change', 'FontSize', 16,'FontName', 'Times New Roman');
    set(gca,'TickLength',[0 0])
    ylim([-50 100])
    s12 = subplot(4,3,12);
    bar(100*(eigen2./eigen1-1), 'm');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[])
    %title('Eigenvector Centrality Change', 'FontSize', 16,'FontName', 'Times New Roman');
    set(gca,'TickLength',[0 0])
    ylim([-50 100])
    linkaxes([s1 s5], 'xy');
    linkaxes([s2 s6], 'xy');
    linkaxes([s3 s7], 'xy');
    linkaxes([s4 s8], 'xy');
    linkaxes([s9 s10 s11 s12], 'xy');
end

%% Figure 4b, Supplementary Figure 1
figure
for num = 1:22
figure
A1 = H.s100kb.oeTrim{num}(:,:,1);
s1 = subplot(1,3,1);
imagesc(log2(A1), [-5.5 5.5]);
hold on 
TAD_interval1 = TAD_Laplace_Sijia(A1,0.75,3);
    for i=1:size(TAD_interval1, 1)-2
        line([TAD_interval1(i),TAD_interval1(i)],[TAD_interval1(i),TAD_interval1(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
        line([TAD_interval1(i),TAD_interval1(i+1)],[TAD_interval1(i),TAD_interval1(i)], 'Color', 'blue', 'LineWidth', 1.5);
         line([TAD_interval1(i+1),TAD_interval1(i+1)],[TAD_interval1(i),TAD_interval1(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
        line([TAD_interval1(i),TAD_interval1(i+1)],[TAD_interval1(i+1),TAD_interval1(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
    end

axis square   
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', []);
%title("Chromosome " + num2str(num) + " (Healthy)")
xlabel("# TADs = " + num2str(size(TAD_interval1, 1)), 'FontWeight', 'Bold');
set(gca, 'Colormap', hot);
temp = get(gca,'Colormap');
temp = flipud(temp);
set(gca, 'Colormap', temp);
A2 = H.s100kb.oeTrim{num}(:,:,2);
s2 = subplot(1,3,2);
imagesc(log2(A2),  [-5.5 5.5]);
hold on 
TAD_interval2 = TAD_Laplace_Sijia(A2,0.75,3);
    for i=1:size(TAD_interval2, 1)-2
        line([TAD_interval2(i),TAD_interval2(i)],[TAD_interval2(i),TAD_interval2(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
        line([TAD_interval2(i),TAD_interval2(i+1)],[TAD_interval2(i),TAD_interval2(i)], 'Color', 'blue', 'LineWidth', 1.5);
         line([TAD_interval2(i+1),TAD_interval2(i+1)],[TAD_interval2(i),TAD_interval2(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
        line([TAD_interval2(i),TAD_interval2(i+1)],[TAD_interval2(i+1),TAD_interval2(i+1)], 'Color', 'blue', 'LineWidth', 1.5);
    end
set(gca, 'Colormap', hot)
temp = get(gca,'Colormap');
temp = flipud(temp);
set(gca, 'Colormap', temp)
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', []);
%title("Chromosome " + num2str(num) + " (Trisomy)")
xlabel("# TADs = " + num2str(size(TAD_interval2, 1)), 'FontWeight', 'Bold')
axis square
s3 = subplot(1,3,3);
imagesc(log2(A2./A1), [-5.5 5.5]);
set(gca, 'Colormap', hot);
temp = get(gca,'Colormap');
temp = flipud(temp);
set(gca, 'Colormap', temp);
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman','XTick',[], 'YTick', []);
%title("Difference Hi-C of Chromosome" + num)
axis square
end

%% Figure 4c

mean1 = zeros(22,1);
mean2 = zeros(22,1);
v1 = [];
v2 = [];
for num = 1:22
A1 = H.s100kb.oeTrim{num}(:,:,1);
TAD_interval1 = TAD_Laplace_Sijia(A1,0.75 ,3);
v1 = [v1; length(TAD_interval1)];
A2 = H.s100kb.oeTrim{num}(:,:,2);
TAD_interval2 = TAD_Laplace_Sijia(A2,0.75,3);
v2 = [v2; length(TAD_interval2)];
TAD1 = TAD_interval1(2:end)-TAD_interval1(1:end-1);
TAD2 = TAD_interval2(2:end)-TAD_interval2(1:end-1);
mean1(num) = sum(TAD1,1)/length(TAD1);
mean2(num) = sum(TAD2,1)/length(TAD2);
end
v1 = transpose(v1(:, 1));
v2 = transpose(v2(:, 1));
grapher = [];
for int = 1:22
    grapher = [grapher; 100*(mean2(int)-mean1(int))/mean1(int) 100*(v2(int)-v1(int))/v1(int)];
end
figure
b = bar(grapher, 'stacked')
hold on
color1raw = '#77D8C0';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#FF8C69';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
b(1).FaceColor = color1;
b(2).FaceColor = color2;
y1 = yline(100*mean((v2-v1)./v1), 'LineWidth', 3);
y2 = yline(100*mean((mean2-mean1)./mean1), 'LineWidth', 3);
y1.Color = color2;  
y2.Color = color1; 
legend(["Size", "Number"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0], 'Color', 'k')
xlabel("Chromosome", 'FontWeight', 'Bold')
ylabel("Percent Change in TAD", 'FontWeight', 'Bold')
ylim([-35 35])

%% Figures 5b and 5c 
v1 = [];
v2 = [];
v3 = [];
v4 = [];
for num = 1:22
    r1 = R.s100kb.tpmMeanTrim{num}(:,1);
    r2 = R.s100kb.tpmMeanTrim{num}(:,2);
    t1 = graph(H.s100kb.krTrim{num}(:,:,1));
    t2 = graph(H.s100kb.krTrim{num}(:,:,2));
    bet1 = centrality(t1, 'betweenness');
    bet2 = centrality(t2, 'betweenness');
    close1 = centrality(t1, 'closeness');
    close2 = centrality(t2, 'closeness');
    deg1 = centrality(t1, 'degree');
    deg2 = centrality(t2, 'degree');
    eigen1 = centrality(t1, 'eigenvector');
    eigen2 = centrality(t2, 'eigenvector');
    v1 = [v1; 100*[mean(bet2)/mean(bet1)-1]];
    v2 = [v2; 100*[mean(close2)/mean(close1)-1]];
    v3 = [v3; 100*[mean(deg2)/mean(deg1)-1]];
    v4 = [v4; 100*[mean(eigen2)/mean(eigen1)-1]];

end


grapher1 = zeros(22, 2);
for i = 1:22
    grapher1(i, :) = [v3(i), v1(i)];
end

figure
b = bar(grapher1, 'stacked')
hold on
color1raw = '#77D8C0';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#FF8C69';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
b(1).FaceColor = color1;
b(2).FaceColor = color2;
y1 = yline(mean(v1), 'LineWidth', 3);
y2 = yline(mean(v3), 'LineWidth', 3);
y1.Color = color2;  
y2.Color = color1; 
legend(["Degree", "Betweenness"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0], 'Color', 'k')
xlabel("Chromosome", 'FontWeight', 'Bold')
ylabel("Average Percent Change in Centralities", 'FontWeight', 'Bold')
ylim([-50 50])

grapher2 = zeros(22, 2);
for i = 1:22
    grapher2(i, :) = [v2(i), v4(i)];
end

figure
b = bar(grapher2, 'stacked')
hold on
color1raw = '#77D8C0';
color1 = sscanf(color1raw(2:end),'%2x%2x%2x',[1 3])/255;
color2raw = '#FF8C69';
color2 = sscanf(color2raw(2:end),'%2x%2x%2x',[1 3])/255;
b(1).FaceColor = color1;
b(2).FaceColor = color2;
y1 = yline(mean(v4), 'LineWidth', 3);
y2 = yline(mean(v2), 'LineWidth', 3);
y1.Color = color2;  
y2.Color = color1; 
legend(["Closeness", "Eigenvector"])
hold off
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman')
set(gca,'TickLength',[0 0], 'Color', 'k')
xlabel("Chromosome", 'FontWeight', 'Bold')
ylabel("Average Percent Change in Centralities", 'FontWeight', 'Bold')
ylim([-50 50])

