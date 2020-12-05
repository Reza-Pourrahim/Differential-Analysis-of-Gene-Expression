close all
clear all
clc
%% Load data
adjCN_pos = csvread('adj_matrix_pos.csv',1,0);
adjCN_neg = csvread('adj_matrix_neg.csv',1,0);
adjCN_unsign = csvread('adj_matrix_unsigned.csv',1,0);
[num,gene_id] = xlsread('TCGA-BLCA_selected_gene_ids.csv');
gene_id(1,:) = [];
CentrMeasures = {'degree','closeness','betweenness','eigenvector'};
%% Centrality Measures for Binary and signed network: Subnetwork with only Positive Links
gene_symDe= gene_id;
K_pos = sum(adjCN_pos,2);
adjCN_pos((K_pos==0),:) = []; adjCN_pos(:,(K_pos==0)) = [];

G_pos = graph(adjCN_pos);
gene_symDe((K_pos==0)) = [];
clear adjCN_pos K_pos
for i=1:length(CentrMeasures)
    clear index_pos Y_pos
    index_pos = centrality(G_pos,CentrMeasures{i});
    Y_pos = prctile(index_pos,95);
    figure(i); clf; h2_pos = plot(G_pos,'Layout','force','UseGravity',true);
    highlight(h2_pos,(index_pos>Y_pos),'NodeColor','r','Marker','h','MarkerSize',4)
    type=CentrMeasures{i};
    title(sprintf('%s %s',type,'- Positive subnetwork'));
    Hubs_pos.(CentrMeasures{i}) = gene_symDe(index_pos>Y_pos);
end
%% Centrality Measures for Binary and signed network: Subnetwork with only Negative Links
gene_symDe= gene_id;
K_neg = sum(adjCN_neg,2);
adjCN_neg((K_neg==0),:) = []; adjCN_neg(:,(K_neg==0)) = [];

G_neg = graph(adjCN_neg);
gene_symDe((K_neg==0)) = [];
clear adjCN_neg K_neg

for j=1:length(CentrMeasures)
    clear index_neg Y_neg
    index_neg = centrality(G_neg,CentrMeasures{j});
    Y_neg = prctile(index_neg,95);
    figure(j+4);clf; h2_neg = plot(G_neg,'Layout','force','UseGravity',true);
    highlight(h2_neg,(index_neg>Y_neg),'NodeColor','r','Marker','h','MarkerSize',4)
    type=CentrMeasures{j};
    title(sprintf('%s %s',type,'- Negative subnetwork'));
    Hubs_neg.(CentrMeasures{j}) = gene_symDe(index_neg>Y_neg);
end
%% Centrality Measures for Binary and unsigned network
gene_symDe= gene_id;
K_unsign = sum(adjCN_unsign,2);
adjCN_unsign((K_unsign==0),:) = []; adjCN_unsign(:,(K_unsign==0)) = [];

G_unsign = graph(adjCN_unsign);
gene_symDe((K_unsign==0)) = [];
% mFCDe((K==0)) = [];
clear adjCN_unsign K_unsign

for k=1:length(CentrMeasures)
    clear index_unsign Y_unsign
    index_unsign = centrality(G_unsign,CentrMeasures{k});
    Y_unsign = prctile(index_unsign,95);
    figure(k+8); h2_unsign = plot(G_unsign,'Layout','force','UseGravity',true);
    highlight(h2_unsign,(index_unsign>Y_unsign),'NodeColor','r','Marker','h','MarkerSize',4)
    type=CentrMeasures{k};
    title(sprintf('%s %s',type,'- Unsigned Subnetwork'));
    Hubs_unsign.(CentrMeasures{k}) = gene_symDe(index_unsign>Y_unsign);
end