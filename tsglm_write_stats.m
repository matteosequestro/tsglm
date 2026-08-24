function clusters_cut  = tsglm_write_stats(modelout)


modelout.alpha = 0.05;
clusters = modelout.obs_clusters_sum;
clust_par_means = tsglm_clusters_parmeans(modelout);
bfs = tsglm_plot_clust_par_means(clust_par_means(:, 2:end), "noplot", 0);

clusters.BF01 = 1./bfs.BF10;
clusters.BF10 = bfs.BF10;

clusters_cut = clusters(:, ["pred", "clusterstat", "prob", "BF01", "BF10", "length", "first", "last"]);
clusters_cut.Properties.VariableNames = ["pred", "clusterstat", "pcorr", "BF01", "BF10", "length", "first", "last"];

% print rsults
disp(clusters_cut)





end