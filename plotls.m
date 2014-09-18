function h_plot = plotls(x,y,inds)
linespec = {'-','--','-.','..'};
for ip = 1:length(inds)
    h_plot(ip) = plot(x(:,inds(ip)),y,linespec{ip});
end