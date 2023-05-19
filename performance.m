n = [1e+2,1e+3,1e+4,1e+5,1e+6,1e+7,1e+8,1e+9];
m = [1e+2,1e+3,1e+4,1e+5,1e+6,1e+7,1e+8,1e+9];
xh = 0.25; xs = 0.5; xd = 1; xq = 2;
r_tot = [5,10,20,50,70,73,74,100];

for rr = 1:length(r_tot)
    r = r_tot(rr);
    p = ((r+1)*(r+2)/2)-1;
    
    for i = 1:length(m)
        for j = 1:i
            costq = 2*m(i)*n(j)+m(i)+2*m(i)*n(j)*r+m(i)*r;
            costs = 2*r*(n(j)^2+2*n(j)-1)+(r^2+3*r)*(2*n(j)^2+14*n(j)-3);%2*(r*n(j)^2+p*(2*n(j)^2+14*n(j)-3));
            costh = 2*m(i)*n(j)^2-2*(n(j)^3)/3; %2*m(i)*n(j)^2-n^2+(n^3)/3; %for cholesky+A^TA
            costd = 2*m(i)+4*n(j)+2*n(j)^2-2+r*(2*m(i)*n(j)+5*m(i)+11*n(j)-5);
            cost_total_mp= xq*costq+xd*costd+xs*costs+xh*costh;
            cost_total_uniform= xd*(costq+costd+costs+costh);
            speedup(i,j) = cost_total_uniform/cost_total_mp;
        end
    end
    xvalues = {'10^2','10^3','10^4','10^5','10^6','10^7','10^8','10^9'};
    yvalues = {'10^9','10^8','10^7','10^6','10^5','10^4','10^3','10^2'};
    h = heatmap(xvalues,yvalues,flipud(speedup))
    k = '$(u_r,u,u_p,u_q)$ = (quad,double,single,half),';
    h.Title = strcat(k,' $r =$',num2str(r));
    h.XLabel = 'n';
    h.YLabel = 'm';
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).Title.Interpreter = 'latex';
    h.ColorLimits = [0,4];
    set(gca,'FontSize',14)
    
    snbase = strcat('figs/performance_');
    savename = strcat(snbase,num2str(r));
    savefig(strcat(savename,'.fig'));
    saveas(gcf, strcat(savename,'.pdf'));
    close all
end