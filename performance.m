% performance.m
% Calculates speedup for RQI-PCGTLS-MP with mxn sized matrices and number r
% of Rayleigh quotgient iterations performed in total. Thsi theoretical
% calculation assume sthat PCGTLS performs exactly k+1 iterations in kth
% Rayleigh quotient itertaion.

% Size of mxn matrices
n = [1e+2,1e+3,1e+4,1e+5,1e+6,1e+7,1e+8,1e+9];
m = [1e+2,1e+3,1e+4,1e+5,1e+6,1e+7,1e+8,1e+9];

% Coefficients for the precision setting. xh for fp16, xs for fp32, xd for fp64.
xh = 0.25; xs = 0.5; xd = 1;

% Number r of total Rayleigh quotient iterations performed
r_tot = [20,50,100,500,1000,2000];

for rr = 1:length(r_tot)
    r = r_tot(rr);
    
    for i = 1:length(m)
        for j = 1:i

            % cost of the parts calculated using u_p
            costu_p = 2*r*(n(j)^2+2*n(j)-1)+(r^2+3*r)*(2*n(j)^2+14*n(j)-3);    

            % cost of the parts calculated using u_q (using householder qr factorization)
            costu_q = 2*m(i)*n(j)^2-2*(n(j)^3)/3;

            % cost of the parts calculated using u
            costu = 2*m(i)*n(j)+3*m(i)+4*n(j)+2*n(j)^2-2+r*(4*m(i)*n(j)+6*m(i)+11*n(j)-5);

            % total cost (multiply costu_p, costu_q, and costu with the corresponding constant xh, xs, or xd)
            cost_total_mp= xd*costu+xs*costu_p+xh*costu_q;

            % total cost of the uniform precision algorithm calculated in fp64
            cost_total_uniform= xd*(costu_q+costu+costu_p+costu_q);

            % calculate speedup
            speedup(i,j) = cost_total_uniform/cost_total_mp;
        end
    end

    % Plot heatmap
    xvalues = {'10^2','10^3','10^4','10^5','10^6','10^7','10^8','10^9'};
    yvalues = {'10^9','10^8','10^7','10^6','10^5','10^4','10^3','10^2'};
    h = heatmap(xvalues,yvalues,flipud(speedup))
    k = '$(u,u_p,u_q)$ = (double,single,half),';
    h.Title = strcat(k,' $r =$',num2str(r));
    h.XLabel = 'n';
    h.YLabel = 'm';
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).Title.Interpreter = 'latex';
    h.ColorLimits = [0,4];
    set(gca,'FontSize',14)
end