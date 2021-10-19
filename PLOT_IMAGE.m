function [h]=PLOT_IMAGE(IM,L)
load dom
load scat

subplot(1,1,1)
colormap(jet)
imax=max(abs(IM(:)));

h=pcolor(dom.y-L/Lambda,dom.x,(abs(IM)/imax));
hold
for i=1:n_scat
  plot((y_scat(i,2)-L)/Lambda+H_y/(2*Lambda),(y_scat(i,1))/Lambda+H_x/(2*Lambda),'wx','linewidth',1.3)
end
hold
axis tight; shading flat;
set(gca,'Fontsize',16)
xlabel('range in \lambda_0')
ylabel('cross-range in \lambda_0')
xlabel('range in \lambda_0')
colorbar
shading flat
return

