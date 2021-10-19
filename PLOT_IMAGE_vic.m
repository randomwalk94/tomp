function [h]=PLOT_IMAGE_vic(dat,d,dx,dy)
load dom
load scat

ld = length(dat);

colorp = [rand(1),rand(1),rand(1)];

for i=1:ld
    if dat(i)>0
        [xk,yk] = index_to_image(i,d);
        plot(-dom.Lxd/2+(xk-1)*dx,-dom.Lyd/2+(yk-1)*dy,'+','Color',[colorp(1),colorp(2),colorp(3)],'linewidth',1)
        if dat(i)>1
            plot(-dom.Lxd/2+(xk-1)*dx,-dom.Lyd/2+(yk-1)*dy,'*','Color',[colorp(1),colorp(2),colorp(3)],'linewidth',2)
        end
    end
end
return

