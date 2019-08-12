function plot_bumpy(Nc,n,N,Dn,x,y,th,r_shape,th_shape,Lx,Ly)

color_list=zeros(Nc,3);
for colork=1:Nc
    color_list(colork,:)=[rand rand rand];
end

x_vert=flat((repmat(x',[1,n])+r_shape.*cos(repmat(th',[1,n])+th_shape))');
y_vert=flat((repmat(y',[1,n])+r_shape.*sin(repmat(th',[1,n])+th_shape))');
D_vert=flat((repmat(Dn',[1,n]))');

figure(1);
clf;
h=zeros(N,3,3);
cc=jet(N);

for np=1:N
    n1=ceil(np/n);
    for nn=-1:1
        for mm=-1:1;
            h(np,nn+2,mm+2)=rectangle('Position',[x_vert(np)-.5*D_vert(np)+nn*Lx+0*mm*Ly y_vert(np)-.5*D_vert(np)+mm*Ly D_vert(np) D_vert(np)],'Curvature',[1 1],'edgecolor',color_list(n1,:),'facecolor',color_list(n1,:));
        end
    end
end

for np=1:Nc
    for nn=-1:1
        for mm=-1:1;
            he(np,nn+2,mm+2)=viscircles([x(np),y(np)],r_shape(np,1),'Color',color_list(np,:),'LineStyle','-','EnhanceVisibility',false);
%                 he(np,nn+2,mm+2)=patch(x_vert((1:n)+(np-1)*n), y_vert((1:n)+(np-1)*n),color_list(np,:),'EdgeColor','none');
        end
    end
end
hold on
%rectangle('Position',[0 0 Lx Ly],'Curvature',[0 0],'EdgeColor','k');
axis('equal'); box on; set(gca,'XTick',[],'YTick',[])
%     axis([-.5*Lx 1.5*Lx -.5*Ly 1.5*Ly]);
hb=plot([0 Lx Lx+0 +0 0],[0 0 Ly Ly 0],'k');
%     hc=plot([0 Lx Lx+0 +0 0],[0 0 Ly Ly 0],'k');

axis('equal');
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

% for i=1:size(contact_list,1)
%     line([x(contact_list(i,1)),x(contact_list(i,2))],[y(contact_list(i,1)),y(contact_list(i,2))],'Color','black','LineStyle','--');
% end

cr=false;
if (cr)
    mf='jamming_movie.avi';
    mov=VideoWriter(mf);
    open(mov);
end



end

