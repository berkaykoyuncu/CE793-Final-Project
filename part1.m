%QUESTION 1 

x=linspace(0,1,50);
y=linspace(0,1,50);
nx=50;
ny=50;     

nt=1000;

dt=0.001;

dx=x(2)-x(1);
dy=y(2)-y(1);

vis=0.015;

u=zeros(nx,ny);
un=zeros(nx,ny);
v=zeros(nx,ny);
vn=zeros(nx,ny);

%Initial Conditions
for i = 1:nx;
    for j = 1:ny;
        u(i,j) = 0;
        v(i,j) = 0;
    end
end

%Plus, boundary conditions are mandated
u(1,:)=0;
u(nx,:)=0;
u(:,1)=sin(2*pi.*x);
u(:,ny)=sin(2*pi.*x);
v(:,1)=0;
v(:,ny)=1;
v(1,:)=1-y;
v(nx,:)=1-y;

i=2:nx-1;
j=2:ny-1;

for it=0:nt
    un=u;
    vn=v;
    %plotting the velocity field
    h=quiver(x,y,u,v,'k');
    axis([0 1 0 1])
    axis square
    title({['2-D Burgers'' equation with {\nu} = ',num2str(vis)];'Vector field = (u and v)';['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Spatial co-ordinate (y) \rightarrow')
    drawnow;
    refreshdata(h)
    u(i,j)=un(i,j)-(dt*(un(i,j)-un(i-1,j)).*un(i,j)/dx)-(dt*(un(i,j)-un(i,j-1)).*vn(i,j)/dy)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j-1)-2*un(i,j)+un(i,j+1))/(dy*dy));
    v(i,j)=vn(i,j)-(dt*(vn(i,j)-vn(i-1,j)).*un(i,j)/dx)-(dt*(vn(i,j)-vn(i,j-1)).*vn(i,j)/dy)+(vis*dt*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))/(dx*dx))+(vis*dt*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/(dy*dy));
end

figure(2);
plot(u(20,:), 'LineWidth', 2)
hold on
plot(v(20,:), 'LineWidth', 2)
title('u vs v at the middle slice (25,:)')
legend('u explicit', 'v explicit')
hold off

% N=100;
% x = linspace(0,1,N);
% y = linspace(0,1,N);
% dt = 0.001;
% vis = 0.015;
% del_x = x(2)-x(1);            %Width of space step(x)
% del_y = y(2)-y(1);            %Width of space step(y)
% 
% u=zeros(N,N);
% v=zeros(N,N);
% 
% u_new = zeros(N,N);
% v_new = zeros(N,N);
% 
% %BOUNDARY CONDITIONS
% u(1,:)=0;
% u(N,:)=0;
% u(:,1)=sin(2*pi.*x);
% u(:,N)=sin(2*pi.*x);
% v(:,1)=0;
% v(:,N)=1;
% v(1,:)=1-y;
% v(N,:)=1-y;
% 
% total_error = 100;
% for iteration=1:1000
%     for i=2:N-1
%         for j=2:N-1
%             u_new(i,j) = u(i,j) - (1/(2*del_x))*(dt*((u(i+1,j) - u(i-1,j)))*u(i,j)) - (1/(2*del_y))*(dt*((u(i,j+1)-u(i,j-1)))*v(i,j))+(vis*dt*(u(i+1,j)-2*u(i,j)+u(i-1,j))/(del_x^2))+(vis*dt*(u(i,j-1)-2*u(i,j)+u(i,j+1))/(del_y^2));
%             v_new(i,j) = v(i,j) - (1/(2*del_x))*(dt*((v(i+1,j) - v(i-1,j)))*u(i,j)) - (1/(2*del_y))*(dt*((v(i,j+1)-v(i,j-1)))*v(i,j))+(vis*dt*(v(i+1,j)-2*v(i,j)+v(i-1,j))/(del_x^2))+(vis*dt*(v(i,j-1)-2*v(i,j)+v(i,j+1))/(del_y^2));
%         end
%         for i=2:N-1
%             for j=2:N-1
%                 u(i,j) = u_new(i,j);
%                 v(i,j) = v_new(i,j);
%             end
%         end
%     end
% end
% 
% 
% % for i=2:N-1
% %     for j=2:N-1
% %         u_new(i,j) = u(i,j);
% %         v_new(i,j) = v(i,j);
% %     end
% % end
% 
% % time = 0.0005;
% 
% % for i = 1:100
% %     for j = 1:100
% %         u_exact(i,j) = (x(i)+y(j)-(2.*x(i).*time)) ./ (1-(2.*time.^2));
% %         v_exact(i,j) = (x(i)-y(j)-(2.*y(j).*time)) ./ (1-(2.*time.^2));
% %     end
% % end
% 
% for it=0:100
%     %plotting the velocity field
%     h=quiver(x,y,u',v','k');
%     axis([0 1 0 1])
%     axis square
%     title({['2-D Burgers'' equation with {\nu} = ',num2str(vis)];'Transport property vector field {\bfu}=(u_x,u_y)';['time(\itt) = ',num2str(dt*it)]})
%     xlabel('Spatial co-ordinate (x) \rightarrow')
%     ylabel('Spatial co-ordinate (y) \rightarrow')
%     drawnow;
%     refreshdata(h)
% end
% 
% % figure(1);
% % quiver(x,y,u,v, 'k')
% % hold on
% % title('Surface plot of u and v [quiver]')
% % xlabel('x')
% % ylabel('y')
% % hold off
% 
% % figure(2);
% % quiver(x,y,u_exact,v_exact, 'k')
% % hold on
% % title('Surface plot of u and v [Exact]')
% % xlabel('x')
% % ylabel('y')
% % hold off