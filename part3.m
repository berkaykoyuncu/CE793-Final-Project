clear all
%FINAL PROJECT
%QUESTION 3 - Time-Accurate Manner

x=linspace(0,1,100);
y=linspace(0,1,100);
nx=100;
ny=100;     

nt=1000;

dt=0.001;

dx=x(2)-x(1);
dy=y(2)-y(1);

vis=0.015;

u=zeros(nx,ny);
un=zeros(nx,ny);
v=zeros(nx,ny);
vn=zeros(nx,ny);

nt = 12;   %FINAL TIME STEP

%INITIAL CONDITIONS
for i = 1:nx;
    for j = 1:ny;
        u(i,j) = 0;
        v(i,j) = 0;
    end
end

%BOUNDARY CONDITIONS
for time=0:dt:nt
    for j=1:ny
        u(1,j) = y(j)./(1 - 2.*time.^2);
    end
    for j=1:ny
        u(ny,j) = (1 + y(j) - (2.*time ))/(1 - 2.*time.^2);
    end
    for i=1:nx
        u(i,1) = (x(i) - 2.*x(i).*time) / (1 - 2.*time.^2);
    end
    for i=1:nx
        u(i,nx) = (x(i) + 1 - (2.*x(i).*time)) / (1 - 2.*time.^2);
    end
    for j=1:ny
        v(1,j) = (-y(j) - 2.*y(j).*time) / (1 - 2.*time.^2);
    end
    for j=1:nx
        v(nx,j) = (1 - y(j) - 2.*y(j).*time) / (1 - 2.*time.^2);
    end
    for i=1:nx
        v(i,1) = x(i) / (1 - 2.*time^2);
    end
    for i=1:ny
        v(i,ny) = (x(i) - 1 - 2.*time) / (1 - 2.*time.^2);
    end
end

%EXPLICIT SCHEME

i=2:nx-1;
j=2:ny-1;

for it=0:dt:nt
    un=u;
    vn=v;
    u(i,j)=un(i,j)-(dt*(un(i,j)-un(i-1,j)).*un(i,j)/dx)-(dt*(un(i,j)-un(i,j-1)).*vn(i,j)/dy)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j-1)-2*un(i,j)+un(i,j+1))/(dy*dy));
    v(i,j)=vn(i,j)-(dt*(vn(i,j)-vn(i-1,j)).*un(i,j)/dx)-(dt*(vn(i,j)-vn(i,j-1)).*vn(i,j)/dy)+(vis*dt*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))/(dx*dx))+(vis*dt*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/(dy*dy));
end

%EXACT SOLUTION

for time=0:dt:nt
    for i = 1:nx
        for j = 1:ny
            u_exact(i,j) = (x(i)+y(j)-(2.*x(i).*time)) ./ (1-(2.*time.^2));
            v_exact(i,j) = (x(i)-y(j)-(2.*y(j).*time)) ./ (1-(2.*time.^2));
        end
    end
end

figure(1);
plot(u(50,:), 'LineWidth', 2)
hold on
plot(u_exact(50,:), 'LineWidth', 2)
title('u vs u_{exact} at the 50^{th} slice')
legend('u (explicit)', 'u (exact)')
hold off

figure(2);
plot(v(50,:), 'LineWidth', 2)
hold on
plot(v_exact(50,:), 'LineWidth', 2)
title('v vs v_{exact} at the 50^{th} slice')
legend('v (explicit)', 'v (exact)')
hold off

% u_new = u;
% v_new = v;
% 
% 
% total_error = 100;
% 
% for time=0:dt:1
%     for j=1:N
%         u(1,j) = y(j)./(1 - 2.*time.^2);
%     end
%     for j=1:N
%         u(N,j) = (1 + y(j) - (2.*time ))/(1 - 2.*time.^2);
%     end
%     for i=1:N
%         u(i,1) = (x(i) - 2.*x(i).*time) / (1 - 2.*time.^2);
%     end
%     for i=1:N
%         u(i,N) = (x(i) + 1 - (2.*x(i).*time)) / (1 - 2.*time.^2);
%     end
%     for j=1:N
%         v(1,j) = (-y(j) - 2.*y(j).*time) / (1 - 2.*time.^2);
%     end
%     for j=1:N
%         v(N,j) = (1 - y(j) - 2.*y(j).*time) / (1 - 2.*time.^2);
%     end
%     for i=1:N
%         v(i,1) = x(i) / (1 - 2.*time^2);
%     end
%     for i=1:N
%         v(i,N) = (x(i) - 1 - 2.*time) / (1 - 2.*time.^2);
%     end
% end
% 
% u_new = u;
% v_new = v;
% 
% for time=0:dt:1
%     for i=2:N-1
%         for j=2:N-1
%             u_new(i,j) = u(i,j) - (1/(2*del_x))*(dt*((u(i+1,j) - u(i-1,j)))*u(i,j)) - (1/(2*del_y))*(dt*((u(i,j+1)-u(i,j-1)))*v(i,j))+(vis*dt*(u(i+1,j)-2*u(i,j)+u(i-1,j))/(del_x^2))+(vis*dt*(u(i,j-1)-2*u(i,j)+u(i,j+1))/(del_y^2));
%             v_new(i,j) = v(i,j) - (1/(2*del_x))*(dt*((v(i+1,j) - v(i-1,j)))*u(i,j)) - (1/(2*del_y))*(dt*((v(i,j+1)-v(i,j-1)))*v(i,j))+(vis*dt*(v(i+1,j)-2*v(i,j)+v(i-1,j))/(del_x^2))+(vis*dt*(v(i,j-1)-2*v(i,j)+v(i,j+1))/(del_y^2));
%         end
%         for i=2:N-1
%             for j=2:N-1
%                 u(i,j)=u_new(i,j);
%                 v(i,j)=v_new(i,j);
%             end
%         end
%     end
% end
%                 %u_new(i,j)=u(i,j)-(dt*(u(i,j)-u(i-1,j)).*u(i,j)/del_x)-(dt*(u(i,j)-u(i,j-1)).*v(i,j)/del_y)+(vis*dt*(u(i+1,j)-2*u(i,j)+u(i-1,j))/(del_x^2))+(vis*dt*(u(i,j-1)-2*u(i,j)+u(i,j+1))/(del_y^2));
%                 %v_new(i,j)=v(i,j)-(dt*(v(i,j)-v(i-1,j)).*u(i,j)/del_x)-(dt*(v(i,j)-v(i,j-1)).*v(i,j)/del_y)+(vis*dt*(v(i+1,j)-2*v(i,j)+v(i-1,j))/(del_x^2))+(vis*dt*(v(i,j-1)-2*v(i,j)+v(i,j+1))/(del_y^2));
% %                 if u(i,j) > 0
% %                     du_del_x = (u(i,j)-u(i-1,j))/del_x;
% %                 else
% %                     du_del_x = (u(i-1,j) - u(i,j))/del_x;
% %                 end
% %                 if v(i,j) > 0
% %                     dv_del_x = (v(i,j)-v(i-1,j))/del_x;
% %                 else
% %                     dv_del_x = (v(i-1,j) - v(i,j))/del_x;
% %                 end
% %                 if u(i,j) > 0
% %                     du_del_y = (u(i,j)-u(i,j-1))/del_y;
% %                 else
% %                     du_del_y = (u(i,j-1) - u(i,j))/del_y;
% %                 end
% %                 if v(i,j) > 0
% %                     dv_del_y = (v(i,j)-v(i,j-1))/del_y;
% %                 else
% %                     dv_del_y = (v(i,j-1) - v(i,j))/del_y;
% %                 end
%                 %u_new(i,j)=u(i,j)-(dt*du_del_x).*u(i,j)-(dt*du_del_y).*v(i,j)+(vis*dt*(u(i+1,j)-2*u(i,j)+u(i-1,j))/(del_x^2))+(vis*dt*(u(i,j-1)-2*u(i,j)+u(i,j+1))/(del_y^2));
%                 %v_new(i,j)=v(i,j)-(dt*dv_del_x).*u(i,j)-(dt*dv_del_y).*v(i,j)+(vis*dt*(v(i+1,j)-2*v(i,j)+v(i-1,j))/(del_x^2))+(vis*dt*(v(i,j-1)-2*v(i,j)+v(i,j+1))/(del_y^2));
%                 %-(dt*du_del_y).*v(i,j)
%                 %-(dt*dv_del_y).*v(i,j)
% %         end
% %     end
% %     error_x = u_new - u;
% %     error_y = v_new - v;
% %     total_error = 0;
% %     for i=1:N
% %         for j=1:N
% %             total_error = total_error + abs(error_x(i,j))+abs(error_y(i,j));
% %         end
% %     end
% %     for i=2:N-1
% %         for j=2:N-1
% %             u(i,j)=u_new(i,j);
% %             v(i,j)=v_new(i,j);
% %         end
% %     end
% % end
% 
% 
% %EXACT SOLUTION
% 
% for time=0:dt:1
%     for i = 1:N
%         for j = 1:N
%             u_exact(i,j) = (x(i)+y(j)-(2.*x(i).*time)) ./ (1-(2.*time.^2));
%             v_exact(i,j) = (x(i)-y(j)-(2.*y(j).*time)) ./ (1-(2.*time.^2));
%         end
%     end
% end
% 
% 
% for it=0:1000
%     %plotting the velocity field
%     h = plot(v_exact(50,:), 'LineWidth', 2)
%     axis([0 1 0 1])
%     axis square
%     title({['2-D Burgers'' equation with {\nu} = ',num2str(vis)];'Transport property vector field {\bfu}=(u_x,u_y)';['time(\itt) = ',num2str(dt*it)]})
%     xlabel('Spatial co-ordinate (x) \rightarrow')
%     ylabel('Spatial co-ordinate (y) \rightarrow')
%     drawnow;
%     refreshdata(h)
% end


% figure(1);
% plot(v_exact(50,:), 'LineWidth', 2)
% hold on
% plot(v(50,:), 'LineWidth', 2)
% legend('Exact', 'Explicit')
% title('u - Explicit Scheme vs Exact Value (Time-accurate manner)')
% hold off

% plot(u(25,:))
% hold on
% plot(u_exact(25,:))
% hold off