clear all
%FINAL PROJECT
%QUESTION 2 - ADI for Diffusion term and Explicit for convection term

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

u(1,:)=0;
u(nx,:)=0;
u(:,1)=sin(2*pi.*x);
u(:,ny)=sin(2*pi.*x);
v(:,1)=0;
v(:,ny)=1;
v(1,:)=1-y;
v(nx,:)=1-y;

i = 1:nx;
j = 1:ny;

un(i,j) = u(i,j);
vn(i,j) = v(i,j);

i = 2:nx-1;
j = 2:ny-1;

for it=0:nt
    un=u;
    vn=v;
    u(i,j)=un(i,j)-(dt*(un(i,j)-un(i-1,j)).*un(i,j)/dx)-(dt*(un(i,j)-un(i,j-1)).*vn(i,j)/dy)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j-1)-2*un(i,j)+un(i,j+1))/(dy*dy));
    v(i,j)=vn(i,j)-(dt*(vn(i,j)-vn(i-1,j)).*un(i,j)/dx)-(dt*(vn(i,j)-vn(i,j-1)).*vn(i,j)/dy)+(vis*dt*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))/(dx*dx))+(vis*dt*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/(dy*dy));
end

for i=2:nx-1
    for j=2:ny-1
        A_xU(i,j) = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/dx^2;
        A_yU(i,j) = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dy^2;
        A_xV(i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dx^2;
        A_yV(i,j) = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dy^2;
    end
end
% 
A_xU(:,ny)=0;
A_xU(nx,:)=0;
A_xV(:,ny)=0;
A_xV(nx,:)=0;
A_yU(:,ny)=0;
A_yU(nx,:)=0;
A_yV(:,ny)=0;
A_yV(nx,:)=0;

II = eye(nx,ny);  %Identity matrix

term_x_u = (II - (vis/2)*dt*A_xU);
term_y_u = (II + (vis/2)*dt*A_yU);
term_x_v = (II - (vis/2)*dt*A_xV);
term_y_v = (II + (vis/2)*dt*A_yV);

i=2:nx-1;
j=2:ny-1;

term_x_u = inv(term_x_u);
term_x_v = inv(term_x_v);

% %u^n+1/2 calculation
for iterations=1:1000
    
    u_new(i,j)=u(i,j)-(dt*(u(i,j)-u(i-1,j)).*u(i,j)/dx)-(dt*(u(i,j)-u(i,j-1)).*v(i,j)/dy) + (term_x_u(i,j).*term_y_u(i,j).*un(i,j));
    v_new(i,j)=v(i,j)-(dt*(v(i,j)-v(i-1,j)).*u(i,j)/dx)-(dt*(v(i,j)-v(i,j-1)).*v(i,j)/dy) + (term_x_v(i,j).*term_y_v(i,j).*vn(i,j));
end
% 
% u_new = zeros(nx,ny);
% v_new = zeros(nx,ny);
% 
term_x_uu = (II + (vis/2)*dt*A_xU);
term_y_uu = (II - (vis/2)*dt*A_yU);
term_x_vv = (II + (vis/2)*dt*A_xV);
term_y_vv = (II - (vis/2)*dt*A_yV);
% 
% 
term_y_uu = inv(term_x_u);
term_y_vv = inv(term_x_v);
% 


% %u^n+1 calculation
for iterations=1:1000
    u_neww(i,j)=u(i,j)-(dt*(u(i,j)-u(i-1,j)).*u(i,j)/dx)-(dt*(u(i,j)-u(i,j-1)).*v(i,j)/dy) + (term_y_uu(i,j).*term_x_uu(i,j).*u_new(i,j));
    v_neww(i,j)=v(i,j)-(dt*(v(i,j)-v(i-1,j)).*u(i,j)/dx)-(dt*(v(i,j)-v(i,j-1)).*v(i,j)/dy) + (term_y_vv(i,j).*term_x_vv(i,j).*v_new(i,j));
end

u_neww(:,ny)=0;
u_neww(nx,:)=0;
v_neww(:,ny)=0;
v_neww(nx,:)=0;

figure(1);
quiver(x,y,u_neww, v_neww, 'k')
title('Surface plot of u and v via ADI method')
xlabel('x coordinate \rightarrow')
ylabel('y coordinate \rightarrow')

figure(2);
plot(u_neww(25,:), 'LineWidth', 2)
hold on
plot(v_neww(25,:), 'LineWidth', 2)
title('u vs v at the middle slice (25,:)')
legend('u explicit', 'v explicit')
hold off

%QUESTION 2

% x = linspace(0,1,50);
% y = linspace(0,1,50);
% dt = 0.005;
% time = 0:dt:0.5;
% N = 50;                       %Number of steps in x, y
% vis = 0.015;
% del_x = x(2)-x(1);            %Width of space step(x)
% del_y = y(2)-y(1);            %Width of space step(y)
% 
% u=zeros(N,N);                 %preallocation
% v=zeros(N,N);                 %preallocation
% 
% %ADI METHOD
% 
% %PREPARING THE Ax and Ay MATRIX FOR U and V
% 

% 
% A_xU(:,50)=0;
% A_xU(50,:)=0;
% A_xV(:,50)=0;
% A_xV(50,:)=0;
% A_yU(:,50)=0;
% A_yU(50,:)=0;
% A_yV(:,50)=0;
% A_yV(50,:)=0;
% 
% II = eye(50,50);
% 
% %PREPARING THE TERMS SEPARATELY
% 

% 
% %INVERT THE term_x MATRIX FOR U and V  [n+1/2]
% 
% term_x_u = inv(term_x_u);
% term_x_v = inv(term_x_v);
% 
% 
% for i=1:50
%     for j=1:50
%         u_one_half(i,j) = term_x_u(i,j) * term_y_u(i,j) * u(i,j);
%         v_one_half(i,j) = term_x_v(i,j) * term_y_v(i,j) * v(i,j);
%     end
% end
% 
% for i=1:50
%     for j=1:50
%         u_n_plus_one(i,j) = inv(term_x_u(i,j)) * inv(term_y_u(i,j)) * u_one_half(i,j);
%         v_n_plus_one(i,j) = inv(term_x_v(i,j)) * inv(term_y_v(i,j)) * v_one_half(i,j);
%     end
% end
% 
% %MANDATE THE BOUNDARY CONDITIONS AGAIN
% u_n_plus_one(1,:)=0;
% u_n_plus_one(50,:)=0;
% u_n_plus_one(:,1)=sin(2*pi.*x);
% u_n_plus_one(:,50)=sin(2*pi.*x);
% v_n_plus_one(:,1)=0;
% v_n_plus_one(:,50)=1;
% v_n_plus_one(1,:)=1-y;
% v_n_plus_one(50,:)=1-y;