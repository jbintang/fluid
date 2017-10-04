function [Hx,Hy] = HHP(u,v,dx1,dx2,dy1,dy2,nx,ny,Re)


% X-Momentum
ue2 = (1/2.*(u(4:nx,2:ny-1) + u(3:nx-1,2:ny-1))).^2;
uw2 = (1/2.*(u(3:nx-1,2:ny-1) + u(2:nx-2,2:ny-1))).^2;
du2dx = (ue2-uw2)./dx2(2:nx-2,2:ny-1);

unu = 1/2.*(u(3:nx-1,3:ny) + u(3:nx-1,2:ny-1));
vnu = 1/2.*(v(3:nx-1,3:ny) + v(2:nx-2,3:ny));
usu = 1/2.*(u(3:nx-1,2:ny-1) + u(3:nx-1,1:ny-2));
vsu = 1/2.*(v(3:nx-1,2:ny-1) + v(2:nx-2,2:ny-1));
duvdy = (unu.*vnu - usu.*vsu)./dy1(2:nx-2,2:ny-1);

d2udx2 = ((u(4:nx,2:ny-1) - u(3:nx-1,2:ny-1))./dx1(3:nx-1,2:ny-1) ...
    - (u(3:nx-1,2:ny-1) - u(2:nx-2,2:ny-1))./dx1(2:nx-2,2:ny-1)) ./ dx2(2:end-1,2:ny-1);
d2udy2 = ((u(3:nx-1,3:ny) - u(3:nx-1,2:ny-1))./dy2(2:nx-2,2:ny-1) ...
    - (u(3:nx-1,2:ny-1) - u(3:nx-1,1:ny-2))./dy2(2:nx-2,1:ny-2)) ./ dy1(2:nx-2,2:ny-1);

Hx(3:nx-1,2:ny-1) = -du2dx - duvdy + (1/Re).*(d2udx2 + d2udy2);

% Y-Momentum
vn2 = (1/2.*(v(2:nx-1,4:ny) + v(2:nx-1,3:ny-1))).^2;
vs2 = (1/2.*(v(2:nx-1,3:ny-1) + v(2:nx-1,2:ny-2))).^2;
dv2dy = (vn2-vs2)./dy2(2:nx-1,2:ny-2);

uev = 1/2.*(u(3:nx,3:ny-1) + u(3:nx,2:ny-2));
vev = 1/2.*(v(2:nx-1,3:ny-1) + v(3:nx,3:ny-1));
uwv = 1/2.*(u(2:nx-1,3:ny-1) + u(2:nx-1,2:ny-2));
vwv = 1/2.*(v(1:nx-2,3:ny-1) + v(2:nx-1,3:ny-1));
duvdx = (uev.*vev - uwv.*vwv)./dx1(2:nx-1,3:ny-1);

d2vdx2 = ((v(3:nx,3:ny-1) - v(2:nx-1,3:ny-1))./dx2(2:nx-1,2:ny-2) ...
    - (v(2:nx-1,3:ny-1) - v(1:nx-2,3:ny-1))./dx2(1:nx-2,2:ny-2))./dx1(2:nx-1,2:ny-2);
d2vdy2 = ((v(2:nx-1,4:ny) - v(2:nx-1,3:ny-1))./dy1(2:nx-1,3:ny-1) ...
    - (v(2:nx-1,3:ny-1) - v(2:nx-1,2:ny-2))./dy1(2:nx-1,2:ny-2))./dy2(2:nx-1,2:end-1);

Hy(2:nx-1,3:ny-1) = -dv2dy - duvdx + (1/Re).*(d2vdx2 + d2vdy2);
