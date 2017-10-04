function [Hx,Hy,Hz] = HHP3D(u,v,w,dx1,dx2,dy1,dy2,dz1,dz2,nx,ny,nz,Re,g)

% X-Momentum
ue2 = (1/2.*(u(4:nx,2:ny-1,2:nz-1) + u(3:nx-1,2:ny-1,2:nz-1))).^2;
uw2 = (1/2.*(u(3:nx-1,2:ny-1,2:nz-1) + u(2:nx-2,2:ny-1,2:nz-1))).^2;
du2dx = (ue2-uw2)./dx2(2:nx-2,2:ny-1,2:nz-1);

unu = 1/2.*(u(3:nx-1,3:ny,2:nz-1) + u(3:nx-1,2:ny-1,2:nz-1));
vnu = 1/2.*(v(3:nx-1,3:ny,2:nz-1) + v(2:nx-2,3:ny,2:nz-1));
usu = 1/2.*(u(3:nx-1,2:ny-1,2:nz-1) + u(3:nx-1,1:ny-2,2:nz-1));
vsu = 1/2.*(v(3:nx-1,2:ny-1,2:nz-1) + v(2:nx-2,2:ny-1,2:nz-1));
duvdy = (unu.*vnu - usu.*vsu)./dy1(2:nx-2,2:ny-1,2:nz-1);

utu = 1/2.*(u(3:nx-1,2:ny-1,3:nz) + u(3:nx-1,2:ny-1,2:nz-1));
wtu = 1/2.*(w(2:nx-2,2:ny-1,3:nz) + w(3:nx-1,2:ny-1,3:nz));
ubu = 1/2.*(u(3:nx-1,2:ny-1,1:nz-2) + u(3:nx-1,2:ny-1,2:nz-1));
wbu = 1/2.*(w(2:nx-2,2:ny-1,2:nz-1) + w(3:nx-1,2:ny-1,2:nz-1));
duwdz = (utu.*wtu - ubu.*wbu)./dz1(2:nx-2,2:ny-1,2:nz-1);

d2udx2 = ((u(4:nx,2:ny-1,2:nz-1) - u(3:nx-1,2:ny-1,2:nz-1))./dx1(3:nx-1,2:ny-1,2:nz-1) ...
    - (u(3:nx-1,2:ny-1,2:nz-1) - u(2:nx-2,2:ny-1,2:nz-1))./dx1(2:nx-2,2:ny-1,2:nz-1)) ./ dx2(2:end-1,2:ny-1,2:nz-1);
d2udy2 = ((u(3:nx-1,3:ny,2:nz-1) - u(3:nx-1,2:ny-1,2:nz-1))./dy2(2:nx-2,2:ny-1,2:nz-1) ...
    - (u(3:nx-1,2:ny-1,2:nz-1) - u(3:nx-1,1:ny-2,2:nz-1))./dy2(2:nx-2,1:ny-2,2:nz-1)) ./ dy1(2:nx-2,2:ny-1,2:nz-1);
d2udz2 = ((u(3:nx-1,2:ny-1,3:nz) - u(3:nx-1,2:ny-1,2:nz-1))./dz2(2:nx-2,2:ny-1,2:nz-1) ...
    - (u(3:nx-1,2:ny-1,2:nz-1) - u(3:nx-1,2:ny-1,1:nz-2))./dz2(2:nx-2,2:ny-1,1:nz-2)) ./ dz1(2:nx-2,2:ny-1,2:nz-1);

Hx(3:nx-1,2:ny-1,2:nz-1) = -du2dx - duvdy - duwdz + (1/Re).*(d2udx2 + d2udy2 + d2udz2);

% Y-Momentum
vn2 = (1/2.*(v(2:nx-1,4:ny,2:nz-1) + v(2:nx-1,3:ny-1,2:nz-1))).^2;
vs2 = (1/2.*(v(2:nx-1,3:ny-1,2:nz-1) + v(2:nx-1,2:ny-2,2:nz-1))).^2;
dv2dy = (vn2-vs2)./dy2(2:nx-1,2:ny-2,2:nz-1);

uev = 1/2.*(u(3:nx,3:ny-1,2:nz-1) + u(3:nx,2:ny-2,2:nz-1));
vev = 1/2.*(v(2:nx-1,3:ny-1,2:nz-1) + v(3:nx,3:ny-1,2:nz-1));
uwv = 1/2.*(u(2:nx-1,3:ny-1,2:nz-1) + u(2:nx-1,2:ny-2,2:nz-1));
vwv = 1/2.*(v(1:nx-2,3:ny-1,2:nz-1) + v(2:nx-1,3:ny-1,2:nz-1));
duvdx = (uev.*vev - uwv.*vwv)./dx1(2:nx-1,3:ny-1,2:nz-1);

vtv = 1/2.*(v(2:nx-1,3:ny-1,2:nz-1) + v(2:nx-1,3:ny-1,3:nz));
wtv = 1/2.*(w(2:nx-1,2:ny-2,3:nz) + w(2:nx-1,3:ny-1,3:nz));
vbv = 1/2.*(v(2:nx-1,3:ny-1,1:nz-2) + v(2:nx-1,3:ny-1,2:nz-1));
wbv = 1/2.*(w(2:nx-1,2:ny-2,2:nz-1) + w(2:nx-1,3:ny-1,2:nz-1));
dvwdz = (vtv.*wtv - vbv.*wbv)./dz1(2:nx-1,2:ny-2,2:nz-1);

d2vdx2 = ((v(3:nx,3:ny-1,2:nz-1) - v(2:nx-1,3:ny-1,2:nz-1))./dx2(2:nx-1,2:ny-2,2:nz-1) ...
    - (v(2:nx-1,3:ny-1,2:nz-1) - v(1:nx-2,3:ny-1,2:nz-1))./dx2(1:nx-2,2:ny-2,2:nz-1))./dx1(2:nx-1,2:ny-2,2:nz-1);
d2vdy2 = ((v(2:nx-1,4:ny,2:nz-1) - v(2:nx-1,3:ny-1,2:nz-1))./dy1(2:nx-1,3:ny-1,2:nz-1) ...
    - (v(2:nx-1,3:ny-1,2:nz-1) - v(2:nx-1,2:ny-2,2:nz-1))./dy1(2:nx-1,2:ny-2,2:nz-1))./dy2(2:nx-1,2:end-1,2:nz-1);
d2vdz2 = ((v(2:nx-1,3:ny-1,3:nz) - v(2:nx-1,3:ny-1,2:nz-1))./dz2(2:nx-1,2:ny-2,2:nz-1) ...
    - (v(2:nx-1,3:ny-1,2:nz-1) - v(2:nx-1,3:ny-1,1:nz-2))./dz2(2:nx-1,2:ny-2,1:nz-2))./dz1(2:nx-1,2:ny-2,2:nz-1);

Hy(2:nx-1,3:ny-1,2:nz-1) = -dv2dy - duvdx - dvwdz + (1/Re).*(d2vdx2 + d2vdy2 +d2vdz2);

% Z-Momentum
wt2 = (1/2.*(w(2:nx-1,2:ny-1,4:nz) + w(2:nx-1,2:ny-1,3:nz-1))).^2;
wb2 = (1/2.*(w(2:nx-1,2:ny-1,3:nz-1) + w(2:nx-1,2:ny-1,2:nz-2))).^2;
dw2dz = (wt2 - wb2)./dz2(2:nx-1,2:ny-1,2:nz-2);

uew = 1/2.*(u(3:nx,2:ny-1,3:nz-1) + u(3:nx,2:ny-1,2:nz-2));
wew = 1/2.*(w(2:nx-1,2:ny-1,3:nz-1) + w(3:nx,2:ny-1,3:nz-1));
uww = 1/2.*(u(2:nx-1,2:ny-1,3:nz-1) + u(2:nx-1,2:ny-1,2:nz-2));
www = 1/2.*(w(1:nx-2,2:ny-1,3:nz-1) + w(2:nx-1,2:ny-1,3:nz-1));
duwdx = (uew.*wew - uww.*www)./dx1(2:nx-1,2:ny-1,3:nz-1);

vnw = 1/2.*(v(2:nx-1,3:ny,3:nz-1) + v(2:nx-1,3:ny,2:nz-2));
wnw = 1/2.*(w(2:nx-1,2:ny-1,3:nz-1) + w(2:nx-1,3:ny,3:nz-1));
vsw = 1/2.*(v(2:nx-1,2:ny-1,3:nz-1) + v(2:nx-1,2:ny-1,2:nz-2));
wsw = 1/2.*(w(2:nx-1,1:ny-2,3:nz-1) + w(2:nx-1,2:ny-1,3:nz-1));
dvwdy = (vnw.*wnw - vsw.*wsw)./dy1(2:nx-1,2:ny-1,3:nz-1);

d2wdx2 = ((w(3:nx,2:ny-1,3:nz-1) - w(2:nx-1,2:ny-1,3:nz-1))./dx2(2:nx-1,2:ny-1,3:nz-1) ...
    - (w(2:nx-1,2:ny-1,3:nz-1) - w(1:nx-2,2:ny-1,3:nz-1))./dx2(1:nx-2,2:ny-1,3:nz-1))./dx1(2:nx-1,2:ny-1,3:nz-1);
d2wdy2 = ((w(2:nx-1,3:ny,3:nz-1) - w(2:nx-1,2:ny-1,3:nz-1))./dy2(2:nx-1,2:ny-1,3:nz-1) ...
    - (w(2:nx-1,2:ny-1,3:nz-1) - w(2:nx-1,1:ny-2,3:nz-1))./dy2(2:nx-1,1:ny-2,3:nz-1))./dy1(2:nx-1,2:ny-1,3:nz-1);
d2wdz2 = ((w(2:nx-1,2:ny-1,4:nz) - w(2:nx-1,2:ny-1,3:nz-1))./dz1(2:nx-1,2:ny-1,3:nz-1) ...
    - (w(2:nx-1,2:ny-1,3:nz-1) - w(2:nx-1,2:ny-1,2:nz-2))./dz1(2:nx-1,2:ny-1,2:nz-2))./dz2(2:nx-1,2:ny-1,2:nz-2);

Hz(2:nx-1,2:ny-1,3:nz-1) = -dw2dz - duwdx - dvwdy + (1/Re).*(d2wdx2 + d2wdy2 + d2wdz2) + g;