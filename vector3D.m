function [] = vector3D(mx,my,mz,uc,vc,wc,L,tstring,pdir)
cd(pdir);
[MX,MY,MZ] = meshgrid(mx,my,mz);
figure(2)
quiver3(MX,MY,MZ,ipermute(uc,[2 1 3]),ipermute(vc,[2 1 3]),ipermute(wc,[2 1 3]))
axis([0,L(1),0,L(2),0,L(3)])
xlabel('X')
ylabel('Y')
zlabel('Z')
title(['Vektor Kecepatan 3D, ',tstring])
view(3)
fname = 'VectorXYZ';
        
saveas(2,fname,'jpeg');
saveas(2,fname,'fig');
cd ..