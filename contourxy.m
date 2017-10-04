function [] = contourxy(mx,my,mz,uc,vc,L,tstring,pdir)

cd(pdir);
mkdir('ContourXY');
pdir = [pwd,'\','ContourXY'];
cd(pdir);

[MX,MY] = meshgrid(mx,my);
[~,j] = size(mz);
prompt = 'Pilih dari ';
for i=2:j-1
    Z = ['(',num2str(i),') ',num2str(mz(i))];
    prompt = [prompt char(10) Z];
end
k = input([prompt,char(10),'all',char(10)],'s');
if strcmp(k,'all') == 1
    for i=2:j-1
        u = squeeze(uc(:,:,i));
        v = squeeze(vc(:,:,i));
        uu = sqrt(u.^2 + v.^2);
        figure(100+i)

        contourf(MX,MY,uu',20)
        colorbar
        axis([0,L(1),0,L(2)])
        xlabel('X')
        ylabel('Y')
        title(['Kontur Bidang XY, Z = ',num2str(mz(i)),' ',tstring])

        fname = ['ContourXY',num2str(i)];
        
        
        saveas(100+i,fname,'jpeg');
        saveas(100+i,fname,'fig');
        
    end
else
    k = str2double(k);
    u = squeeze(uc(:,:,k));
    v = squeeze(vc(:,:,k));
    uu = sqrt(u.^2 + v.^2);
    figure(100+k)

    contourf(MX,MY,uu',20)
    colorbar
    axis([0,L(1),0,L(2)])
    xlabel('X')
    ylabel('Y')
    title(['Kontur Bidang XY, Z = ',num2str(mz(k)),' ',tstring])

    fname = ['ContourXY',num2str(k)];
        
        
    saveas(100+k,fname,'jpeg');
    saveas(100+k,fname,'fig');
        
end
cd ..
cd ..