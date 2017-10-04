function [] = contourxz(mx,my,mz,uc,wc,L,tstring,pdir)

cd(pdir);
mkdir('ContourXZ');
pdir = [pwd,'\','ContourXZ'];
cd(pdir);

[MX,MZ] = meshgrid(mx,mz);
[~,j] = size(my);
prompt = 'Pilih dari ';
for i=2:j-1
    Z = ['(',num2str(i),') ',num2str(my(i))];
    prompt = [prompt char(10) Z];
end
k = input([prompt,char(10),'all',char(10)],'s');
if strcmp(k,'all') == 1
    for i=2:j-1
        u = squeeze(uc(:,i,:));
        w = squeeze(wc(:,i,:));
        uu = sqrt(u.^2 + w.^2);
        figure(200+i)

        contourf(MX,MZ,uu',20)
        colorbar
        axis([0,L(1),0,L(3)])
        xlabel('X')
        ylabel('Z')
        title(['Kontur Bidang XZ, Y = ',num2str(my(i)),' ',tstring])

        fname = ['ContourXZ',num2str(i)];
        
        
        saveas(200+i,fname,'jpeg');
        saveas(200+i,fname,'fig');
        
    end
else
    k = str2double(k);
    u = squeeze(uc(:,k,:));
    w = squeeze(wc(:,k,:));
    uu = sqrt(u.^2 + w.^2);
    figure(200+k)

    contourf(MX,MZ,uu',20)
    colorbar
    axis([0,L(1),0,L(3)])
    xlabel('X')
    ylabel('Z')
    title(['Kontur Bidang XZ, Y = ',num2str(my(k)),' ',tstring])

    fname = ['ContourXZ',num2str(k)];
        
    saveas(200+k,fname,'jpeg');
    saveas(200+k,fname,'fig');
        
end
cd ..
cd ..