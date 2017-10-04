function [] = contouryz(mx,my,mz,vc,wc,L,tstring,pdir)

cd(pdir);
mkdir('ContourYZ');
pdir = [pwd,'\','ContourYZ'];
cd(pdir);

[MY,MZ] = meshgrid(my,mz);
[~,j] = size(mx);
prompt = 'Pilih dari ';
for i=2:j-1
    Z = ['(',num2str(i),') ',num2str(mx(i))];
    prompt = [prompt char(10) Z];
end
k = input([prompt,char(10),'all',char(10)],'s');
if strcmp(k,'all') == 1
    for i=2:j-1
        v = squeeze(vc(i,:,:));
        w = squeeze(wc(i,:,:));
        uu = sqrt(v.^2 + w.^2);
        figure(300+i)

        contourf(MY,MZ,uu',20)
        colorbar
        axis([0,L(2),0,L(3)])
        xlabel('Y')
        ylabel('Z')
        title(['Kontur Bidang YZ, X = ',num2str(mx(i)),' ',tstring])

        fname = ['ContourYZ',num2str(i)];
        
        saveas(300+i,fname,'jpeg');
        saveas(300+i,fname,'fig');
        
    end
else
    k = str2double(k);
    v = squeeze(vc(k,:,:));
    w = squeeze(wc(k,:,:));
    uu = sqrt(v.^2 + w.^2);
    figure(300+k)

    contourf(MY,MZ,uu',20)
    colorbar
    axis([0,L(2),0,L(3)])
    xlabel('Y')
    ylabel('Z')
    title(['Kontur Bidang YZ, Y = ',num2str(mx(k)),' ',tstring])

    fname = ['ContourYZ',num2str(k)];
        
    saveas(300+k,fname,'jpeg');
    saveas(300+k,fname,'fig');
        
end
cd ..
cd ..