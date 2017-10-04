function [] = vectorxy(mx,my,mz,uc,vc,L,tstring,pdir)

cd(pdir);
mkdir('VectorXY');
pdir = [pwd,'\','VectorXY'];
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
        
        figure(700+i)

        quiver(MX,MY,u',v')
        axis([0,L(1),0,L(2)])
        xlabel('X')
        ylabel('Y')
        title(['Vektor Kecepatan Bidang XY, Z = ',num2str(mz(i)),' ',tstring])

        fname = ['VectorXY',num2str(i)];
        
        
        saveas(700+i,fname,'jpeg');
        saveas(700+i,fname,'fig');
        
    end
else
    k = str2double(k);
    u = squeeze(uc(:,:,k));
    v = squeeze(vc(:,:,k));
    
    figure(700+k)

    quiver(MX,MY,u',v')
    axis([0,L(1),0,L(2)])
    xlabel('X')
    ylabel('Y')
    title(['Vektor Kecepatan Bidang XY, Z = ',num2str(mz(k)),' ',tstring])

    fname = ['VectorXY',num2str(k)];
        
        
    saveas(700+k,fname,'jpeg');
    saveas(700+k,fname,'fig');
        
end
cd ..
cd ..