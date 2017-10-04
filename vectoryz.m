function [] = vectoryz(mx,my,mz,vc,wc,L,tstring,pdir)

cd(pdir);
mkdir('VectorYZ');
pdir = [pwd,'\','VectorYZ'];
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
        
        figure(900+i)

        quiver(MY,MZ,v',w')
        axis([0,L(2),0,L(3)])
        xlabel('Y')
        ylabel('Z')
        title(['Vektor Kecepatan Bidang YZ, X = ',num2str(mx(i)),' ',tstring])

        fname = ['VectorYZ',num2str(i)];
        
        saveas(900+i,fname,'jpeg');
        saveas(900+i,fname,'fig');
        
    end
else
    k = str2double(k);
    v = squeeze(vc(k,:,:));
    w = squeeze(wc(k,:,:));
    
    figure(900+k)

    quiver(MY,MZ,v',w')
    axis([0,L(2),0,L(3)])
    xlabel('Y')
    ylabel('Z')
    title(['Vektor Kecepatan Bidang YZ, Y = ',num2str(mx(k)),' ',tstring])

    fname = ['VectorYZ',num2str(k)];
        
    saveas(900+k,fname,'jpeg');
    saveas(900+k,fname,'fig');
        
end
cd ..
cd ..