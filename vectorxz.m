function [] = vectorxz(mx,my,mz,uc,wc,L,tstring,pdir)

cd(pdir);
mkdir('VectorXZ');
pdir = [pwd,'\','VectorXZ'];
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
        
        figure(800+i)

        quiver(MX,MZ,u',w')
        axis([0,L(1),0,L(3)])
        xlabel('X')
        ylabel('Z')
        title(['Vektor Kecepatan Bidang XZ, Y = ',num2str(my(i)),' ',tstring])

        fname = ['VectorXZ',num2str(i)];
        
        
        saveas(800+i,fname,'jpeg');
        saveas(800+i,fname,'fig');
    end
else
    k = str2double(k);
    u = squeeze(uc(:,k,:));
    w = squeeze(wc(:,k,:));
    
    figure(800+k)

    quiver(MX,MZ,u',w')
    axis([0,L(1),0,L(3)])
    xlabel('X')
    ylabel('Z')
    title(['Vektor Kecepatan Bidang XZ, Y = ',num2str(my(k)),' ',tstring])

    fname = ['VectorXZ',num2str(k)];
        
    saveas(800+k,fname,'jpeg');
    saveas(800+k,fname,'fig');
    
end
cd ..
cd ..