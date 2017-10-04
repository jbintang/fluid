function [] = streamxz(mx,my,mz,uc,wc,L,tstring,pdir)

cd(pdir);
mkdir('StreamlineXZ');
pdir = [pwd,'\','StreamlineXZ'];
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
        
        figure(500+i)

        streamslice(MX,MZ,u',w',2)
        axis([0,L(1),0,L(3)])
        xlabel('X')
        ylabel('Z')
        title(['Garis Arus Bidang XZ, Y = ',num2str(my(i)),' ',tstring])

        fname = ['StreamlineXZ',num2str(i)];
        
        
        saveas(500+i,fname,'jpeg');
        saveas(500+i,fname,'fig');
        
    end
else
    k = str2double(k);
    u = squeeze(uc(:,k,:));
    w = squeeze(wc(:,k,:));
    
    figure(500+k)

    streamslice(MX,MZ,u',w',2)
    axis([0,L(1),0,L(3)])
    xlabel('X')
    ylabel('Z')
    title(['Garis Arus Bidang XZ, Y = ',num2str(my(k)),' ',tstring])

    fname = ['StreamlineXZ',num2str(k)];
        
    saveas(500+k,fname,'jpeg');
    saveas(500+k,fname,'fig');
        
end
cd ..
cd ..