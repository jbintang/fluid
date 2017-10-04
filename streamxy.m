function [] = streamxy(mx,my,mz,uc,vc,L,tstring,pdir)

cd(pdir);
mkdir('StreamlineXY');
pdir = [pwd,'\','StreamlineXY'];
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
        
        figure(400+i)

        streamslice(MX,MY,u',v',2)
        axis([0,L(1),0,L(2)])
        xlabel('X')
        ylabel('Y')
        title(['Garis Arus Bidang XY, Z = ',num2str(mz(i)),' ',tstring])

        fname = ['StreamlineXY',num2str(i)];
        
        
        saveas(400+i,fname,'jpeg');
        saveas(400+i,fname,'fig');
        
    end
else
    k = str2double(k);
    u = squeeze(uc(:,:,k));
    v = squeeze(vc(:,:,k));
    
    figure(400+k)

    streamslice(MX,MY,u',v',2)
    axis([0,L(1),0,L(2)])
    xlabel('X')
    ylabel('Y')
    title(['Garis Arus Bidang XY, Z = ',num2str(mz(k)),' ',tstring])

    fname = ['StreamlineXY',num2str(k)];
        
        
    saveas(400+k,fname,'jpeg');
    saveas(400+k,fname,'fig');
        
end
cd ..
cd ..