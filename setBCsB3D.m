function [u,v,w] = setBCsB3D(u,v,w,V,ibw,ibb,ibt,nx,ny,nz,vel,Vv,arahv,ifan,fil) 

%% Baffle
for i = ibb : (nz-2-ibt)
    if i == nz-2-ibt
        w(2+ibw,:,2+i) = 0;
        v(2+ibw,:,2+i) = 0;
    else
        u(2+ibw:3+ibw,:,2+i) = 0;
        v(2+ibw,:,2+i) = 0;
        w(2+ibw,:,2+i) = 0;
    end
end

%% Sumber Kecepatan dan Dinding
if strcmp(vel,'SD') == 1
    if strcmp(V,'T') == 1
        if strcmp(Vv,'X') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = 2*arahv - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Y') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = 2*arahv - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        end
    elseif strcmp(V,'B') == 1
        if strcmp(Vv,'X') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = 2*arahv - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Y') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = 2*arahv - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        end
    elseif strcmp(V,'N') == 1
        if strcmp(Vv,'X') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = 2*arahv - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Z') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = 2*arahv - w(:,end-1,:);
        end
    elseif strcmp(V,'S') == 1
        if strcmp(Vv,'X') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = 2*arahv - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Z') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = 2*arahv - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        end
    elseif strcmp(V,'E') == 1
        if strcmp(Vv,'Y') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = 2*arahv -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Z') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = 2*arahv -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        end
    elseif strcmp(V,'W') == 1
        if strcmp(Vv,'Y') == 1
            % West
            w(1,:,:) = - w(2,:,:);
            v(1,:,:) = 2*arahv - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        elseif strcmp(Vv,'Z') == 1
            % West
            w(1,:,:) = 2*arahv - w(2,:,:);
            v(1,:,:) = - v(2,:,:);
            u(2,:,:) = 0;

            % East
            w(end,:,:) = -w(end-1,:,:);
            v(end,:,:) = -v(end-1,:,:);
            u(end,:,:) = 0;

            % Bottom
            u(:,:,1) = - u(:,:,2);
            v(:,:,1) = - v(:,:,2);
            w(:,:,2) = 0;

            % Top
            u(:,:,end) = - u(:,:,end-1);
            v(:,:,end) = - v(:,:,end-1);
            w(:,:,end) = 0;

            % South
            u(:,1,:) = - u(:,2,:);
            v(:,2,:) = 0;
            w(:,1,:) = - w(:,2,:);

            % North
            u(:,end,:) = - u(:,end-1,:);
            v(:,end,:) = 0;
            w(:,end,:) = - w(:,end-1,:);
        end
    end
elseif strcmp(vel,'F') == 1
    for i=1:ibw
        w(1+i,:,2+ifan) = arahv;
    end
    % West
    w(1,:,:) = - w(2,:,:);
    v(1,:,:) = - v(2,:,:);
    u(2,:,:) = 0;

    % East
    w(end,:,:) = -w(end-1,:,:);
    v(end,:,:) = -v(end-1,:,:);
    u(end,:,:) = 0;

    % Bottom
    u(:,:,1) = - u(:,:,2);
    v(:,:,1) = - v(:,:,2);
    w(:,:,2) = 0;

    % Top
    u(:,:,end) = - u(:,:,end-1);
    v(:,:,end) = - v(:,:,end-1);
    w(:,:,end) = 0;

    % South
    u(:,1,:) = - u(:,2,:);
    v(:,2,:) = 0;
    w(:,1,:) = - w(:,2,:);

    % North
    u(:,end,:) = - u(:,end-1,:);
    v(:,end,:) = 0;
    w(:,end,:) = - w(:,end-1,:);
end

%% Fillet
if strcmp(fil,'Y') == 1
    fy = (ny-2)/10;
    fx = (nx-3)/10;
    fz = (nz-2)/10;
    for i=1:fx
        for j=1:fy
            for k=1:fz
                if fx >= 2 
                    if (i==fx && j==fy) || (i==fx && k==fz) || (j==fy && k==fz)
                    else
                        % XY
                        w(1+i,1+j,:) = 0;
                        w(end-i,1+j,:) = 0;
                        w(1+i,end-1,:) = 0;
                        w(end-i,end-j,:) = 0;                    

                        v(1+i,2+j,:) = 0;
                        v(end-i,2+j,:) = 0;                    
                        v(1+i,end-j,:) = 0;                    
                        v(end-i,end-j,:) = 0;

                        u(2+i,1+j,:) = 0;                     
                        u(end-i,1+j,:) = 0;                       
                        u(2+i,end-j,:) = 0;                      
                        u(end-i,end-j,:) = 0;  

                        % XZ                      
                        w(i+1,:,2+k) = 0;                     
                        w(end-i,:,2+k) = 0;                       
                        w(i+1,:,end-k) = 0;                        
                        w(end-i,:,end-k) = 0;   

                        v(1+i,:,1+k) = 0;                       
                        v(end-i,:,1+k) = 0;                         
                        v(1+i,:,end-k) = 0;                       
                        v(end-i,:,end-k) = 0;   

                        u(2+k,:,1+k) = 0;                        
                        u(end-i,:,1+k) = 0;                        
                        u(2+k,:,end-k) = 0;                       
                        u(end-i,:,end-k) = 0;    

                        % YZ                       
                        w(:,1+j,2+k) = 0;                       
                        w(:,end-j,2+k) = 0;                        
                        w(:,1+j,end-k) = 0;                        
                        w(:,end-j,end-k) = 0;   

                        v(:,2+j,1+k) = 0;                
                        v(:,end-j,1+k) = 0;                      
                        v(:,2+j,end-k) = 0;                       
                        v(:,end-j,end-k) = 0;   

                        u(:,1+j,1+k) = 0;                      
                        u(:,end-j,1+k) = 0;                      
                        u(:,1+j,end-k) = 0;                    
                        u(:,end-j,end-k) = 0;                    
                    end
                else
                    % XY
                    w(1+i,1+j,:) = 0;
                    w(end-i,1+j,:) = 0;
                    w(1+i,end-1,:) = 0;
                    w(end-i,end-j,:) = 0;

                    v(1+i,2+j,:) = 0;
                    v(end-i,2+j,:) = 0;
                    v(1+i,end-j,:) = 0;
                    v(end-i,end-j,:) = 0;

                    u(2+i,1+j,:) = 0; 
                    u(end-i,1+j,:) = 0;   
                    u(2+i,end-j,:) = 0;  
                    u(end-i,end-j,:) = 0;

                    % XZ
                    w(i+1,:,2+k) = 0;
                    w(end-i,:,2+k) = 0;
                    w(i+1,:,end-k) = 0;
                    w(end-i,:,end-k) = 0;

                    v(1+i,:,1+k) = 0;
                    v(end-i,:,1+k) = 0;
                    v(1+i,:,end-k) = 0;
                    v(end-i,:,end-k) = 0;

                    u(2+k,:,1+k) = 0;
                    u(end-i,:,1+k) = 0;
                    u(2+k,:,end-k) = 0; 
                    u(end-i,:,end-k) = 0;

                    % YZ
                    w(:,1+j,2+k) = 0;
                    w(:,end-j,2+k) = 0;
                    w(:,1+j,end-k) = 0;
                    w(:,end-j,end-k) = 0;
   
                    v(:,2+j,1+k) = 0;
                    v(:,end-j,1+k) = 0;
                    v(:,2+j,end-k) = 0;
                    v(:,end-j,end-k) = 0;

                    u(:,1+j,1+k) = 0;  
                    u(:,end-j,1+k) = 0;   
                    u(:,1+j,end-k) = 0; 
                    u(:,end-j,end-k) = 0;
                end
            end
        end
    end
end