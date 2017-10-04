function [u,v,w] = setBCs3D(u,v,w,V,vel,Vv,arahv,fil) 

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
end


%% Fillet    
if strcmp(fil,'Y') == 1
    f = (ny-2)/10;
    for i=1:f
        for j=1:f
            for k=1:f
                if f >= 2 
                    if (i==f && j==f) || (i==f && k==f) || (j==f && k==f) 
                        % XY
                        w(1+i,1+j,:) = 0;
                        w(end-i,1+j,:) = 0;
                        w(1+i,endj1,:) = 0;
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
                    w(1+i,endj1,:) = 0;
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
