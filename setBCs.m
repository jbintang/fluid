function [u,v] = setBCs(u,v,V,ny,vel,arahv,fil)

%% Sumber Kecepatan dan Dinding
if strcmp(vel,'SD') == 1
    if strcmp(V,'N') == 1
        % North
        u(:,end) = 2*arahv - u(:,end-1);
        v(:,end) = 0;
        % South
        u(:,1) = - u(:,2);
        v(:,2) = 0;
        % West
        v(1,:) = - v(2,:);
        u(2,:) = 0;
        % East
        v(end,:) = -v(end-1,:);
        u(end,:) = 0;
    elseif strcmp(V,'S') == 1
        % North
        u(:,end) = - u(:,end-1);
        v(:,end) = 0;
        % South
        u(:,1) = 2*arahv - u(:,2);
        v(:,2) = 0;
        % West
        v(1,:) = - v(2,:);
        u(2,:) = 0;
        % East
        v(end,:) = -v(end-1,:);
        u(end,:) = 0;
    elseif strcmp(V,'E') == 1
        % North
        u(:,end) = - u(:,end-1);
        v(:,end) = 0;
        % South
        u(:,1) = - u(:,2);
        v(:,2) = 0;
        % West
        v(1,:) = - v(2,:);
        u(2,:) = 0;
        % East
        v(end,:) = 2*arahv -v(end-1,:);
        u(end,:) = 0;
    elseif strcmp(V,'W') == 1
        % North
        u(:,end) = - u(:,end-1);
        v(:,end) = 0;
        % South
        u(:,1) = - u(:,2);
        v(:,2) = 0;
        % West
        v(1,:) = 2*arahv - v(2,:);
        u(2,:) = 0;
        % East
        v(end,:) = -v(end-1,:);
        u(end,:) = 0;
    end
end

%% Fillet
if strcmp(fil,'Y') == 1
    f = (ny-2)/10;
    for i = 1:f
        for j = 1:f
            if f >= 2 
                if j~=f && i~=f
                    % SW
                    u(2+i,1+j) = 0;
                    v(1+i,2+j) = 0;
                    % SE
                    u(end-i,1+j) = 0;
                    v(end-i,2+j) = 0;
                    % NW
                    u(2+i,end-j) = 0;
                    v(1+i,end-j) = 0;
                    % SW
                    u(end-i,end-j) = 0;
                    v(end-i,end-j) = 0;
                end
            else
                 % SW
                 u(2+i,1+j) = 0;
                 v(1+i,2+j) = 0;
                 % SE
                 u(end-i,1+j) = 0;
                 v(end-i,2+j) = 0;
                 % NW
                 u(2+i,end-j) = 0;
                 v(1+i,end-j) = 0;
                 % SW
                 u(end-i,end-j) = 0;
                 v(end-i,end-j) = 0;
            end
        end
    end
end