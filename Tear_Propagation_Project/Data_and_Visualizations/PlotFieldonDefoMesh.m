function PlotFieldonDefoMesh(coordinates,nodes,factor,components,magnitude,clim,cbar_type,cbar_label) 
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on deformed mesh
% Synopsis :
%           ProfileonDefoMesh(coordinates,nodes,component)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y Z] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%           factor - Amplification factor (Change accordingly, trial)
%           depl -  Nodal displacements
%           -----> depl = [UX UY UZ]
%           component - The components whose profile to be plotted
%           -----> components = a column vector in the order of node
%                               numbers
%
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   
%
% version 1: 28 August 2011
% Version 2: 16 September 2016
%--------------------------------------------------------------------------
dimension = size(coordinates,2) ;  % Dimension of the mesh
nel = size(nodes,1) ;                       % number of elements
nnode = length(coordinates);               % total number of nodes in system
nnel = size(nodes,2);                     % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ; UX = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ; UY = zeros(nnel,nel) ;
Z = zeros(nnel,nel) ; UZ = zeros(nnel,nel) ;
profile = zeros(nnel,nel) ;
%
if dimension == 3   % For 3D plots
    ux = components(:,1) ;
    uy = components(:,2) ;
    uz = components(:,3) ;
    if nnel == 4 % surface in 3D 
        for iel=1:nel   
            nd=nodes(iel,:);         % extract connected node for (iel)-th element
            X(:,iel)=coordinates(nd,1);    % extract x value of the node
            Y(:,iel)=coordinates(nd,2);    % extract y value of the node
            Z(:,iel)=coordinates(nd,3) ;   % extract z value of the node
            
            UX(:,iel) = ux(nd') ;         % extract displacement value's of the node 
            UY(:,iel) = uy(nd') ;
            UZ(:,iel) = uz(nd') ;
            profile(:,iel) = magnitude(nd') ;  
        end
        % Plotting the profile of a property on the deformed mesh
        defoX = X+factor*UX ;
        defoY = Y+factor*UY ;
        defoZ = Z+factor*UZ ;
        %figure
        fill3(defoX,defoY,defoZ,profile)
        title('Profile of component on deformed Mesh') ;       
        rotate3d on ;
        axis off ;
        % Colorbar Setting
        if clim == [0 0]
            clim = caxis;
        end
        SetColorbar(clim, cbar_type, cbar_label)
    elseif nnel==8  % solid in 3D
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,nel) ;
        profile = XYZ ;
        for e=1:nel
            nd=nodes(e,:);
            X = coordinates(nd,1)+factor*ux(nd) ;
            Y = coordinates(nd,2)+factor*uy(nd) ;
            Z = coordinates(nd,3)+factor*uz(nd) ;
            XYZ{e} = [X Y Z] ;
            profile{e} = magnitude(nd);
        end
        %figure
        cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,.......
            repmat({'Faces'},1,nel),repmat({fm},1,nel),......
            repmat({'FaceVertexCdata'},1,nel),profile,......
            repmat({'FaceColor'},1,nel),repmat({'flat'},1,nel));
        view(3)
        %set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;
        % Colorbar Setting
        if clim == [0 0]
            clim = caxis;
        end
        SetColorbar(clim, cbar_type, cbar_label)
    elseif nnel==20  % solid in 3D
        fm = [1 9 2 18 6 13 5 17; 2 10 3 19 7 14 6 18; 3 11 4 20 8 15 7 19; 4 12 1 17 5 16 8 20; 1 9 2 10 3 11 4 12; 5 13 6 14 7 15 8 16];
        XYZ = cell(1,nel) ;
        profile = XYZ ;
        for e=1:nel
            nd=nodes(e,:);
            X = coordinates(nd,1)+factor*ux(nd) ;
            Y = coordinates(nd,2)+factor*uy(nd) ;
            Z = coordinates(nd,3)+factor*uz(nd) ;
            XYZ{e} = [X Y Z] ;
            profile{e} = magnitude(nd); 
        end
        
        %figure
        cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,.......
            repmat({'Faces'},1,nel),repmat({fm},1,nel),......
            repmat({'FaceVertexCdata'},1,nel),profile,......
            repmat({'FaceColor'},1,nel),repmat({'interp'},1,nel),...
            repmat({'CDataMapping'},1,nel),repmat({'scaled'},1,nel));
        view(3)
        %set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;
        % Colorbar Setting
        if clim == [0 0]
            clim = caxis;
        end
        SetColorbar(clim, cbar_type, cbar_label)
    end
elseif dimension == 2           % For 2D plots
    ux = components(:,1) ;
    uy = components(:,2) ;
    for iel=1:nel   
        nd=nodes(iel,:);         % extract connected node for (iel)-th element
        X(:,iel)=coordinates(nd,2);    % extract x value of the node
        Y(:,iel)=coordinates(nd,3);    % extract y value of the node
        
        UX(:,iel) = ux(nd') ;
        UY(:,iel) = uy(nd') ;
        profile(:,iel) = magnitude(nd') ;      
    end
    % Plotting the profile of a property on the deformed mesh
    defoX = X+factor*UX ;
    defoY = Y+factor*UY ;  
    %figure
    plot(defoX,defoY,'k')
    fill(defoX,defoY,profile)
    title('Profile of UX on deformed Mesh') ;      
    axis off ;
    % Colorbar Setting
    if clim == [0 0]
        clim = caxis;
    end
    SetColorbar(clim, cbar_type, cbar_label)
end
