function hila12

% upload an STL file
clc
close all
F = findall(0,'type','figure');
delete(F);
warning off
[FileName,path]=uigetfile('*.stl','Select the STL code file');
[F,V]=cad2matdemo(FileName,path);
grid on
Cont='Yes';
k=0;
[center,rads]=find_sphere(V(1:3,:)',F,zeros(length(V),1))

while strcmp(Cont,'Yes')
    
task = questdlg('pick a minimal plane by','Pick a cut','one point','two points','end','one point');

switch task
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'end'
        return       
    case 'one point'
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=helpdlg(['Pick an approximated plane by selecting three points (1st pt will remain fixed)']);
waitfor(h); 
i=0;x=0;y=0;z=0; 
for pts=1:3
            pause
            [P VV VI] = select3d;
            i=i+1;
            hold on;
            plot3(VV(1),VV(2),VV(3),'+b');
            x(i)=VV(1);
            y(i)=VV(2);
            z(i)=VV(3);
end
plot3(x(1),y(1),z(1),'ob');
%text(x(1),y(1),z(1),'anchore point')
k=k+1;
find_plane1(V,x,y,z,k)
case 'two points'
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=helpdlg(['Pick an approximated plane by selecting three points (1st and 2nd pts will remain fixed)']);
waitfor(h); 
i=0;x=0;y=0;z=0; 
for pts=1:3
            pause
            [P VV VI] = select3d;
            i=i+1;
            hold on;
            plot3(VV(1),VV(2),VV(3),'+b');
            x(i)=VV(1);
            y(i)=VV(2);
            z(i)=VV(3);
end
plot3(x(1),y(1),z(1),'ob');
plot3(x(2),y(2),z(2),'ob');
%text(x(1),y(1),z(1),'anchore point')
k=k+1;
find_plane2(V,x,y,z,k)

Cont=questdlg('Shall we procceed to the next cut ?','Next cut','Yes','No','Yes');        
end
end
end


function find_plane1(V,x,y,z,k);
    % V is the bone vertices
        % (x,y,z) is the selected cut points
        % the function creates an obj fule of the path passing through the extermal points
    
    
    kind = questdlg('What do you need ? ','min or max area','minimal','maximal','minimal');        
    waitfor(kind)
    
    
    % finds indices of V of chosen waypoints 
    for jj=1:3
        ind=find(sum((V(1:3,:)'-ones(length(V),1)*[x(jj),y(jj),z(jj)]).^2,2)==0);
        I(jj)=ind(1);
    end
    p1=V(1:3,I(1));
    p2=V(1:3,I(2));
    p3=V(1:3,I(3)); 
    
    % finding a 90deg vectors eminanting from p1
    [mid,rad,v1,v2] = circlefit3d(p1',p2',p3'); % v1,v2 are perp |rad| length vectrs (can be thought of starting from mid)
    alpha1=sum((p1'-mid).*v1);
    beta1=sum((p1'-mid).*v2);
    beta2=-alpha1;%(-2*beta1-sqrt(beta1^2+4*alpha1^2))/(2*alpha1);
    alpha2=beta1;%sqrt(rad^2-beta2^2);
   
    p22=v2'*alpha2+v1'*beta2+mid'; % a  vector (starting from the origin) ending on the circle
    p33=-v2'*alpha2-v1'*beta2+mid';% a  vector (starting from the origin)  ending on the circle s.t. p22-mid and p33-mid are perp
    
    search=30;
    ds=search/5;
    iii=0;
    for a1=-search:ds:search
        iii=iii+1;
        jjj=0;
        for a2=-search:ds:search
            jjj=jjj+1;
            [AA,pp,VV]=cut_bone(p22,p33,p1,V,a1,a2);
    
    Area(iii,jjj)=AA;
    projx(iii,jjj,:)=pp(1,:);
    projy(iii,jjj,:)=pp(2,:);
    V_circx(iii,jjj,:)=VV(1,:);
    V_circy(iii,jjj,:)=VV(2,:);
    V_circz(iii,jjj,:)=VV(3,:);
      [a1,a2]
     end
    end
 figure
 
    [tx,ty]=meshgrid(-search:ds:search);
    [tx_fine,ty_fine]=meshgrid(-search:ds/10:search);
    Area_fine=interp2(tx,ty,Area,tx_fine,ty_fine);
   contour3(tx_fine,ty_fine,Area_fine,300,'g')
      axis square
hold on
switch kind
        case 'minimal'
   [Area_m,I]=min(Area_fine(:));
   [iii,jjj] = ind2sub(size(Area_fine),I);
   m_ty_est=(iii-1)*ds/10-search;
   m_tx_est=(jjj-1)*ds/10-search;
   plot3(m_tx_est,m_ty_est,Area_m,'or')

        case 'maximal'
   [Area_m,I]=max(Area_fine(:));
   [iii,jjj] = ind2sub(size(Area_fine),I);
   m_ty_est=(iii-1)*ds/10-search;
   m_tx_est=(jjj-1)*ds/10-search;
   plot3(m_tx_est,m_ty_est,Area_m,'or')
end
   
 

 figure(2)
 title(['selected plane. Area=',num2str(Area_m)]);
 
            [AA,pp,VV]=cut_bone(p22,p33,p1,V,m_ty_est,m_tx_est);
 figure(1)
 hold on
    % smoothing the path
   % plot3(VV(1,:),VV(2,:),VV(3,:),'r')

        M=length(VV);
    m=8;
  %  plot3(VV(1,1:m:M),VV(2,1:m:M),VV(3,1:m:M),'g')
    V_circ(1,:)=interp1(1:m:(M+2*m),[VV(1,1:m:M),VV(1,m),VV(1,2*m)],1:M,'spline');
    V_circ(2,:)=interp1(1:m:(M+2*m),[VV(2,1:m:M),VV(2,m),VV(2,2*m)],1:M,'spline');
    V_circ(3,:)=interp1(1:m:(M+2*m),[VV(3,1:m:M),VV(3,m),VV(3,2*m)],1:M,'spline');
   hold on;
   plot3(V_circ(1,:),V_circ(2,:),V_circ(3,:),'r')
	 title(['selected plane. Area=',num2str(Area_m)]);

    
   % figure(2)
    %plot(proj(1,:),proj(2,:))
    %axis equal
    %title(['The area is ',num2str(Area), 'which is ',txt, ' of all the ancored cross sections'])
    %figure(1)
    [FileName,PathName] = uiputfile('*.obj');
    fname=[PathName,FileName];
    V_circ=round(V_circ*1000)/1000;
    writeobj(fname,V_circ(:,1:end-1));
end


function find_plane2(V,x,y,z,k);
    % V is the bone vertices
        % (x,y,z) is the selected cut points
        % the function creates an obj fule of the path passing through the extermal points
    
    
    kind = questdlg('What do you need ? ','min or max area','minimal','maximal','minimal');        
    waitfor(kind)
    
    
    % finds indices of V of chosen waypoints 
    for jj=1:3
        ind=find(sum((V(1:3,:)'-ones(length(V),1)*[x(jj),y(jj),z(jj)]).^2,2)==0);
        I(jj)=ind(1);
    end
    p1=V(1:3,I(1));
    p2=V(1:3,I(2));
    p3=V(1:3,I(3)); 
    
    % finding a vector perp to p1 and p2
        [mid,rad,v1,v2] = circlefit3d(p1',p2',p3');
aa1=p1;
aa2=p2;
aa3=p3;

p1=aa1;
p2=aa2;
p3=aa3;

        perp12=(p3)-(sum((p2-p1).*(p3))/norm(p1-p2))*(p2-p1)/norm(p2-p1);
        p33=perp12/norm(perp12)*rad+mid'; % a  vector (starting from the origin)  ending on the circle 
        p22=p2;      % a  vector (starting from the origin)  ending on the circle      s.t. p22-mid and p33-mid are perp

     
    search=20;
    ds=search/20;
    iii=0;
    for ang=-search:ds:search
        iii=iii+1;
         [AA,pp,VV]=cut_bone(p22,p33,p1,V,0,ang);
    %plot3(VV(1,:),VV(2,:),VV(3,:),'r')
    Area(iii)=AA;
    projx(iii,:)=pp(1,:);
    projy(iii,:)=pp(2,:);
    V_circx(iii,:)=VV(1,:);
    V_circy(iii,:)=VV(2,:);
    V_circz(iii,:)=VV(3,:);
    end
 figure
 
    tx_fine=-search:ds:search;
    
    %tx_fine=-search:ds/10:search;
    %Area_fine=interp1(tx,Area,tx_fine);
        Area_fine=Area;
       plot(tx_fine,Area_fine,'g')
   
      axis square
hold on
switch kind
        case 'minimal'
   [Area_m,I]=min(Area_fine(:));
   jjj = ind2sub(size(Area_fine),I);
 m_tx_est=(jjj-1)*ds-search;
 plot(m_tx_est,Area_m,'or')

        case 'maximal'
   [Area_m,I]=max(Area_fine(:));
   jjj = ind2sub(size(Area_fine),I);
 m_tx_est=(jjj-1)*ds-search;
 plot(m_tx_est,Area_m,'or')
end
    
 figure(2)
 title(['selected plane. Area=',num2str(Area_m)]);
 
            [AA,pp,VV]=cut_bone(p22,p33,p1,V,m_tx_est,0);
 figure(1)
 hold on
% smooth path

                 M=length(VV);
                m=8;
                  %plot3(VV(1,1:m:M),VV(2,1:m:M),VV(3,1:m:M),'g')
                V_circ(1,:)=interp1(1:m:(M+2*m),[VV(1,1:m:M),VV(1,m),VV(1,2*m)],1:M,'spline');
                V_circ(2,:)=interp1(1:m:(M+2*m),[VV(2,1:m:M),VV(2,m),VV(2,2*m)],1:M,'spline');
                V_circ(3,:)=interp1(1:m:(M+2*m),[VV(3,1:m:M),VV(3,m),VV(3,2*m)],1:M,'spline');
                  plot3(V_circ(1,:),V_circ(2,:),V_circ(3,:),'r')
   hold on;
   %plot3(VV(1,:),VV(2,:),VV(3,:),'r')            
	 title(['selected plane. Area=',num2str(Area_m)]);

    
   % figure(2)
    %plot(proj(1,:),proj(2,:))
    %axis equal
    %title(['The area is ',num2str(Area), 'which is ',txt, ' of all the ancored cross sections'])
    %figure(1)
    [FileName,PathName] = uiputfile('*.obj');
    fname=[PathName,FileName];
        writeobj(fname,V_circ(:,1:end-1));
end








    function [AA,pp,VV]=cut_bone(p22,p33,p1,V,a1,a2)
        
% finds points on the plane
  
    nrm=cross(p22-p1,p33-p1)/norm(cross(p22-p1,p33-p1));
    p2=p22+nrm*norm(p1-p22)*tan(a1*pi/180);
    p3=p33+nrm*norm(p1-p33)*tan(a2*pi/180);
    circ=[];
    [mid,rad,v1,v2] = circlefit3d(p1',p2',p3');
    t=1:361; aa = t/180*pi; 
    circ(1,:) =sin(aa)*rad.*v1(1)+cos(aa)*rad.*v2(1); 
    circ(2,:) =sin(aa)*rad.*v1(2)+cos(aa)*rad.*v2(2); 
    circ(3,:) =sin(aa)*rad.*v1(3)+cos(aa)*rad.*v2(3); 
        
    % finds directions of all vertices
        hypotenuse_V=sqrt(sum((V(1:3,:)-mid'*ones(1,length(V))).^2));
        x_V=(V(1,:)-mid(1))./hypotenuse_V;
        y_V=(V(2,:)-mid(2))./hypotenuse_V;
        z_V=(V(3,:)-mid(3))./hypotenuse_V;
    
    pp=[];
    for kk=1:length(circ) % finds vertices that match a point on the circlular path
            dotV_vec=(circ(1,kk)*x_V+circ(2,kk)*y_V+circ(3,kk)*z_V); 
            ind_rel=find(dotV_vec==min(dotV_vec(:)));
            ind_rel=ind_rel(1); 
            VV(1:3,kk)=V(1:3,ind_rel);
            pp(1,kk)=sum(VV(1:3,kk).*v1');
            pp(2,kk)=sum(VV(1:3,kk).*v2');
    end
    AA=polyarea(pp(1,:),pp(2,:));
    end
    function writeobj(name,V_circ);
    fid = fopen(name,'w');
for i=1:length(V_circ)
%fprintf(fid,'v %f %f\n',proj(1,i),proj(2,i)); % prints the projected curve  (to the plane) into the obj file
fprintf(fid,'v %f %f %f \n',V_circ(1,i),V_circ(2,i),V_circ(3,i));
end
for i=1:length(V_circ)-1
%fprintf(fid,'v %f %f\n',v(1,i),v(2,i));
fprintf(fid,'e %d %d  \n',i,i+1);
end
fprintf(fid,'e %d %d  \n',i+1,1);
%fprintf(fid,'g\n');
x=num2str(mean(V_circ(1,:)));
y=num2str(mean(V_circ(2,:)));
z=num2str(mean(V_circ(3,:)));
disp(['The curve center is=[',x,',',y,',',z,']'])
fclose(fid);
    end
function cut_obj(V,F,x,y,z,k);
    % V is the bone vertices
    % F is the bone faces
        % (x,y,z) is the selected cut points
        % k is the cut number created.
        % the function creates an obj file of the path passing through the extermal points
    
    % creating weight on the bone
    N=length(x);
    ax=sum(x)/N;
    ay=sum(y)/N;
    az=sum(z)/N;
    mid=[ax,ay,az];
    plot3(ax,ay,az,'ro')
    R=(V(1:3,:)'-ones(length(V),1)*[ax,ay,az]).^2;
    W=sqrt(sum(R')');
    kind = questdlg('what kind of cut should it be ? concave/convex = waste/chest','kind of cut','concave','convex','convex');        
    waitfor(kind)
    switch kind
        case 'concave'
            W=max(W)+1-W;
    end
    
    h=msgbox(['Creating file SpatialCrossSection',num2str(k),'.OBJ.      wait.....  Next we shall procceed to the next cut.']);

    % finds indices of V of chosen waypoints 
    for jj=1:N
        ind=find(sum((V(1:3,:)'-ones(length(V),1)*[x(jj),y(jj),z(jj)]).^2,2)==0);
        I(jj)=ind(1);
    end
    
    % finds paths between waypoints
    path=[];
    t=0:0.1:1;
    for ii=1:N
    WP1=V(1:3,I(ii));
    WP2=V(1:3,I(round(mod(ii+1,N+0.1))));
    p=WP2*t+WP1*(1-t);  % a path = straight line between WP2 and WP1
    nrm=cross((mid'-WP1),(mid'-WP2)); % normal of mid-WP1-WP2 plane
    nrm=nrm/norm(nrm);
    
    
    for jj=1:length(t)  % for each point on the path  between way points performs
        
        tp=[-.1:.01:.1];  % outlines a perpendicular segment to that point 
        vec=p(1:3,jj)*ones(1,length(tp))+nrm*norm(WP1-WP2)*tp; % to visual: 
        plot3(vec(1,:),vec(2,:),vec(3,:),'r')
        
        % finds the directions perpendicular to the perpendicular path
        hypotenuse_vec=sqrt(sum((vec-mid'*ones(1,length(vec))).^2)); % the radial distance of vec from mid (the circle's center)
        x_vec=(vec(1,:)-mid(1))./hypotenuse_vec; % (x_vec,y_vec,z_vec) is the unit vector from the center to vec
        y_vec=(vec(2,:)-mid(2))./hypotenuse_vec;
        z_vec=(vec(3,:)-mid(3))./hypotenuse_vec;
        
        % finds directions of all vertices of the bone (no need to draw its
        % the locations of all vertices with the origin=mid
        hypotenuse_V=sqrt(sum((V(1:3,:)-mid'*ones(1,length(V))).^2));
        x_V=(V(1,:)-mid(1))./hypotenuse_V;
        y_V=(V(2,:)-mid(2))./hypotenuse_V;
        z_V=(V(3,:)-mid(3))./hypotenuse_V;
        
        W0=-100;
        for kk=1:length(x_vec) % for each point in the perpendicular path performs
            % finds vertices that match a point on the perpendicular path
            dotV_vec=(x_vec(kk)*x_V+y_vec(kk)*y_V+z_vec(kk)*z_V); 
            ind_rel=find(dotV_vec==min(dotV_vec(:)));
            ind_rel=ind_rel(1); 
            V_veckk=V(1:3,ind_rel);
            
            if W(ind_rel)>W0  % selects an extermal vertex out of these vertices
                RelV=V_veckk;
                W0=W(ind_rel);
            end
        end
        path=[path, RelV];
    end
    end
    
    %hold on;
    plot3(path(1,:),path(2,:),path(3,:),'LineWidth',4);
    fname=['SpatialCrossSection',num2str(k),'.OBJ'];
    writeobj(fname,path);
    close(h);
end

function [F, V]=cad2matdemo(filename,path)
% CAD2MATDEMO, a demonstration of importing 3D CAD data into Matlab.
% To get CAD data into Matlab, the process is:
%
% 1) Export the 3D CAD data as an ASCII STL (or Pro/E render SLP) file.
% 2) This Matlab routine reads the CAD data
% 3) Once read, the CAD data is rotated around a bit.
%
% Program has been tested with: AutoCAD, Cadkey, and Pro/Engineer.
% Should work with most any CAD programs that can export STL.
% 
% Format Details:  STL is supported, and the color version of STL
% that Pro/E exports, called 'render.'  The render (SLP) is just 
% like STL but with color added.
% 
% Note: This routine has both the import function and some basic
% manipulation for testing.  The actual reading mechanism is located
% at the end of this file.
  
if nargin == 0    
   filename = 'hook.stl'; % a simple demo part
    %filename = 'Femur_BLVI_9215_DecX3.stl';
   warning(['No file specified, using demo file: ' filename]);
end
%
% Read the CAD data file:
[F, V, C] = rndread(filename,path);
clf;
  p = patch('faces', F, 'vertices' ,V);
    %set(p, 'facec', 'b');              % Set the face color (force it)
    set(p, 'facec', 'flat');            % Set the face color flat
    set(p, 'FaceVertexCData', C);       % Set the color (from file)
    %set(p, 'facealpha',.4)             % Use for transparency
    set(p, 'EdgeColor','none');         % Set the edge color
    %set(p, 'EdgeColor',[1 0 0 ]);      % Use to see triangles, if needed.
    light                               % add a default light
    light('Position',[1 0 0],'Style','infinite');
    light('Position',[0 1 0],'Style','infinite');
    light('Position',[0 0 1],'Style','infinite');
    daspect([1 1 1])                    % Setting the aspect ratio
    view(3)                             % Isometric view
    xlabel('X'),ylabel('Y'),zlabel('Z')
    title(['Imported CAD data from ' filename])
    drawnow                             %, axis manual
    %
    %
    % Move it around.
    % To use homogenous transforms, the n by 3 Vertices will be turned to 
    % n by 4 vertices, then back to 3 for the set command.
    % Note: n by 4 needed for translations, not used here, but could, using tl(x,y,z)
    V = V';
    V = [V(1,:); V(2,:); V(3,:); ones(1,length(V))];
    %
    mid=sum(V')/length(V);
    R=(V(1:3,:)'-ones(length(V),1)*mid(1:3)).^2;
    vsize=max(sqrt(sum(R')'));
    %vsize = maxv(V); %attempt to determine the maximum xyz vertex. 
    axis([-vsize+mid(1) vsize+mid(1) -vsize+mid(2) vsize+mid(2) -vsize+mid(3) vsize+mid(3)]);
    axis([]);
    %
    grid on;
for ang = 0:5:360
    nv = rx(ang)*(V-mean(V')'*ones(1,length(V)))+mean(V')'*ones(1,length(V));
   set(p,'Vertices',nv(1:3,:)')       
   drawnow
end
for ang = 0:5:360
    nv = ry(ang)*(V-mean(V')'*ones(1,length(V)))+mean(V')'*ones(1,length(V));
   set(p,'Vertices',nv(1:3,:)')       
   drawnow
end
    set(p,'Vertices',V(1:3,:)')      
    drawnow
end
%
% End of main routine, here are the functions used:
% Homogeneous manipulation functions follow:
%

function Rx = rx(THETA)
% ROTATION ABOUT THE X-AXIS
%
% Rx = rx(THETA)
%
% This is the homogeneous transformation for
% rotation about the X-axis.
%
%	    NOTE:  The angle THETA must be in DEGREES.
%
THETA = THETA*pi/180;  % Note: THETA in radians.
c = cos(THETA);
s = sin(THETA);
Rx = [1 0 0 0; 0 c -s 0; 0 s c 0; 0 0 0 1];
end
%
function Ry = ry(THETA)
% ROTATION ABOUT THE Y-AXIS
%
% Ry = ry(THETA)
%
% This is the homogeneous transformation for
% rotation about the Y-axis.
%
%		NOTE: The angel THETA must be in DEGREES.
%
THETA = THETA*pi/180;  %Note: THETA is in radians.
c = cos(THETA);
s = sin(THETA);
Ry = [c 0 s 0; 0 1 0 0; -s 0 c 0; 0 0 0 1];
end
%
function Rz = rz(THETA)
% ROTATION ABOUT THE Z-AXIS
%
% Rz = rz(THETA)
%
% This is the homogeneous transformation for
% rotation about the Z-axis.
%
%		NOTE:  The angle THETA must be in DEGREES.
%
THETA = THETA*pi/180;  %Note: THETA is in radians.
c = cos(THETA);
s = sin(THETA);
Rz = [c -s 0 0; s c 0 0; 0 0 1 0; 0 0 0 1];
end
%
function T = tl(x,y,z)
% TRANSLATION ALONG THE X, Y, AND Z AXES
%
% T = tl(x,y,z)
%
% This is the homogeneous transformation for
% translation along the X, Y, and Z axes.
%
T = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
end
%
function vsize = maxv(V)
%
% Look at the xyz elements of V, and determine the maximum
% values during some simple rotations.
    vsize = max(max(V));
    % Rotate it a bit, and check for max and min vertex for viewing.
    for ang = 0:10:360
        vsizex = rx(ang)*V;
        maxv = max(max(vsizex));
        if maxv > vsize, vsize = maxv; end
        vsizey = ry(ang)*V;
        maxv = max(max(vsizey));
        if maxv > vsize, vsize = maxv; end
        vsizez = rz(ang)*V;
        maxv = max(max(vsizez));
        if maxv > vsize, vsize = maxv; end
        vsizev = rx(ang)*ry(ang)*rz(ang)*V;
        maxv = max(max(vsizev));
        if maxv > vsize, vsize = maxv; end
    end
end
    
    function openedge=surfedge(f)
%
% openedge=surfedge(f)
%
% find the edge of an open surface or surface of a volume
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/21
%
% input:
%      f: input, surface face element list, dimension (be,3)
%
% output:
%      openedge: list of edges of the specified surface
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(isempty(f))
    openedge=[];
    return;
end

if(size(f,2)==3)
    edges=[f(:,[1,2]);
           f(:,[2,3]);
           f(:,[3,1])];             % create all the edges
elseif(size(f,2)==4)
    edges=[f(:,[1,2,3]);
           f(:,[2,1,4]);
           f(:,[1,3,4]);
           f(:,[2,4,3])];             % create all the edges
else
    error('surfedge only supports 2D and 3D elements');
end
% node4=[f(:,3);f(:,2);f(:,1)];   % node idx concatinated
edgesort=sort(edges,2);
[foo,ix,jx]=unique(edgesort,'rows');

if(isoctavemesh)
        u=unique(jx);
        qx=u(hist(jx,u)==1);
else
        vec=histc(jx,1:max(jx));
        qx=find(vec==1);
end
openedge=edges(ix(qx),:);
    end
% node4=node4(ix(qx));
    function [isoctave verinfo]=isoctavemesh
%
% [isoctave verinfo]=isoctavemesh
%
% determine whether the code is running in octave
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% output:
%   isoctave: 1 if in octave, otherwise 0
%   verinfo: a string, showing the version of octave (OCTAVE_VERSION)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
verinfo='';
isoctave=(exist('OCTAVE_VERSION')~=0);
if(nargout==2 && isoctave)
    verinfo=OCTAVE_VERSION;
end
    end
    
function [fout, vout, cout] = rndread(filename,path)

% Reads CAD STL ASCII files, which most CAD programs can export.
% Used to create Matlab patches of CAD 3D data.
% Returns a vertex list and face list, for Matlab patch command.
% 
% filename = 'hook.stl';  % Example file.
%
dirInfo = dir(path); 
index = strcmp({dirInfo.name},filename);  %
fileSize = dirInfo(index).bytes;
hwait = waitbar(0,'Hila, I am reading the file now. It may take few minutes, drink somthing');


fid = fopen([path,filename]);
fseek(fid, 0, 'eof');
lines_fid  = ftell(fid);
fclose(fid);


fid=fopen([path,filename], 'r'); %Open the file, assumes STL ASCII format.
if fid == -1 
    error('File could not be opened, check name or path.')
end
%
% Render files take the form:
%   
%solid BLOCK
%  color 1.000 1.000 1.000
%  facet
%      normal 0.000000e+00 0.000000e+00 -1.000000e+00
%      normal 0.000000e+00 0.000000e+00 -1.000000e+00
%      normal 0.000000e+00 0.000000e+00 -1.000000e+00
%    outer loop
%      vertex 5.000000e-01 -5.000000e-01 -5.000000e-01
%      vertex -5.000000e-01 -5.000000e-01 -5.000000e-01
%      vertex -5.000000e-01 5.000000e-01 -5.000000e-01
%    endloop
% endfacet
%
% The first line is object name, then comes multiple facet and vertex lines.
% A color specifier is next, followed by those faces of that color, until
% next color line.
%
CAD_object_name = sscanf(fgetl(fid), '%*s %s');  %CAD object name, if needed.
%                                                %Some STLs have it, some don't.   
vnum=0;       %Vertex number counter.
report_num=0; %Report the status as we go.
VColor = 0;
%
dt=0;
while feof(fid) == 0                   % test for end of file, if not then do stuff
    
    

tic    
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
% Check for color
    if strncmpi(fword, 'c',1) == 1;    % Checking if a "C"olor line, as "C" is 1st char.
       VColor = sscanf(tline, '%*s %f %f %f'); % & if a C, get the RGB color data of the face.
    end                                % Keep this color, until the next color is used.
    if strncmpi(fword, 'v',1) == 1;    % Checking if a "V"ertex line, as "V" is 1st char.
       vnum = vnum + 1;                % If a V we count the # of V's
       report_num = report_num + 1;    % Report a counter, so long files show status
       if report_num > 249;
         

    Perc=vnum / (lines_fid/70);
       Trem=lines_fid/250*dt*(1-Perc); %Calculate the time remaining
   Hrs=floor(Trem/3600);Min=floor((Trem-Hrs*3600)/60);
   waitbar(Perc,hwait,[sprintf('%0.1f',Perc*100) ' %, '...
     sprintf('%03.0f',Hrs) ':'...
     sprintf('%02.0f',Min) ':'...
     sprintf('%02.0f',rem(Trem,60)) ' remaining']);
    
        
                   report_num = 0;
                      dt=toc;       

       end
       v(:,vnum) = sscanf(tline, '%*s %f %f %f'); % & if a V, get the XYZ data of it.
       c(:,vnum) = VColor;              % A color for each vertex, which will color the faces.
    end                                 % we "*s" skip the name "color" and get the data.                                          
end
%   Build face list; The vertices are in order, so just number them.
%
fnum = vnum/3;      %Number of faces, vnum is number of vertices.  STL is triangles.
flist = 1:vnum;     %Face list of vertices, all in order.
F = reshape(flist, 3,fnum); %Make a "3 by fnum" matrix of face list data.
%
%   Return the faces and vertexs.
%
fout = F';  %Orients the array for direct use in patch.
vout = v';  % "
cout = c';
%
delete(hwait);
fclose(fid);
end

