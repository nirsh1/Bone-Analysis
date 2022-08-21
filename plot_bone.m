function p=plot_bone(V,F,C,Rotate)
%
% Dr. Nir Shvalb, Ariel University, Nirsh@ariel.ac.il
% Bone toolbox July 2016
%
% V,F,C are vertices Faces and color, Rotate=1 or 0
% usage:
% hh = msgbox('Thats the plane you picked'); 
%child = get(hh,'Children');
%delete(child(2))
%pause(1)
%close(hh)
%clf
%planeplot(VVV4,VVV5,VVV6)
%p=plot_bone(V,F,C,1);
%%
    hold on
    p = patch('faces', F, 'vertices' ,V);
    set(p, 'facec', 'flat');            % Set the face color flat
    set(p, 'FaceVertexCData', C);       % Set the color (from file)
    set(p, 'EdgeColor','none');         % Set the edge color
    light                               % add a default light
    light('Position',[1 0 0],'Style','infinite');
    light('Position',[0 1 0],'Style','infinite');
    light('Position',[0 0 1],'Style','infinite');
    daspect([1 1 1])                    % Setting the aspect ratio
    view(3)                             % Isometric view
    xlabel('X'),ylabel('Y'),zlabel('Z')
    title(['Imported CAD data from file'])
    set(gcf,'units','normalized','position',[0 0 1 .8])
    grid on
    mid=sum(V)'/length(V);
    R=(V'-mid*ones(1,length(V(:,1)))).^2;
    vsize=max(sqrt(sum(R)));
       axis([-vsize+mid(1) vsize+mid(1) -vsize+mid(2) vsize+mid(2) -vsize+mid(3) vsize+mid(3)]);
    axis([]);
    set(gcf,'units','normalized','position',[0 0 1 0.8])
    view(0,0)
    if Rotate==1
            for i=1:2:90 % rotate the bone to view from all angles.
            view([i,0])
            pause(0.0001)
            end   
            for i=90:-2:0 % rotate the bone to view from all angles.
            view([i,0])
            pause(0.0001)
            end 
    end
end
