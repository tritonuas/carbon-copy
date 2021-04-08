function ClCdplot(iter, A, B, CC, ClCd, Amax, Bmax, CCmax)
    
    aSet(iter) = str2double(A);
    bSet(iter) = str2double(B);
    ccSet(iter) = str2double(CC);  %OK if doesn't have leading zero
    ClCd(iter) = ClCd;
    
    hold on
%     if iter > 1 %Cover last point's red outline
%         scatter3(aSet(iter-1),bSet(iter-1),ccSet(iter-1),40,ClCd(iter-1),'filled'); 
%     end
    
    scatter3(aSet(iter),bSet(iter),ccSet(iter),40,ClCd(iter),'filled');
%     scatter3(aSet(iter),bSet(iter),ccSet(iter),25,ClCd(iter),'ro');
    
    view(-31,14);
    axis([0 Amax 0 Bmax 0 CCmax]);
    xlabel('Max Camber (A)');
    ylabel('Location of Max Camber (B)');
    zlabel('Max Thickness (CC)');
    title('ClCd vs. NACA Variables');
    
    c = colorbar;   % create and label the colorbar
    caxis([0,100]); 
    c.Label.String = 'ClCd';
    
end