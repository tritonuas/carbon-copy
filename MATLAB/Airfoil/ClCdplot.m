function ClCdplot(iter,aSet,bSet,ccSet,ClCd,Amax,Bmax,CCmax)
    hold on
    if iter > 1
        scatter3(aSet(iter-1),bSet(iter-1),ccSet(iter-1),40,ClCd(iter-1),'filled');
    end
    scatter3(aSet(iter),bSet(iter),ccSet(iter),40,ClCd(iter),'filled');
    scatter3(aSet(iter),bSet(iter),ccSet(iter),25,ClCd(iter),'ro');
    
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