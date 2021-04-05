function ClCdplot(A,B,CC,ClCd)
    hold on
    scatter3(A,B,CC,40,ClCd,'filled');
    view(-31,14);
    xlabel('Max Camber');
    ylabel('Location of Max Camber');
    zlabel('Max Thickness');
    
    cb = colorbar;   % create and label the colorbar
    caxis([0,100]); 
    clabel = 'ClCd';
    
    
end