function ClCdplot(Cl, GuessAirfoil, iter, aSet, bSet, ccSet, clcdSet, maxA, maxB, maxCC)
    
    hold on
    if iter > 1 %Cover last point's red outline
        scatter3(aSet(iter-1),bSet(iter-1),ccSet(iter-1),40,clcdSet(iter-1),'filled'); 
    end
    
    scatter3(aSet(iter),bSet(iter),ccSet(iter),40,clcdSet(iter),'filled');
    scatter3(aSet(iter),bSet(iter),ccSet(iter),20,clcdSet(iter),'ro');
%     scatter3(aSet(iter),bSet(iter),ccSet(iter),40,clcdSet(iter),'filled');
    
    view(-31,14);
    axis([0 maxA 0 maxB 0 maxCC]);
    xlabel('Max Camber (A)');
    ylabel('Location of Max Camber (B)');
    zlabel('Max Thickness (CC)');
    title('ClCd vs. NACA Variables');
    subtxt = ['For Cl = ' num2str(Cl) ', GuessAirfoil = ' GuessAirfoil];
    subtitle(subtxt);
    
    c = colorbar;   % create and label the colorbar
    caxis([0,100]); 
    c.Label.String = 'ClCd';
    
    print(['C:\TUAS\carbon-copy\MATLAB\Airfoil\Plot_images\plotimg' GuessAirfoil 'cl' num2str(Cl) 'iter' num2str(iter)  '.jpg'],'-djpeg')
end