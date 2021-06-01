function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

    figure;
    disp('plot the transmission loss field!');
    pcolor( r./1000, z, tl); colorbar( 'YDir', 'Reverse' );
    caxis( [tlmin tlmax] ); colormap( flipud(jet) );
    shading flat; view( 0, -90 );title(casename);
    xlabel( 'Range (km)')  ; ylabel( 'Depth (m)');  
    set(gca,'FontSize',20,'FontName','Times New Roman');

end