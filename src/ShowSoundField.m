function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

    figure;
    disp('plot the transmission loss field!');
    pcolor( r./1000, z, tl); hold on;
    caxis( [tlmin tlmax] ); colormap( flipud(jet) );
    shading flat; colorbar; view( 0, -90 );
    xlabel( 'Range (km)')  ; ylabel( 'Depth (m)');
    colorbar( 'YDir', 'Reverse' );title(casename);
    set(gca,'FontSize',20,'FontName','Times New Roman');

end