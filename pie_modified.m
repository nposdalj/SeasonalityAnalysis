function pie_modified(dat,tilecolor)
%PIE_MODIFIED    Pie chart.
%   PIE_MODIFIED(DAT,TILECOLOR) draws a pie plot of the data in the vector
%   DAT using the colors in TILECOLOR. It does not slip colors for zero
%   values inputs. Note that this function will work only if the length of
%   the dat and tilecolor are the same.
    No = length(dat);
    if( No ~= length(tilecolor))
        error('Length of tilecolor and dat should match')
    end
    h = pie(dat);
    set(gca,'clim',[0.5 No+0.5])
    colormap(tilecolor);
    z =~(dat == 0);
    j = 1:2:2*No;
    index = 0;
    for i = 0:No-1
        if(z(i+1) ~= 0)
            index = index +1;
            set(h(j(index)),'Cdata',i+1)
         end
    end
end
