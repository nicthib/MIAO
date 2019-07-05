function genchoose_PLOT(x,y,i,col,lw)
line([x(i),x(i)],[y(i),y(i+1)],'Color',col,'LineWidth',lw); hold on
line([x(i),x(i+1)],[y(i),y(i)],'Color',col,'LineWidth',lw);
line([x(i),x(i+1)],[y(i+1),y(i+1)],'Color',col,'LineWidth',lw);
line([x(i+1),x(i+1)],[y(i),y(i+1)],'Color',col,'LineWidth',lw);
