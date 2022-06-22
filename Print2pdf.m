%%
%%
%% I print to working directory
%%
str = 'Reproduced';
filename = ['ReWorkedFigures/',str];
%
L = gcf; 
set(L,'Units','Inches');
pos = get(L,'Position');
set(L,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(L,filename,'-dsvg','-r0')
print(L,filename,'-dpdf','-r0')
saveas(L,filename,'fig')

filenamePdf = ['pdfcrop ', filename,'.pdf']
system(filenamePdf)

%%
h=gca; 
h.YAxis.FontSize =15
h.XAxis.FontSize =15
h.Legend.FontSize =15
h.YLabel.FontSize=22
h.XLabel.FontSize=22
