fig_name = 'C:\Users\SON\Desktop\BSI\Algo\ML\new\ZF';

fig = openfig([fig_name '.fig']);

h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles t
dataObjs = dataObjs{2};

xdata = get(dataObjs, 'XData'); 
ydata = get(dataObjs, 'YData');

save('C:\Users\SON\Desktop\BSI\Algo\ML\new\ZF_x.mat', 'xdata');
save('C:\Users\SON\Desktop\BSI\Algo\ML\new\ZF_y.mat', 'ydata');