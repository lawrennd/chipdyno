function chipPlotter(TF, TFerror,npts)
%CHIPPLOTTER plots errorbars in a printable way
figure, f=errorbar([1:npts],TF,TFerror);
set(gca,'linewidth',2)
set(f,'linewidth',2)
set(gca,'fontsize',15)