function [x]=Geometric_Mesh(Xmin,Xmax,DeltaXmax,DeltaXmin,alfa)

Nmax=log10(DeltaXmax/DeltaXmin)/log10(alfa);
Xm=Xmin+0.5*(Xmax-Xmin);


xb=Xmin;
i=1;
while xb<Xm
        x1(i)=xb;
        if i<Nmax
            delta=DeltaXmin*(alfa^i);
        else
            delta=DeltaXmin*(alfa^Nmax);
        end
            xb=xb+delta;
            i=i+1;         
end
x1(i)=Xm;
delta_xb=fliplr(diff(x1));
xb=max(x1);
for i=1:length(delta_xb)
    x2(i)=xb+delta_xb(i);
    xb=x2(i);
end

x=[x1 x2];