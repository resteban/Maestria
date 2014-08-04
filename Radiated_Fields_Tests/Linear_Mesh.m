function [x]=Linear_Mesh(Xmin,Xmax,DeltaXmax,DeltaXmin,N)

DeltaXincrement=(DeltaXmax-DeltaXmin)/N;    %%%% definition of the increment mesh
Xm=Xmin+0.5*(Xmax-Xmin);                    %%%% definition of the middle of the interval to mesh



Xo=Xmin;

xb=Xo;
x1(1)=Xo;
x1(2)=Xo+DeltaXmin;
i=2;

while xb<Xm-DeltaXmax
    
% %         x1(i)=x1(i-1)+(DeltaXmin+DeltaXincrement)*i;
        x1(i)=Xo+(i-1)*DeltaXmin+((i-1)*(i-2)*DeltaXincrement)/2;

            if (x1(i)-x1(i-1))>DeltaXmax
                x1(i)=x1(i-1)+DeltaXmax;
            end

        xb=x1(i);      
        i=i+1;  
end

x1(length(x1))=x1(length(x1)-1)+0.5*(Xm-x1(length(x1)-1));
x1(length(x1)+1)=Xm;

delta_xb=fliplr(diff(x1));
xb=max(x1);
for i=1:length(delta_xb)
    x2(i)=xb+delta_xb(i);
    xb=x2(i);
end

x=[x1 x2];