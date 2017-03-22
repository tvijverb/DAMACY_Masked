function giffy( filename )
%GIFFY Summary of this function goes here
%   Detailed explanation goes here
campos([2.5 -1 180]);
axis off;
for n = 1:1:181
camorbit(2,0,'camera')
drawnow;
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

end

