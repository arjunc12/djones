clear; 

load rs.d

nx=74; 
ny=49; 
xend=3; 
yend=2; 
sl=9;           % no of frames

x=0:xend/nx:xend; 
y=0:yend/ny:yend; 
%x=0:75:3
%y=0:50:2

sx=size(x); 
sy=size(y); 

for k=1,sl
for i=1:sx(2)
for j=1:sy(2)
h(i,j,k)=rs(j+sy(2)*(i-1)); 
end 
end 
end

%figure(1); 

%p=pcolor(y,x,h(:,:,k));
%shading interp 
%theaxis2=axis;            % uses the axis of s always
%axis image 
%title('Height') 
%xlabel('r - axis') 
%ylabel('x - axis') 
%colorbar('vertic')
%colorbar('horiz') 

     M = moviein(sl);          % sl=number of frames
    for k=1:sl, pcolor(h(:,:,k)'), axis('equal'),shading('interp'),  
    caxis([74,49]); colorbar, M(:,k) = getframe; end

%     set(p,'ZData',h(:,:,k),'CData',h(:,:,k));
%     axis(theaxis2);
%     Mh(:,k) = getframe;
%    end

     movie(M,2)           % Shows the movie 2 times
