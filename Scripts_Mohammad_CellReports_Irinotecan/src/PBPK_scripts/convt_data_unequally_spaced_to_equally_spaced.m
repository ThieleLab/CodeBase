function [X,Y] = convt_data_unequally_spaced_to_equally_spaced(x,y,n)
X = linspace(0,x(end),n);
Y = zeros(size(y,1),length(X));
for i=1:size(Y,1)
   for j =1:size(Y,2)
       dummy = find(x>=X(j));
       Y(i,j) = y(i,dummy(1));
   end
end 
end
