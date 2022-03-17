function L = localization_matrix(n,r)

for i = 1:n
   for j = 1:n
      d = min([abs(i-j),abs(n+i-j),abs(n+j-i)]);
      L(i,j) = exp(-(d^2)/(2*r^2)); 
   end
end

end