
%normal function (save this in calc_RMSE.m)
function rmse=calc_RMSE(a,b)
rmse=sqrt(mean((a(:)-b(:)).^2));