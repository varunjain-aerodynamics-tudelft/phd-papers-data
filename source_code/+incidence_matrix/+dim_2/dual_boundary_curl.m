function [curl_lambda_bndry] = dual_boundary_curl(p)

curl_lambda_bndry = zeros((p+1)^2,4*p);

%---------- bttm edges

for i = 2:p
    j = 1;
    
    surf_id = (j-1)*(p+1) + i;
    
    lambda_left = i-1;
    lambda_rght = i;
    
    curl_lambda_bndry(surf_id,lambda_left) = +1;
    curl_lambda_bndry(surf_id,lambda_rght) = -1;
end

%---------- topp edges

for i = 2:p
    j = p+1;
    
    surf_id = (j-1)*(p+1) + i;
    
    lambda_left = 3*p + i-1;
    lambda_rght = 3*p + i;
    
    curl_lambda_bndry(surf_id,lambda_left) = +1;
    curl_lambda_bndry(surf_id,lambda_rght) = -1;
end

%---------- right edges

for j = 2:p
    i = p+1;
    
    surf_id = (j-1)*(p+1) + i;
    
    lambda_up = p + j;
    lambda_dn = p + j-1;
    
    curl_lambda_bndry(surf_id,lambda_up) = -1;
    curl_lambda_bndry(surf_id,lambda_dn) = +1;
end
    
%---------- left edges

for j = 2:p
    i = 1;
    
    surf_id = (j-1)*(p+1) + i;
    
    lambda_up = 2*p + j;
    lambda_dn = 2*p + j-1;
    
    curl_lambda_bndry(surf_id,lambda_up) = -1;
    curl_lambda_bndry(surf_id,lambda_dn) = +1;
end

%---------- bottom-left corner

surf_id = 1;

lambda_rght = 1;
lambda_up = 2*p + 1;

curl_lambda_bndry(surf_id,lambda_rght) = -1;
curl_lambda_bndry(surf_id,lambda_up) = +1;

%---------- bottom-rght corner

surf_id = p+1;

lambda_left = p;
lambda_up = p+1;

curl_lambda_bndry(surf_id,lambda_up) = -1;
curl_lambda_bndry(surf_id,lambda_left) = +1;

%---------- topp-left corner

surf_id = p*(p+1) +1;

lambda_bttm = 3*p;
lambda_rght = 3*p +1;

curl_lambda_bndry(surf_id,lambda_bttm) = -1;
curl_lambda_bndry(surf_id,lambda_rght) = +1;

%---------- topp-rght corner

surf_id = (p+1)^2;

lambda_bttm = 2*p;
lambda_topp = 4*p;

curl_lambda_bndry(surf_id,lambda_topp) = -1;
curl_lambda_bndry(surf_id,lambda_bttm) = +1;

end

