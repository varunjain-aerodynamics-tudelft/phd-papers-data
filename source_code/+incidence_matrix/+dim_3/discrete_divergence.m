function [E32] = discrete_divergence(p)

nr_surfc = 3*p^2*(p+1);
nr_volum = p^3;

E32 = zeros(nr_volum,nr_surfc);

for i = 1:p
    for j = 1:p
        for k =1:p
            vol_id = (k-1)*p^2 + (j-1)*p + i;
            face_yz_1 = (k-1)*p*(p+1) + (j-1)*(p+1) + i;
            face_yz_2 = (k-1)*p*(p+1) + (j-1)*(p+1) + i + 1;
            face_xz_1 = p^2*(p+1) + (k-1)*p*(p+1) + (j-1)*p + i;
            face_xz_2 = p^2*(p+1) + (k-1)*p*(p+1) + (j-1)*p + i + p;
            face_xy_1 = 2*p^2*(p+1) + (k-1)*p^2 + (j-1)*p + i;
            face_xy_2 = 2*p^2*(p+1) + (k-1)*p^2 + (j-1)*p + i + p^2;
            
            E32(vol_id,face_yz_1) = -1;
            E32(vol_id,face_yz_2) = +1;
            E32(vol_id,face_xz_1) = -1;
            E32(vol_id,face_xz_2) = +1;
            E32(vol_id,face_xy_1) = -1;
            E32(vol_id,face_xy_2) = +1;
        end
    end
end

E32 = sparse(E32);

end

