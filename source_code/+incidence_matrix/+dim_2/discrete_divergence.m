function [E21] = discrete_divergence(p)

nr_surfc = 2*p*(p+1);
nr_volum = p^2;

E21 = zeros(nr_volum,nr_surfc);

for i = 1:p
    for j = 1:p
        vol_id = (j-1)*p + i;
        
        face_bttm = (j-1)*p + i;
        face_topp = j*p + i;
        face_left = p*(p+1) + (j-1)*(p+1) + i;
        face_rght = p*(p+1) + (j-1)*(p+1) + i + 1;
        
        E21(vol_id,face_bttm) = -1;
        E21(vol_id,face_topp) = +1;
        E21(vol_id,face_left) = -1;
        E21(vol_id,face_rght) = +1;
    end
end

E21 = sparse(E21);

end

