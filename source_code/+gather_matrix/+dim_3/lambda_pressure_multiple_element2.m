function [GM] = lambda_pressure_multiple_element2(K1x,K1y,K1z,K2x,K2y,K2z,p)

ttl_nr_el = K1x*K1y*K1z;

nr_boundary_faces = 2*K2x*K2y*p^2 + 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2;

GM = zeros(ttl_nr_el,nr_boundary_faces);

Kx = K1x;
Ky = K1y;
Kz = K1z;

% right faces

count1 = 1;

for k = 1:Kz
    for j = 1:Ky
        for i = 1:Kx-1
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
            faceid = 1:K2y*K2z*p^2;
            GM(eleeid,faceid) = (count1-1)*K2y*K2z*p^2 + 1 : (count1-1)*K2y*K2z*p^2 + K2y*K2z*p^2;
            count1 = count1+1;
        end
    end
end

% left faces

for i = 2:Kx
    for j = 1:Ky
        for k = 1:Kz
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
%             faceid = K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2+1: 2*K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2;
            faceid = K2y*K2z*p^2 + 1:2*K2y*K2z*p^2;
            facert = 1:K2y*K2z*p^2;
            GM(eleeid,faceid) = GM(eleeid-1,facert);
        end
    end
end

% top faces

count2 = count1;

for k = 1:Kz
    for j = 1:Ky-1
        for i = 1:Kx
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
%             faceid = K2y*K2z*p^2+1:K2y*K2z*p^2 + K2x*K2z*p^2;
            faceid = 2*K2y*K2z*p^2 + 1 : 2*K2y*K2z*p^2 + K2x*K2z*p^2;
            GM(eleeid,faceid) = (count2-1)*K2x*K2z*p^2 + 1 : (count2-1)*K2x*K2z*p^2 + K2x*K2z*p^2;
            count2 = count2+1;
        end
    end
end

% bottom faces

for i = 1:Kx
    for j = 2:Ky
        for k = 1:Kz
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
%             faceid = 2*K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2+1:2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2;
            faceid = 2*K2y*K2z*p^2 + K2x*K2z*p^2 + 1 : 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2;
%             facetp = K2y*K2z*p^2+1:K2y*K2z*p^2 + K2x*K2z*p^2;
            facetp = 2*K2y*K2z*p^2 + 1 : 2*K2y*K2z*p^2 + K2x*K2z*p^2;
            GM(eleeid,faceid) = GM(eleeid-Kx,facetp);
        end
    end
end

count3 = count2;

% front faces

for k = 1:Kz-1
    for j = 1:Ky
        for i = 1:Kx
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
%             faceid = K2y*K2z*p^2 + K2x*K2z*p^2 +1:K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2;
            faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 +1 : 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2;
            GM(eleeid,faceid) = (count3-1)*K2x*K2y*p^2 + 1 : (count3-1)*K2x*K2y*p^2 + K2x*K2y*p^2;
            count3 = count3+1;
        end
    end
end

% back faces

for i = 1:Kx
    for j = 1:Ky
        for k = 2:Kz
            eleeid = (k-1)*Kx*Ky + (j-1)*Kx + i;
%             faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2+1:2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + 2*K2x*K2y*p^2;
            faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2 + 1: 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + 2*K2x*K2y*p^2;
%             faceft = K2y*K2z*p^2 + K2x*K2z*p^2 +1:K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2;
            faceft = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 +1 : 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2;
            GM(eleeid,faceid) = GM(eleeid-Kx*Ky,faceft);
        end
    end
end

% Right domain faces

count4 = count3;

for k=1:Kz
    for j=1:Ky
        eleeid = (k-1)*Kx*Ky + (j-1)*Kx + Kx;
        faceid = 1:K2y*K2z*p^2;
        GM(eleeid,faceid) = (count4-1)*K2y*K2z*p^2 + 1 : (count4-1)*K2y*K2z*p^2 + K2y*K2z*p^2;
        count4 = count4 + 1;
    end
end

% Left domain faces

count5 = count4;

for k=1:Kz
    for j=1:Ky
        eleeid = (k-1)*Kx*Ky + (j-1)*Kx + 1;
%         faceid = K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2+1: 2*K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2;
        faceid = K2y*K2z*p^2 + 1:2*K2y*K2z*p^2;
        GM(eleeid,faceid) = (count5-1)*K2y*K2z*p^2 + 1 : (count5-1)*K2y*K2z*p^2 + K2y*K2z*p^2;
        count5 = count5 + 1;
    end
end

% Top domain faces

count6 = count5;

for k=1:Kz
    for i=1:Kx
        eleeid = (k-1)*Kx*Ky + (Ky-1)*Kx + i;
%         faceid = K2y*K2z*p^2+1:K2y*K2z*p^2 + K2x*K2z*p^2;
        faceid = 2*K2y*K2z*p^2 + 1 : 2*K2y*K2z*p^2 + K2x*K2z*p^2;
        GM(eleeid,faceid) = (count6-1)*K2x*K2z*p^2 + 1 : (count6-1)*K2x*K2z*p^2 + K2x*K2z*p^2;
        count6 = count6 + 1;
    end
end

% Bottom domain faces

count7 = count6;

for k=1:Kz
    for i=1:Kx
        eleeid = (k-1)*Kx*Ky + (1-1)*Kx + i;
%         faceid = 2*K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2+1:2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2;
        faceid = 2*K2y*K2z*p^2 + K2x*K2z*p^2 + 1 : 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2;
        GM(eleeid,faceid) = (count7-1)*K2x*K2z*p^2 + 1 : (count7-1)*K2x*K2z*p^2 + K2x*K2z*p^2;
        count7 = count7 + 1;
    end
end

% Front domain faces

count8 = count7;

for j=1:Ky
    for i=1:Kx
        eleeid = (Kz-1)*Kx*Ky + (j-1)*Kx + i;
%         faceid = K2y*K2z*p^2 + K2x*K2z*p^2 +1:K2y*K2z*p^2 + K2x*K2z*p^2 + K2x*K2y*p^2;
        faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 +1 : 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2;
        GM(eleeid,faceid) = (count8-1)*K2x*K2y*p^2 + 1 : (count8-1)*K2x*K2y*p^2 + K2x*K2y*p^2;
        count8 = count8 + 1;
    end
end

% Back domain faces

count9 = count8;

for j=1:Ky
    for i=1:Kx
        eleeid = (1-1)*Kx*Ky + (j-1)*Kx + i;
%         faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2+1:2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + 2*K2x*K2y*p^2;
        faceid = 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + K2x*K2y*p^2 + 1: 2*K2y*K2z*p^2 + 2*K2x*K2z*p^2 + 2*K2x*K2y*p^2;
        GM(eleeid,faceid) = (count9-1)*K2x*K2y*p^2 + 1 : (count9-1)*K2x*K2y*p^2 + K2x*K2y*p^2;
        count9 = count9 + 1;
    end
end

GM = GM';

end

