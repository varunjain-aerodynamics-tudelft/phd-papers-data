function output_naca = get_y_naca_alpha(x, side, derivative, alpha);
            a = [0.2969 -0.126 -0.3516 0.2843 -0.1036];
            exp = [0.5 1 2 3 4];
            t2c = 0.18;

            for j = 1:5
                upper_collect(:,j) = a(j).*x.^(exp(j));
                d_upper_collect(:,j) = exp(j).*a(j).*x.^(exp(j)-1);
                lower_collect(:,j) = a(j).*x.^(exp(j));
                d_lower_collect(:,j) = exp(j).*a(j).*x.^(exp(j)-1);
            end
            if derivative == 0
                y_upper = -x.'.*sind(alpha)+sum(upper_collect,2).*t2c/0.2.*cosd(alpha);
                y_lower = x.'.*sind(alpha)+sum(lower_collect,2).*t2c/0.2.*cosd(alpha);     
            else
                y_upper = -sind(alpha)+sum(d_upper_collect,2).*t2c/0.2.*cosd(alpha);
                y_upper(y_upper==Inf)=100;
                y_lower = sind(alpha)+sum(d_lower_collect,2).*t2c/0.2.*cosd(alpha);
                y_lower(y_lower==-Inf)=-100;
            end
            switch side
                case 'up'
                    output_naca = y_upper;
                case 'lo'
                    output_naca = -y_lower;
            end
            % end of program
        end