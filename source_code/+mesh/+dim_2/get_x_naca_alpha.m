function output_naca = get_x_naca_alpha(x, side, derivative, alpha);
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
                x_upper = x.'.*cosd(alpha)+sum(upper_collect,2).*t2c/0.2.*sind(alpha);
                x_lower = x.'.*cosd(alpha)-sum(lower_collect,2).*t2c/0.2.*sind(alpha);     
            else
                x_upper = cosd(alpha)+sum(d_upper_collect,2).*t2c/0.2.*cosd(alpha);
                x_upper(x_upper==Inf)=100;
                x_lower = cosd(alpha)-sum(d_lower_collect,2).*t2c/0.2.*cosd(alpha);
                x_lower(x_lower==-Inf)=-100;
            end
            switch side
                case 'up'
                    output_naca = x_upper;
                case 'lo'
                    output_naca = x_lower;
            end
            % end of program
        end