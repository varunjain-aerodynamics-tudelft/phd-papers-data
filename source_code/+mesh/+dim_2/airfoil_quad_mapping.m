function [domain] = airfoil_quad_mapping(x1,x2,x3,x4,y1,y2,y3,y4,typeL,typeR,typeB,typeT,alpha)

% --- substitution below
% lineqn = @(s,z1,z2) s*(z2-z1) + z1;
% dlineqn_ds = @(s,z1,z2) ones(size(s)).*(z2-z1);

mapp.F_L_x = @(s,t) get_x_airfoil_alpha(t,x1,x4,typeL,alpha);
mapp.F_L_y = @(s,t) get_y_airfoil_alpha(t,y1,y4,x1,x4,typeL,alpha);

mapp.F_R_x = @(s,t) get_x_airfoil_alpha(t,x2,x3,typeR,alpha);
mapp.F_R_y = @(s,t) get_y_airfoil_alpha(t,y2,y3,x2,x3,typeR,alpha);

mapp.F_B_x = @(s,t) get_x_airfoil_alpha(s,x1,x2,typeB,alpha);
mapp.F_B_y = @(s,t) get_y_airfoil_alpha(s,y1,y2,x1,x2,typeB,alpha);

mapp.F_T_x = @(s,t) get_x_airfoil_alpha(s,x4,x3,typeT,alpha);
mapp.F_T_y = @(s,t) get_y_airfoil_alpha(s,y4,y3,x4,x3,typeT,alpha);

derr.dF_L_x = @(s,t) get_dx_airfoil_alpha(t,x1,x4,typeL,alpha);
derr.dF_L_y = @(s,t) get_dy_airfoil_alpha(t,y1,y4,x1,x4,typeL,alpha);

derr.dF_R_x = @(s,t) get_dx_airfoil_alpha(t,x2,x3,typeR,alpha);
derr.dF_R_y = @(s,t) get_dy_airfoil_alpha(t,y2,y3,x2,x3,typeR,alpha);

derr.dF_B_x = @(s,t) get_dx_airfoil_alpha(s,x1,x2,typeB,alpha);
derr.dF_B_y = @(s,t) get_dy_airfoil_alpha(s,y1,y2,x1,x2,typeB,alpha);

derr.dF_T_x = @(s,t) get_dx_airfoil_alpha(s,x4,x3,typeT,alpha);
derr.dF_T_y = @(s,t) get_dy_airfoil_alpha(s,y4,y3,x4,x3,typeT,alpha);

domain.mapping = @(xii,eta) mesh.transfinite_mesh.mapping(xii,eta,mapp);
domain.jacobian = @(xii,eta) mesh.transfinite_mesh.jacobian(xii,eta,mapp,derr);
    
    %% --- substitution of 'lineqn' in x
    function lineqnx = get_x_airfoil_alpha(s,z1,z2,type,alpha)
        switch type
            case 'li'
                lineqnx = s*(z2-z1) + z1;
            case 'af_up'
                lineqnx = get_x_naca_alpha((s*(z2-z1)+z1), 'up', 0, alpha)';
            case 'af_lo'
                lineqnx = get_x_naca_alpha((s*(z2-z1)+z1), 'lo', 0, alpha)';
        end
        % --- function for x naca
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
    end
    
    %% --- substitution of 'lineqn' in y
    function lineqny = get_y_airfoil_alpha(s,z1,z2,ref1,ref2,type,alpha)
        switch type
            case 'li'
                lineqny = s*(z2-z1) + z1;
            case 'af_up'
                lineqny = get_y_naca_alpha((s*(ref2-ref1)+ref1), 'up', 0, alpha)';
            case 'af_lo'
                lineqny = get_y_naca_alpha((s*(ref2-ref1)+ref1), 'lo', 0, alpha)';
        end
        % --- function for y naca
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
    end
    
    %% --- substitution of 'dlineqn' in x
    function dlineqnx = get_dx_airfoil_alpha(s,z1,z2,type,alpha)
        switch type
            case 'li'
                dlineqnx = (z2-z1).*ones(1,length(s));
            case 'af_up'
                dlineqnx = (z2-z1).*get_x_naca_alpha((s*(z2-z1)+z1), 'up', 1, alpha)';
            case 'af_lo'
                dlineqnx = (z2-z1).*get_x_naca_alpha((s*(z2-z1)+z1), 'lo', 1, alpha)';
        end
        % --- function for dx naca
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
    end

    %% --- substitution of 'dlineqn' in y
    function dlineqny = get_dy_airfoil_alpha(s,z1,z2,ref1,ref2,type,alpha)
        switch type
            case 'li'
                dlineqny = (z2-z1).*ones(1,length(s));
            case 'af_up'
                dlineqny = (ref2-ref1).*get_y_naca_alpha((s*(ref2-ref1)+ref1), 'up', 1, alpha)';
            case 'af_lo'
                dlineqny = (ref2-ref1).*get_y_naca_alpha((s*(ref2-ref1)+ref1), 'lo', 1, alpha)';
        end
        % --- function for dy naca
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
    end

end