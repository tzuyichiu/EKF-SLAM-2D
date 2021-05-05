classdef VictoriaEKF < handle
    properties
        % vehicle constant parameters
        a = 3.78;
        b = 0.50;
        L = 2.83;
        H = 0.76;
        
        % errors
        % sigma_x = sigma_y = 5 cm, sigma_phi = 0.5 degree
        sigma_x = 0.05;
        sigma_y = 0.05;
        sigma_phi = 0.5*pi/180;
        
        % sigma_gps = 3 m
        sigma_gps = 3;

        % sigma_beta = 0.5 m, sigma_gamma = 5 degrees
        sigma_beta = 0.5;
        sigma_gamma = 5*pi/180;

        % others
        num_trees = 0;
        state_length = 3;
        trees_diameter = [];
        new_trees
        matched_trees

        % EKF
        cov_odo
        cov_gps
        cov_lidar
        state
        cov
    end

    methods
        function self = VictoriaEKF(init_state, init_cov)
            self.state = init_state;
            self.cov = init_cov;
            self.cov_odo = diag([self.sigma_x^2, self.sigma_y^2, self.sigma_phi^2]);
            self.cov_gps = (self.sigma_gps^2)*eye(2);
            self.cov_lidar = diag([self.sigma_beta^2, self.sigma_gamma^2]);
        end

        %%% get methods to access properties %%%
        function vehicle_x = get_vehicle_x(self)
            vehicle_x = self.state(1);
        end
        
        function vehicle_y = get_vehicle_y(self)
            vehicle_y = self.state(2);
        end

        function vehicle_phi = get_vehicle_phi(self)
            vehicle_phi = self.state(3);
        end

        function vehicle_cov = get_vehicle_cov(self)
            vehicle_cov = self.cov(1:2, 1:2);
        end

        function tree_x = get_tree_x(self, idx)
            tree_x = self.state(2*idx+2);
        end

        function tree_y = get_tree_y(self, idx)
            tree_y = self.state(2*idx+3);
        end

        function trees_x = get_trees_x(self)
            trees_x = self.state(4:2:end);
        end

        function trees_y = get_trees_y(self)
            trees_y = self.state(5:2:end);
        end

        function tree_cov = get_tree_cov(self, idx)
            tree_idx = 2*idx+2:2*idx+3;
            tree_cov = self.cov(tree_idx, tree_idx);
        end
        %%%

        function self = clamp(self)
            % make within [-pi, pi[
            self.state(3) = mod(self.state(3) + pi, 2*pi) - pi;
            
            % make symetric
            self.cov = (self.cov + self.cov')/2;
        end

        function self = time_update(self, odo, dt) % odo --- [ve; alpha]
            ve = odo(1);
            alp = odo(2);
            phi = self.state(3);
            c = cos(phi);
            s = sin(phi);
            ta = tan(alp);
            vc = ve / (1 - ta*self.H/self.L);
            
            % state update
            dxdt = [vc*(c - ta*(self.a*s + self.b*c)/self.L);
                    vc*(s + ta*(self.a*c - self.b*s)/self.L);
                    vc*ta/self.L];
            self.state(1:3) = self.state(1:3) + dt*dxdt;

            % cov update
            dfdx = eye(self.state_length);
            dfdx(1, 3) = -dt*vc*(s + ta*(self.a*c - self.b*s)/self.L);
            dfdx(2, 3) = dt*vc*(c - ta*(self.a*s + self.b*c)/self.L);
            
            self.cov = dfdx*self.cov*(dfdx');
            self.cov(1:3, 1:3) = self.cov(1:3, 1:3) + self.cov_odo;
            
            % clamp
            self.clamp();
        end

        function self = gps_update(self, gps) % gps --- [x; y]
            dhdx = zeros(2, self.state_length);
            dhdx(1:2, 1:2) = eye(2);
            phi = self.state(3);
            c = cos(phi);
            s = sin(phi);
            dhdx(1, 3) = self.a*s + self.b*c;
            dhdx(2, 3) = -self.a*c + self.b*s;
            tmp = self.cov*dhdx';
            S = dhdx*tmp + self.cov_gps;
            K = tmp/S;
            residual = gps - self.state(1:2);

            % throw out unlikely measurements
            mahalanobis = residual'*(S\residual);
            if mahalanobis < 13.8
                self.state = self.state + K*residual;
                self.cov = (eye(self.state_length) - K*dhdx)*self.cov;
                self.clamp();
            else
                disp('Unlikely GPS measurement')
            end
        end
        
        function self = init_trees(self, lidars, assoc) 
            % lidars --- [[beta; gamma; diameter], ...]
            new_detections = lidars(:, assoc == 0);
            for lidar = new_detections
                self.trees_diameter(end+1) = lidar(3);
                delta = lidar(1);
                bearing = lidar(2) + self.state(3);
                c = cos(bearing);
                s = sin(bearing);
                tree_state = delta.*[c; s] + self.state(1:2);
                W = [c, -delta*s; s, delta*c];
                tree_cov = W*self.cov_lidar*W' + self.cov(1:2, 1:2);
                self.state = [self.state; tree_state];
                self.cov = [self.cov, zeros(self.state_length, 2); 
                            zeros(2, self.state_length), tree_cov];
                self.num_trees = self.num_trees + 1;
                self.state_length = self.state_length + 2;
            end
            self.new_trees = (self.num_trees-size(new_detections, 2)+1):self.num_trees;
        end
        
        % data association
        function [mahalanobis, dhdx, residual] = cost(self, idx, lidar) 
            % lidar --- [beta; gamma; diameter]
            tree_idx = 2*idx+2:2*idx+3;
            tree_state = self.state(tree_idx);
            delta_x = tree_state(1) - self.state(1);
            delta_y = tree_state(2) - self.state(2);
            delta2 = delta_x^2 + delta_y^2;
            delta = sqrt(delta2);
            bearing = atan2(delta_y, delta_x) - self.state(3);
            bearing = mod(bearing + pi, 2*pi) - pi;

            dhdx_vehicle = [-delta_x/delta, -delta_y/delta, 0;
                            delta_y/delta2, -delta_x/delta2, -1]; % 2x3
            dhdx = zeros(2, self.state_length);
            dhdx(:, 1:3) = dhdx_vehicle;
            dhdx(:, tree_idx) = -dhdx_vehicle(:, 1:2);
            tmp = self.cov*dhdx'; % 5x2
            S = dhdx*tmp + self.cov_lidar;
            residual = lidar(1:2) - [delta; bearing];
            mahalanobis = residual'*(S\residual);
        end

        function assoc = data_association(self, lidars)
            % lidars --- [[beta; gamma; diameter], ...]
            % returns:
            %   assoc same length as lidars (measurements 2xn)
            %   assoc(j) = 0 if j is new detection; 
            %   assoc(j) = i if tree i corresponds to j;
            %   assoc(j) = -1 if j is unmatched (unlikely measurement)
            num_meas = size(lidars, 2);
            assoc = zeros(1, num_meas);
            cost_matrix = 5.99*ones(num_meas, self.num_trees + num_meas);
            for j = 1:num_meas
                for i = 1:self.num_trees
                    [mahalanobis, ~, ~] = self.cost(i, lidars(:, j));
                    cost_matrix(j, i) = mahalanobis;
                end
            end

            [assign, ~] = munkres(cost_matrix);
            for j = 1:num_meas
                i = find(assign(j, :));
                if i > self.num_trees
                    if min(cost_matrix(j, 1:self.num_trees)) < 13.81
                        assoc(j) = -1;
                    end
                elseif cost_matrix(j, i) < 5.99
                    assoc(j) = i;
                else
                    assoc(j) = -1;
                end
            end
        end
        
        function self = lidar_update(self, lidars, assoc) 
            % lidars --- [[beta; gamma; diameter], ...]
            self.matched_trees = find(assoc > 0);
            num_matched = sum(assoc > 0);
            dhdx = zeros(2*num_matched, self.state_length);
            residual = zeros(2*num_matched, 1);
            K = zeros(self.state_length, 2*num_matched);
            for idx = 1:num_matched
                matched_idx = self.matched_trees(idx);
                tree_idx = assoc(matched_idx);
                lidar = lidars(1:2, matched_idx);
                [~, tree_dhdx, tree_residual] = self.cost(tree_idx, lidar);
                dhdx(idx:idx+1, :) = tree_dhdx;
                tmp = self.cov*tree_dhdx';
                S = tree_dhdx*tmp + self.cov_lidar;
                K(:, idx:idx+1) = tmp/S;
                residual(idx:idx+1, :) = tree_residual;
            end
            
            self.state = self.state + K*residual;
            self.cov = (eye(self.state_length) - K*dhdx)*self.cov;
            self.clamp();
        end

        function r_ellipse = ellipse(self, covariance, pos)
            % plot ellipse from cov
            [eigenvec, eigenval] = eig(covariance);
            [largest_eigenvec_ind_c, rR] = find(eigenval == max(max(eigenval)));
            largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
            largest_eigenval = max(max(eigenval));
            if largest_eigenvec_ind_c == 1
                smallest_eigenval = max(eigenval(:,2));
                smallest_eigenvec = eigenvec(:,2);
            else
                smallest_eigenval = max(eigenval(:,1));
                smallest_eigenvec = eigenvec(1,:);
            end
            
            angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
            angle = mod(angle, 2*pi);
            
            theta_grid = linspace(0, 2*pi);
            phi = angle;
            chisquare_val = 5.99;
            a = sqrt(chisquare_val*largest_eigenval);
            b = sqrt(chisquare_val*smallest_eigenval);
            
            ellipse_x_r = a*cos(theta_grid);
            ellipse_y_r = b*sin(theta_grid);
            c = cos(phi);
            s = sin(phi);
            R = [c s; -s c];
            r_ellipse = [ellipse_x_r; ellipse_y_r]' * R' + pos;
        end
    end
end
