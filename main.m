close all
clear all

load aa3_dr.mat
load aa3_gpsx.mat
load aa3_lsr2.mat

do_plot = true;
do_plot_ellipse = false;
do_plot_lidar = false;
do_cal = false;

if do_cal
    pkg load statistics
end

ndata = size(time, 1) + size(timeGps, 1) + size(TLsr, 1);
manager = zeros(ndata, 4);
manager(1:size(time, 1), 1) = time;
manager(1:size(time, 1), 2) = 0;
manager(1:size(time, 1), 3) = speed;
manager(1:size(time, 1), 4) = steering;
manager(size(time, 1)+1:size(time, 1)+size(timeGps, 1), 1) = timeGps;
manager(size(time, 1)+1:size(time, 1)+size(timeGps, 1), 2) = 1;
manager(size(time, 1)+1:size(time, 1)+size(timeGps, 1), 3) = Lo_m;
manager(size(time, 1)+1:size(time, 1)+size(timeGps, 1), 4) = La_m;
manager(size(time, 1)+size(timeGps, 1)+1:end, 1) = TLsr;
manager(size(time, 1)+size(timeGps, 1)+1:end, 2) = 2;
manager(size(time, 1)+size(timeGps, 1)+1:end, 3) = 1:size(TLsr, 1);

manager = sortrows(manager);
timeline = manager(:, 1);

ekf = VictoriaEKF([Lo_m(1); La_m(1); 36*pi/180], diag([0.1^2, 0.1^2, 1^2]));

if do_plot
    hold on;
    grid on;
    box on;
end

last_t = timeGps(1);
gps_idx = 1;
for i = 1:ndata
    if ~mod(i, 100)
        disp([num2str(i) '/' num2str(ndata)])
    end
    t = manager(i, 1);
    switch round(manager(i, 2))
        case 0 % odo
            ekf.time_update(manager(i, 3:4)', (t - last_t)/1000);
        case 1 % gps
            if ~do_cal
                ekf.gps_update(manager(i, 3:4)');
                if do_plot
                    h_gps = plot(Lo_m(gps_idx), La_m(gps_idx), 
                                'r.', 'MarkerSize', 10);
                    gps_idx = gps_idx + 1;
                end
            end
        case 2 % lidar
            lidars = detect(LASER(round(manager(i, 3)), :));
            assoc = ekf.data_association(lidars);
            ekf.lidar_update(lidars, assoc);
            ekf.init_trees(lidars, assoc);
            if do_plot
                if exist('h_tree')
                    delete(h_tree);
                end
                if exist('h_ellipse')
                    for h = h_ellipse
                        if isgraphics(h)
                            delete(h);
                        end
                    end
                end
                if exist('h_line1')
                    for h = h_line1
                        if isgraphics(h)
                            delete(h);
                        end
                    end
                end
                if exist('h_line2')
                    for h = h_line2
                        if isgraphics(h)
                            delete(h);
                        end
                    end
                end
                
                h_tree = plot(ekf.get_trees_x(), ekf.get_trees_y(), 'gx', 
                            'MarkerSize', 10, 'linewidth', 3);
                
                if do_plot_ellipse
                    for idx = 1:ekf.num_trees
                        covariance = ekf.get_tree_cov(idx);
                        pos = [ekf.get_tree_x(idx), ekf.get_tree_y(idx)];
                        r_ellipse = ekf.ellipse(covariance, pos);
                        h_ellipse(idx) = plot(r_ellipse(:, 1), r_ellipse(:, 2), 'b');
                    end
                    covariance = ekf.get_vehicle_cov();
                    pos = [ekf.get_vehicle_x(), ekf.get_vehicle_y()];
                    r_ellipse = ekf.ellipse(covariance, pos);
                    h_ellipse(ekf.num_trees+1) = plot(r_ellipse(:, 1), r_ellipse(:, 2), 'k');
                end
                
                if do_plot_lidar
                    for j = 1:size(ekf.matched_trees, 2)
                        idx = ekf.matched_trees(j);
                        h_line1(j) = plot([ekf.get_vehicle_x() ekf.get_tree_x(idx)],
                            [ekf.get_vehicle_y() ekf.get_tree_y(idx)], 'g');
                    end
                    for j = 1:size(ekf.new_trees, 2)
                        idx = ekf.new_trees(j);
                        h_line2(j) = plot([ekf.get_vehicle_x() ekf.get_tree_x(idx)],
                            [ekf.get_vehicle_y() ekf.get_tree_y(idx)], 'b');
                    end
                end
            end
    end
    last_t = t;
    
    if do_plot    
        plot(ekf.get_vehicle_x(), ekf.get_vehicle_y(), 'k.', 'MarkerSize', 1);
        if ~do_plot_ellipse
            if exist('h_car')
                delete(h_car);
            end
            h_car = plot(ekf.get_vehicle_x(), ekf.get_vehicle_y(), 
                'ko', 'MarkerSize', 10, 'linewidth', 2);
        end
        title(['t = ' num2str(t)]);
        legend('GPS');
        drawnow;
    end
    
    if do_cal
        vehicle_x(i)= ekf.get_vehicle_x();
        vehicle_y(i)= ekf.get_vehicle_y();
        if (t == timeGps(gps_idx))
            err(gps_idx) = pdist([vehicle_x(i), vehicle_y(i); 
                Lo_m(gps_idx), La_m(gps_idx)], 'euclidean');
            gps_idx = gps_idx + 1;
            occurance = 1:gps_idx-1;
            plot(occurance, err);
            grid on
            ylabel("Euclidean Distance");
            xlabel("Occurence of GPS signal");
            title("Distance between GPs coor. and estimated pose");
            drawnow;
        end
    end
end

pause;
