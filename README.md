# EKF SLAM (Odometry, GPS, LIDAR) on Victoria Park 2D Dataset

## Authors

Tzu-yi Chiu, Narimane Zennaki (Ecole Polytechnique de Montréal)

## Usage

The script to be executed is `main.m` with Matlab/Octave (both compatible).

One can manually set the boolean variables `do_plot`, `do_plot_ellipse` 
(covariance) and `do_plot_lidar` (lidar detections matched with landmarks, and 
initialization of landmarks) for different levels of visualization, at the cost 
of the running speed. 

The figure shows the evolution of the vehicle position and the map of trees 
being constructed, while also showing the received GPS signals.
`do_cal` is set to true if one wants to obtain the error evolution, comparing 
the vehicle position to the GPS data.

## Files

- `aa3_dr.mat`, `aa3_gpsx.mat` and `aa3_lsr2.mat` are the Victoria Park 
datasets, corresponding to odometry, GPS and LIDAR respectively.

- `info.txt` contains the format information relative to the dataset.

- `main.m` is the main script to be executed.

- `detectTreesI16.m` is provided along with the Victoria Park dataset, which
allows to transform LIDAR raw data to measurements. `detect.m` a 
post-processing function for our own purpose.

- `VictoriaEKF.m` is a class which implements the EKF.

- `munkres.m` is a function which implements the Hungarian algorithm, allowing 
to compute data association. Borrowed from: 
https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
