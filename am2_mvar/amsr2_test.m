%
% beam_opt = 1: Pencil beam, multiple sensor_los, single mlock_dlos
% beam_opt = 2: Pencil beam, single sensor_los, multiple mlock_dlos
% beam_opt = 3: Antenna included, multiple sensor_los, single mlock_dlos
% beam_opt = 4: Antenna included, single sensor_los, multiple antenna_dlos
%
% FORMAT [y,y_geo] = amsr2_test(nfoot,angfac,ifreq,pol,beam_opt)
%
% OUT   y       As the ARTS WSV
%       y_geo   As the ARTS WSV
% IN    nfoot   Size of footptint scene to simulate. This number gives both
%               of many scans that will be done. and the number of footprints 
%               per scan that will be considered. T
%       angfac  An integer that control how the angular resoluion of the
%               simulations. The value 1 means to use some hard-coded
%               values. A value of 2 means that the angular spacing is
%               halfed, etc.
%       ifreq   Index of frequency to plot.
%       pol     Polarisation to plot, 'v' or 'h'.
%       beamopt See above.

function [y,y_geo] = amsr2_test(nfoot,angfac,ifreq,pol,beam_opt)
%
if nargin < 5, beam_opt = 1; end

assert( iswhole(nfoot) );
assert( iswhole(angfac) );


wfolder = '~/WORKAREA';


% Lats and lons
%
lat0     = 10;
lon0     = 0;
%
lat_grid = lat0 + [-3:0.25:3]';
lon_grid = lon0 + [-3:0.25:3]';
%
nlat     = length(lat_grid);
nlon     = length(lon_grid);
%
xmlStore( fullfile(wfolder,'lat_grid.xml'), lat_grid, 'Vector' );
xmlStore( fullfile(wfolder,'lon_grid.xml'), lon_grid, 'Vector' );


% Surface properties
%
surface_props_names{1} = 'Water skin temperature';
surface_props_names{2} = 'Wind speed';
surface_props_names{3} = 'Wind direction';
surface_props_names{4} = 'Salinity';
%
surface_props_data     = zeros( length(surface_props_names), ...
                                length(lat_grid), length(lon_grid) );
%
tmean                     = 296;
surface_props_data(1,:,:) = tmean;
surface_props_data(1,14,14) = tmean+5;
surface_props_data(2,:,:) = 5;
surface_props_data(2,12,13) = 10;
surface_props_data(3,:,:) = 45;
surface_props_data(4,:,:) = 0.03;
%
xmlStore( fullfile(wfolder,'surface_props_names.xml'), surface_props_names, ...
          'ArrayOfString' );
xmlStore( fullfile(wfolder,'surface_props_data.xml'), surface_props_data, ...
          'Tensor3' );
xmlStore( fullfile(wfolder,'z_surface.xml'), zeros(nlat,nlon), 'Matrix' );


% Atmosphere
%
z      = [ 0 : 500 : 10e3 ]';
%
p_grid    = z2p_simple( z );
t_field   = repmat( tmean-7e-3*z, [1,nlat,nlon] );
vmr_field = zeros( 4, length(z), nlat, nlon );
%
vmr_field(1,:,:,:) = 0.7808;
vmr_field(2,:,:,:) = 0.2095;
vmr_field(3,:,:,:) = relhum_to_vmr( 20, t_field, repmat(p_grid,[1,nlat,nlon]) );
vmr_field(4,:,:,:) = 0;
%
xmlStore( fullfile(wfolder,'p_grid.xml'), p_grid, 'Vector' );
xmlStore( fullfile(wfolder,'t_field.xml'), t_field, 'Tensor3' );
xmlStore( fullfile(wfolder,'z_field.xml'), repmat(z,[1,nlat,nlon]), 'Tensor3' );
xmlStore( fullfile(wfolder,'vmr_field.xml'), vmr_field, 'Tensor4' );


% Set and calculate some basic variables for antenna and scanning
%
zsat   = 705e3;       % Satellite altitde
vsat   = 6.76e3;      % Satellite velocity
dt     = 2.6e-3;      % Integration time
rpm    = 40;          % Rotations per minute
awidth = 2;           % Max half-width of antenna  
%
m2deg  = 1/(60*constants('NAUTICAL_MILE'));   % Conversion from m to latitude
dang   = dt * 360 * rpm / 60;                 % Angular distance between samples
dlat   = m2deg * vsat * 60 / rpm;             % Latitude distance between scans


% Antenna and f_grid
%
R      = amsr2_antenna( dang/(angfac*5), awidth );
f_grid = R.grids{2};
%
xmlStore( fullfile(wfolder,'f_grid.xml'), f_grid, 'Vector' );


% Determine bore-sight angles to use for one scan 
%
steps   = [-(nfoot-1)/2:(nfoot-1)/2]';
los0    = [ 180-47.5, 0 ];
bsights = [ repmat(los0(1),nfoot,1), los0(2)+dang*steps ];
  

% Store data describing angles
%
if beam_opt >= 3
  % Zenith grid matching mblock_dlos_grid 
  n       = floor( angfac * awidth / dang );
  za_grid = los0(1) + (dang/angfac) * ( -n:n );
end
%
switch beam_opt
  
  case 1
    % Here we set sensor_los to each bore-sight, and selects data to get
    % mblock_dlos_grid = (0,0)
    xmlStore( fullfile(wfolder,'sensor_los.xml'), repmat(bsights,nfoot,1), 'Matrix' );
    xmlStore( fullfile(wfolder,'mblock_reference_los.xml'), los0, 'Vector' );
    xmlStore( fullfile(wfolder,'mblock_target_los.xml'), los0, 'Matrix' );

  case 2
    % Here we set sensor_los to bore-sight of centre footprint, and simulate
    % scanning by mblock_dlos_grid
    xmlStore( fullfile(wfolder,'sensor_los.xml'), repmat( los0, nfoot, 1 ), 'Matrix' );
    xmlStore( fullfile(wfolder,'mblock_reference_los.xml'), los0, 'Vector' );
    xmlStore( fullfile(wfolder,'mblock_target_los.xml'), bsights, 'Matrix' );

  case 3
    % Here we set sensor_los to each bore-sight and mblock_dlos_grid to cover
    % just one antenna diagram and selects data to get antena_dlos = (0,0)
    xmlStore( fullfile(wfolder,'sensor_los.xml'), repmat(bsights,nfoot,1), 'Matrix' );
    xmlStore( fullfile(wfolder,'antenna_response.xml'), R, 'GriddedField4' );
    xmlStore( fullfile(wfolder,'antenna_reference_los.xml'), los0, 'Vector' );
    xmlStore( fullfile(wfolder,'antenna_target_los.xml'), los0, 'Matrix' );
    xmlStore( fullfile(wfolder,'mblock_reference_los.xml'), los0, 'Vector' );
    aa_grid = los0(2) + (dang/angfac) * ( -n:n );
    xmlStore( fullfile(wfolder,'mblock_target_los.xml'), ...
              [repmat(za_grid',length(aa_grid),1), ...
               mat2col(repmat(aa_grid,length(za_grid),1))], 'Matrix' );
    
  case 4
    % Here we set sensor_los to bore-sight of centre footprint, and simulate
    % scanning by antenna_dlos. mblock_dlos_grid must cover the full scan
    xmlStore( fullfile(wfolder,'sensor_los.xml'), repmat( los0, nfoot, 1 ), 'Matrix' );
    xmlStore( fullfile(wfolder,'antenna_response.xml'), R, 'GriddedField4' );
    xmlStore( fullfile(wfolder,'antenna_reference_los.xml'), los0, 'Vector' );
    xmlStore( fullfile(wfolder,'antenna_target_los.xml'), bsights, 'Matrix' );
    xmlStore( fullfile(wfolder,'mblock_reference_los.xml'), los0, 'Vector' );
    n       = 1 + floor( angfac * ( bsights(end,2) + awidth ) / dang );
    aa_grid = los0(2) + (dang/angfac) * ( -n:n );
    xmlStore( fullfile(wfolder,'mblock_target_los.xml'), ...
              [repmat(za_grid',length(aa_grid),1), ...
               mat2col(repmat(aa_grid,length(za_grid),1))], 'Matrix' );

  otherwise
    error( 'Unknown choice for *beam_opt*.' );
end


% Set sensor_pos
%
shift   = 7.5;   % Latitude shift to centre calculations around lat0
%
if isodd( beam_opt )
  spos = [ repmat(zsat,nfoot*nfoot,1), ...
           mat2col( repmat(lat0-shift+dlat*steps, 1, nfoot )' ), ... 
           repmat(lon0,nfoot*nfoot,1) ]; 
else
  spos = [ repmat(zsat,nfoot,1), lat0-shift+dlat*steps, repmat(lon0,nfoot,1) ];
end
%
xmlStore( fullfile(wfolder,'sensor_pos.xml'), spos, 'Matrix' );




% Run ARTS
%
copyfile( 'amsr2_cfile_common.arts', fullfile(wfolder,'cfile_common.arts' ) );
%
cfile = fullfile(wfolder,'cfile.arts' );
if beam_opt <= 2
  copyfile( 'amsr2_cfile_pbeam.arts', cfile );
else
  copyfile( 'amsr2_cfile.arts', cfile );
end
%
[s,r] = system( sprintf( 'arts -r100 -o %s %s', wfolder, cfile ) );
%
if s
  disp(r)
  return
end
%
y = xmlLoad( fullfile(wfolder,'y.xml' ) );
y_geo = xmlLoad( fullfile(wfolder,'y_geo.xml' ) );



% Plot
%
if pol == 'v'
  tb  = y(1:2:end) + y(2:2:end);
elseif pol == 'h'
  tb  = y(1:2:end) - y(2:2:end);
else
  error( '*pol* has to be ''v'' or ''h''.' ); 
end
geo = y_geo(1:2:end,:);
%
nf     = length(f_grid);
ind    = ifreq:nf:length(tb);
%
scatter( geo(ind,3), geo(ind,2), 80, tb(ind), 'filled' );
colorbar
xlabel( 'Longitude' );
ylabel( 'Latitude' );
title( sprintf( '%.1f GHz / %s-pol', f_grid(ifreq)/1e9, upper(pol) ) );

