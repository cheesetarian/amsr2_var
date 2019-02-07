function R = amsr2_antenna(resol,width)

data = [  6.925e9 1.8
          7.3e9   1.8
         10.96e9  1.2
         18.7e9   0.65
         23.8e9   0.75
         36.5e9   0.35
         89.0e9   0.15 ];

x  = symgrid( 0:resol:width );
x2 = x.^2;
nf = size(data,1);
nx = length(x);

R.name      = 'AMSR2 antenna response';
R.gridnames = { 'Polarisation', 'Frequency', 'Zenith angle', 'Azimuth angle' };
R.grids     = { {'1'}, data(:,1), x, x };
R.dataname  = 'Response';
R.data      = zeros( 1, nf, nx, nx );

for i = 1:nf

  si = fwhm2si( data(i,2) );

  R.data(1,i,:,:) = exp( -(repmat(x2,nx,1)+repmat(x2',1,nx))/si^2 );
    
end

