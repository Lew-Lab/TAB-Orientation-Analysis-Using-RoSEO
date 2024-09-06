% 190718 - v3 TD - Excluded Fourier transform calculation from this
% funcition for reducing calculation cost

% 190717 - TD - Fixed roi function and changed the bases generator to
% Nanoscope.computeBases from the output of Nanoscope.computeBasis. Changed
% the normalization factor brightness_scaling.

% 190716 - v2 TD - use Nanoscope object directly. Normalized images on
% camera before multiplying with signal s.

% 190404 Tianben Ding
% Calculate Fisher information matric of second moment
% Using fully pixelated images
% Ober, R.J., Ram, S., and Ward, E.S., “Localization accuracy in single-
% molecule microscopy.,” Biophysical journal 86(2), 1185–200 (2004).
% Shechtman, Y., Sahl, S.J., Backer, A.S., and Moerner, W.E., “Optimal 
% point spread function design for 3D imaging,” Physical Review Letters 
% 113(3), 1–5 (2014).

function FIM = calFIMSecondM_v3(Bx,By,Ix,Iy,s,backg)

backgx = backg( :,1:(size(backg,2)/2) );
backgy = backg( :,(size(backg,2)/2)+1:size(backg,2) );

FIM = zeros(6,6);

FIM(1,1) = 0.5.*(s^2).*(   sum(  sum( (Bx.aa)./(Ix + backgx) + (By.aa)./(Iy + backgy) )  )   );
FIM(2,2) = 0.5.*(s^2).*(   sum(  sum( (Bx.bb)./(Ix + backgx) + (By.bb)./(Iy + backgy) )  )   );
FIM(3,3) = 0.5.*(s^2).*(   sum(  sum( (Bx.cc)./(Ix + backgx) + (By.cc)./(Iy + backgy) )  )   );
FIM(4,4) = 0.5.*(s^2).*(   sum(  sum( (Bx.dd)./(Ix + backgx) + (By.dd)./(Iy + backgy) )  )   );
FIM(5,5) = 0.5.*(s^2).*(   sum(  sum( (Bx.ee)./(Ix + backgx) + (By.ee)./(Iy + backgy) )  )   );
FIM(6,6) = 0.5.*(s^2).*(   sum(  sum( (Bx.ff)./(Ix + backgx) + (By.ff)./(Iy + backgy) )  )   );

FIM(1,2) = (s^2).*(   sum(  sum( (Bx.ab)./(Ix + backgx) + (By.ab)./(Iy + backgy) )  )   );
FIM(1,3) = (s^2).*(   sum(  sum( (Bx.ac)./(Ix + backgx) + (By.ac)./(Iy + backgy) )  )   );
FIM(1,4) = (s^2).*(   sum(  sum( (Bx.ad)./(Ix + backgx) + (By.ad)./(Iy + backgy) )  )   );
FIM(1,5) = (s^2).*(   sum(  sum( (Bx.ae)./(Ix + backgx) + (By.ae)./(Iy + backgy) )  )   );
FIM(1,6) = (s^2).*(   sum(  sum( (Bx.af)./(Ix + backgx) + (By.af)./(Iy + backgy) )  )   );

FIM(2,3) = (s^2).*(   sum(  sum( (Bx.bc)./(Ix + backgx) + (By.bc)./(Iy + backgy) )  )   );
FIM(2,4) = (s^2).*(   sum(  sum( (Bx.bd)./(Ix + backgx) + (By.bd)./(Iy + backgy) )  )   );
FIM(2,5) = (s^2).*(   sum(  sum( (Bx.be)./(Ix + backgx) + (By.be)./(Iy + backgy) )  )   );
FIM(2,6) = (s^2).*(   sum(  sum( (Bx.bf)./(Ix + backgx) + (By.bf)./(Iy + backgy) )  )   );

FIM(3,4) = (s^2).*(   sum(  sum( (Bx.cd)./(Ix + backgx) + (By.cd)./(Iy + backgy) )  )   );
FIM(3,5) = (s^2).*(   sum(  sum( (Bx.ce)./(Ix + backgx) + (By.ce)./(Iy + backgy) )  )   );
FIM(3,6) = (s^2).*(   sum(  sum( (Bx.cf)./(Ix + backgx) + (By.cf)./(Iy + backgy) )  )   );

FIM(4,5) = (s^2).*(   sum(  sum( (Bx.de)./(Ix + backgx) + (By.de)./(Iy + backgy) )  )   );
FIM(4,6) = (s^2).*(   sum(  sum( (Bx.df)./(Ix + backgx) + (By.df)./(Iy + backgy) )  )   );

FIM(5,6) = (s^2).*(   sum(  sum( (Bx.ef)./(Ix + backgx) + (By.ef)./(Iy + backgy) )  )   );

FIM = FIM + FIM.';
end

