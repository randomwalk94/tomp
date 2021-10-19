function [dom] = make_domain(Oxd, Oyd, Lxd, Lyd, Nldx,Nldy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_domain:
% a simple routine that creates a structure 'dom' containing all the informations
% needed to describe an imaging 2D domain
%
% SYNOPSIS [dom] = make_domain(Oxd, Oyd, Lxd, Lyd, Nld)
%          
% INPUT - Oxd to Nld defines the domain to be searched (d stands  for domain !)
%              all those values are given with respect to 
%              the centralwavelength (found in the file name.sism)
%       - Oxd, Oyd = origin 
%       - Lxd, Lyd = size of the domain (also often named Nxd and Nyd in the Mfiles...)
%       - Nldx  = number of points per wavelength --> dx = lambda/Nldx
%       - Nldy  = number of points per wavelength --> dy = lambda/Nldy
%
% OUTPUT a struct containing the following fields:
%          - the inputs (of course)
%          - Ntxd, Ntyd, the size of the image
%          - x, y, the '1d grid' associated to the domain in each direction
%
% NOTE - typically an image will be defined using the coordinates x and y 
%           with a meshgrid [X, Y] = meshgrid(dom.x, dom.y)
%           and it will be defined as IM = zeros(size(X)).... etc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dom.Oxd = Oxd; 
dom.Oyd = Oyd; 
dom.Lxd = Lxd; 
dom.Lyd = Lyd; 
dom.Nldx = Nldx;
dom.Nldy = Nldy;

dom.Ntxd = Lxd*Nldx+1;  % number of steps in x (t stands for total)
dom.Ntyd = Lyd*Nldy+1;  % number of steps in y (t stands for total)

% the grid in both directions
dom.x = linspace(Oxd, Oxd + Lxd, dom.Ntxd);
dom.y = linspace(Oyd, Oyd + Lyd, dom.Ntyd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        that's all folks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
