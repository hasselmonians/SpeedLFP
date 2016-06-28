function savefig_fname(filename)
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will one-by-one bring each currently open
% figure to the foreground.  Then it will supply a text prompt in the main
% console window asking you for a filename.  It will save that figure to
% that filename in the .fig, .emf, .png, and .eps formats.  
%
% The formats that it saves in can be changed by commenting out or adding
% lines below.
%
% Copyright 2010 Matthew Guidry 
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)

set(gcf,'renderer','painters');
set(gcf,'PaperPositionMode','auto')
saveas(gcf, [filename '.fig']) %Matlab .FIG file
%saveas(hfigs(m), [filename '.emf']) %Windows Enhanced Meta-File (best for powerpoints)
saveas(gcf, [filename '.png']) %Standard PNG graphics file (best for web)
saveas(gcf, [filename '.eps'],'epsc2') %Standard PNG graphics file (best for web)
