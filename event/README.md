The folder 'display/' contains routines to plot event hits and reconstructed tracks in 2d or 3d in the detector frame.

Example:

    $python3 display/2d.py -s izen -r 13 -n 10 -max 1000  #Plot 10 random events of the 'izen' run#13 inside the first 1000 events

or for specific eventid(s):

    $python3 display/2d.py -s izen -r 13 -e 98 542 -max 1000  #Plot events '98' & '542' inside the first 1000 events

Same syntax with 'display/3d.py'