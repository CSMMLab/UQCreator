# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'pm3d.6.png'
set bar 1.000000 front
set border 4095 front lt black linewidth 1.000 dashtype solid
set style circle radius graph 0.02, first 0.00000, 0.00000 
set style ellipse size graph 0.05, 0.03, first 0.00000 angle 0 units xy
set style textbox transparent margins  1.0,  1.0 border
#unset logscale
#set view map scale 1
#set isosamples 100, 100
#unset surface 
#set style data pm3d
#set style function pm3d
#set xyplane relative 0
#unset paxis 1 tics
#unset paxis 2 tics
#unset paxis 3 tics
#unset paxis 4 tics
#unset paxis 5 tics
#unset paxis 6 tics
#unset paxis 7 tics
#set title "gray map" 
set xlabel "x" 
#set xrange [ -15.0000 : 15.0000 ] noreverse nowriteback
set ylabel "y" 
#set yrange [ -15.0000 : 15.0000 ] noreverse nowriteback
#set zrange [ -0.250000 : 1.00000 ] noreverse nowriteback
#set paxis 1 range [ * : * ] noreverse nowriteback
#set paxis 2 range [ * : * ] noreverse nowriteback
#set paxis 3 range [ * : * ] noreverse nowriteback
#set paxis 4 range [ * : * ] noreverse nowriteback
#set paxis 5 range [ * : * ] noreverse nowriteback
#set paxis 6 range [ * : * ] noreverse nowriteback
#set paxis 7 range [ * : * ] noreverse nowriteback
set pm3d implicit at b
#set palette positive nops_allcF maxcolors 0 gamma 1.5 gray
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front  noinvert bdefault
splot 'plotInXi' matrix with lines
