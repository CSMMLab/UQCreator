# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'pm3d.33.png'
#set bar 1.000000 front
#set border 895 front lt black linewidth 1.000 dashtype solid
#set style circle radius graph 0.02, first 0.00000, 0.00000 
#set style ellipse size graph 0.05, 0.03, first 0.00000 angle 0 units xy
#set grid nopolar
#set grid xtics nomxtics ytics nomytics noztics nomztics \
# nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
#set grid front   lt 0 linewidth 0.500,  lt 0 linewidth 0.500
#set style line 100  linecolor rgb "#f0e442"  linewidth 0.500 dashtype solid pointtype 5 pointsize default pointinterval 0
#set style textbox transparent margins  1.0,  1.0 border
#unset logscale
#set view 20, 20, 1, 1
#set samples 11, 11
#set isosamples 11, 11
#unset surface 
#set style data pm3d
#set style function pm3d
#set xyplane relative 0
#set nomcbtics
#unset paxis 1 tics
#unset paxis 2 tics
#unset paxis 3 tics
#unset paxis 4 tics
#unset paxis 5 tics
#unset paxis 6 tics
#unset paxis 7 tics
#set title "Using interpolation with datafile; pm3d at s ftriangles interpolate 10,1" 
#set xlabel "X LABEL" 
#set ylabel "Y LABEL" 
#set paxis 1 range [ * : * ] noreverse nowriteback
#set paxis 2 range [ * : * ] noreverse nowriteback
#set paxis 3 range [ * : * ] noreverse nowriteback
#set paxis 4 range [ * : * ] noreverse nowriteback
#set paxis 5 range [ * : * ] noreverse nowriteback
#set paxis 6 range [ * : * ] noreverse nowriteback
#set paxis 7 range [ * : * ] noreverse nowriteback
#set lmargin  0
#set pm3d implicit at s
#set pm3d scansforward
#set pm3d interpolate 10,1 flush begin ftriangles noborder corners2color mean
#set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front  noinvert bdefault
#x = 0.0
## Last datafile plotted: "triangle.dat"
splot 'plotInXi' matrix with lines
