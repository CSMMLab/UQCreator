set terminal gif animate delay 10
set output "T.gif"
do for [i=0:176] {
   j=i*1000+1
   plot sprintf('RadiationHydrodynamics_1D_newExpectedValue_%d', j) using 1:2 with lines
   pause 0.1
   reread
 }
set output

#do for [i=1:1000:9001] {plot sprintf('RadiationHydrodynamics_1DExpectedValue_%d', i) using 1:2 with lines; pause 0.5} 
