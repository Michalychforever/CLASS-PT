set term pdf enhanced dashed
#set termoption dashed
#set output 'PWS.eps'
set style line 1 default
set encoding utf8

set output '1.pdf'


#set xlabel "k (1/Mpc/h)"

set log x
set log y
set xrange [0.0001:1]
#plot "mvs_tCl.dat" u ($1):($2) w l lc rgb "red" lw 3 title "mvs tCl", "v1_tCl.dat" u ($1):($2) w l lc rgb "black" lw 1.2 title "v1 tCl"



plot "< join tCl.dat hector_pk_nl.dat" u ($1):(abs($2-$43)/abs($2)) w l lc rgb "red" lw 3 title "mvs-v1"
#, "< join mvs_tCl.dat v2_tCl.dat" u ($1):(abs($2-$43)/abs($2)) w l lc rgb "blue" lw 3 title "mvs-v2", "< join mvs_tCl.dat v3_tCl.dat" u ($1):(abs($2-$43)/abs($2)) w l lc rgb "green" lw 3 title "mvs-v3", "< join mvs_tCl.dat v3int_tCl.dat" u ($1):(abs($2-$43)/abs($2)) w l lc rgb "black" lw 3 title "mvs-v3int"
#plot "< join mvs_tCl.dat v1_tCl.dat" u ($1):(abs($3-$44)/abs($3)) w l lc rgb "red" lw 3 title "mvs-v1", "< join mvs_tCl.dat v2_tCl.dat" u ($1):(abs($3-$44)/abs($3)) w l lc rgb "blue" lw 3 title "mvs-v2", "< join mvs_tCl.dat v3_tCl.dat" u ($1):(abs($3-$44)/abs($3)) w l lc rgb "green" lw 3 title "mvs-v3"




#, "< join mvs_tCl.dat mvs_mPk.dat" u ($1):(abs($2-$43)/$2) w l lc rgb "green" lw 3 title "mPk-mvs"







