set term png
set output "graphic.png"
plot "latest" using 2  title "Best" with lines, '' using 3 title "Average" with lines
