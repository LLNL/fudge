./strings $1 -v > strings.dat
diff strings.out strings.dat > /dev/null 2>&1
exit $?
