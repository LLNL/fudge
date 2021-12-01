if [ "$2" == '-e' ]; then echo $1; fi
diff $1.out Answers/$1.out > /dev/null
if [ $? -ne 0 ]; then echo "Failure for '$*'"; fi
