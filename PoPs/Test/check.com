if [ "$2" == '-e' ]; then echo $1; fi
diff Answers/$1.out Outputs/$1.out> /dev/null
if [ $? -ne 0 ]; then echo "Failure for '$*'"; fi
