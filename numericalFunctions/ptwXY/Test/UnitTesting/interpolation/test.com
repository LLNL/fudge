./$*
if [ $? -ne 0 ]; then echo "Failure for '$*'"; fi
