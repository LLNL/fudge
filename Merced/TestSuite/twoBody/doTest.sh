#!
# bash for running the checks

for fil in elastic-C12  elastic-C12-Gauss6  inelastic-Kr78  inelastic-Kr78-Gauss6 \
  elastic-C12.rel  inelastic-Kr78.rel
do
  echo "doing $fil"
  ../get_transfer $fil > $fil.info
  diff utfil out.$fil > out.$fil.Diff
  if [ -s out.$fil.Diff ]
  then
    echo out.$fil different
    mv utfil out.$fil.new
  fi
  echo ""
done
rm *Diff
