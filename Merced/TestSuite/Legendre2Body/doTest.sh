#!
# bash for running the checks

for fil in Legendre2Body-Si28  Legendre2Body-Si28-Gauss6  fete2Body-Si28 \
  Lu175.linlin  Lu175.linlog  Legendre2Body-Si28.rel  Lu175.linlog.rel \
  Ag110m Ag110m.rel
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
