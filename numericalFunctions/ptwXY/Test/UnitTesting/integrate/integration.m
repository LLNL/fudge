nn[x_] := ( y2 - y1 ) * ( x - x1 ) / ( x2 - x1 ) + y1
ng[x_] := y1 * Exp[ Log[ y2 / y1 ] * ( x - x1 ) / ( x2 - x1 ) ]
gn[x_] := ( y2 - y1 ) * Log[ x / x1 ] / Log[ x2 / x1 ] + y1
gg[x_] := y1 * (x / x1)^( Log[ y2 / y1 ] / Log[ x2 / x1 ] )
