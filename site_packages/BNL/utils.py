# Big banner, with flowerbox
def bigBanner( x ): 
    sx = x.split( '\n' )
    l = max( map( len, sx ) )+4
    return '\n'.join( [ '+'+l*'-'+'+' ] + [ '|'+y.center( l )+'|' for y in sx ] + [ '+'+l*'-'+'+' ] )


# Small banner, with wings
def smallBanner( x, wingsize=10 ): 
    return wingsize*'*'+' '+x.replace( '\n', '; ' )+' '+wingsize*'*'