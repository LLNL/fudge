<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="reactionSuite">
    <html>   
        <head>
            <!-- ********** basic css + local overrides *********** -->
            <link rel="stylesheet" type="text/css" href="css/base.css" />
            <style type="text/css">
    
                <!-- list with alternating background color -->
                li.reactionList { background: #9EBEE2; }
                li.reactionList:nth-child(odd) { background: #CCE0CB; }
    
                <!-- list with alternating background color -->
                li.particleList { background: #9EBEE2; }
                li.particleList:nth-child(odd) { background: #CCE0CB; }
    
                <!-- list with alternating background color -->
                li.summedReactionList { background: #9EBEE2; }
                li.summedReactionList:nth-child(odd) { background: #CCE0CB; }
    
                span.inChannel { color:red; }
                div.evaluationMeta { background: #CCE0CB; }
    
                <!-- table with alternating background color -->
                .r0 {background-color: #f9f9f9}
                .r1 {background-color: #f7f7e7}
            </style>

            <script type="text/javascript">
              <!-- ********** javascript for collapsing trees *********** -->
              /* CLOSED_IMAGE - the image to be displayed when the sublists are closed
               * OPEN_IMAGE   - the image to be displayed when the sublists are opened
               */
              CLOSED_IMAGE='images/plus.png';
              OPEN_IMAGE='images/minus.png';
  
              /* makeCollapsible - makes a list have collapsible sublists
               * 
               * listElement - the element representing the list to make collapsible
               */
              function makeCollapsible(listElement){
  
                  // removed list item bullets and the space they occupy
                  listElement.style.listStyle='none';
                  listElement.style.marginLeft='0';
                  listElement.style.paddingLeft='0';
  
                  // loop over all child elements of the list
                  var child=listElement.firstChild;
                  while (child!=null){
  
                      // only process li elements (and not text elements)
                      if (child.nodeType==1){
  
                          // build a list of child ol and ul elements and hide them
                          var list=new Array();
                          var grandchild=child.firstChild;
                          while (grandchild!=null){
                              if (grandchild.tagName=='OL' || grandchild.tagName=='UL'){
                                  grandchild.style.display='none';
                                  list.push(grandchild);
                              }
                              grandchild=grandchild.nextSibling;
                          }
  
                          // add toggle buttons
                          var node=document.createElement('img');
                          node.setAttribute('src',CLOSED_IMAGE);
                          node.setAttribute('class','collapsibleClosed');
                          node.onclick=createToggleFunction(node,list);
                          child.insertBefore(node,child.firstChild);
                      }
                      child=child.nextSibling;
                  }
              }
  
              /* createToggleFunction - returns a function that toggles the sublist display
               * 
               * toggleElement  - the element representing the toggle gadget
               * sublistElement - an array of elements representing the sublists that should
               *                  be opened or closed when the toggle gadget is clicked
               */
              function createToggleFunction(toggleElement,sublistElements){
  
                  return function(){
                      // toggle status of toggle gadget
                      if (toggleElement.getAttribute('class')=='collapsibleClosed'){
                          toggleElement.setAttribute('class','collapsibleOpen');
                          toggleElement.setAttribute('src',OPEN_IMAGE);
                      }else{
                          toggleElement.setAttribute('class','collapsibleClosed');
                          toggleElement.setAttribute('src',CLOSED_IMAGE);
                      }
  
                      // toggle display of sublists
                      for (var i=0;i &lt; sublistElements.length;i++){
                          sublistElements[i].style.display=(sublistElements[i].style.display=='block')?'none':'block';
                      }
  
                  }
  
              }

              <!-- for plotting 2-d datasets -->
              var pWidth = 701;
              var pHeight = 501;

              function popupPlot( data ) {
                win = window.open('','popup', 'height=480,width=720', false);
                var doc = win.document;
                doc.write( '&lt;html&gt;&lt;head&gt;' );
                doc.write( '&lt;title&gt; Plot &lt;/title&gt;' );
                doc.write( '&lt;/head&gt;&lt;body&gt;' );
                //doc.write( '&lt;p&gt; Hi there &lt;/p&gt;' );
                doc.write( '&lt;canvas id="plot" width="701" height="501"&gt;' );
                doc.write( '&lt;p&gt;No canvas support&lt;/p&gt;&lt;/canvas&gt;' );
                doc.write( '&lt;/body&gt;&lt;/html&gt;' );

                dataList = data.split(" ");
                var dataList2 = [];
                for (i=0; i&lt;dataList.length / 2; i++) {
                  // hack solution to convert strings to floats:
                  dataList2.push( [ Math.pow( dataList[2*i], 1.0 ),
                    Math.pow( dataList[2*i+1], 1.0 ) ] );
                }
                /* drawPlot( doc, data ); */
                draw( doc, dataList2, true );
                doc.close();
              }

              function drawPlot( doc, data ){
                var canvas = doc.getElementById("plot");
                var ctx = canvas.getContext("2d");
                ctx.lineWidth = 3;
                ctx.strokeStyle = "#FF0000";
                ctx.beginPath();
                ctx.moveTo(5,5);
                ctx.lineTo(395,195);
                ctx.stroke();
                ctx.save();
                ctx.restore();
              }

              function draw(document, points, markers) 
              {
                var canvas = document.getElementById('plot');
                var i, x, y;
                
                if(canvas.getContext) 
                {

                  var height = 500;
                  var width = 700;
                  var loffset = 45; // distance to axis from left side
                  var boffset = 45;

                  var xmin = 0; //points[0][0];
                  var xmax = 2.1e+7; //points[ points.length-1 ][0];
                  var ymin = points[0][1];
                  var ymax = ymin;
                  for (i=1; i&lt;points.length; i++) {
                    if (points[i][1] &lt; ymin) { ymin = points[i][1]; }
                    if (points[i][1] &gt; ymax) { ymax = points[i][1]; }
                    }


                  function calcX( x ) {
                    return ( x - xmin ) * (width - loffset - 15) / (xmax - xmin)
                  }
                  function calcY( y ) {
                    return (y - ymin ) * (height - boffset - 15) / (ymax - ymin)
                  }

                  var ctx = canvas.getContext('2d');
                  ctx.clearRect(0, 0, width+1, height+1);

                  ctx.strokeStyle = "#000000";
                  ctx.strokeRect(0, 0, width+1, height+1);
                  
                  // draw axis
                  ctx.lineWidth = 1;
                  ctx.strokeStyle = "#000000";
                  ctx.beginPath();
                  ctx.moveTo(loffset, 15);
                  ctx.lineTo(loffset, height-boffset+15);
                  ctx.moveTo(loffset-15, height-boffset);
                  ctx.lineTo(width-15, height-boffset);
                  ctx.stroke();
                  
                  // draw points			
                  ctx.save();
                  ctx.translate(loffset,height-boffset);  // move origin
                  ctx.scale(1.0, -1.0); // switch y-axis (start at bottom left instead)
                  ctx.strokeStyle = "#FF0000"; // red
                  ctx.lineWidth = 1;
                  ctx.beginPath();
                  
                  if(points.length > 0) 
                  {
                    x = calcX( points[0][0] );
                    y = calcY( points[0][1] );
                    ctx.moveTo(x, y);
                    for(i = 0; i &lt; points.length; i += 1) 
                    {
                      x = calcX( points[i][0] );
                      y = calcY( points[i][1] );
                      ctx.lineTo(x, y);
                    }
                    ctx.stroke();
                    
                    if (markers)
                    {
                      ctx.fillStyle="#009900";
                      for(i = 0; i &lt; points.length; i += 1) 
                      {
                        x = calcX( points[i][0] );
                        y = calcY( points[i][1] );
                        ctx.beginPath();
                        ctx.arc(x, y, 5, 0, 2*Math.PI, false);
                        ctx.fill();
                        ctx.stroke();
                      }
                    }
                    }

                  ctx.restore();
                  
                  // add text, needs to be last to overlap graph
                  ctx.font = 'italic 20px sans-serif';
                  ctx.textBaseline = 'top';
                  ctx.save();
                  ctx.rotate(-Math.PI/2);
                  ctx.fillText('cross section (b)', -300, 5);
                  ctx.restore();
                  ctx.fillText('Incident energy (eV)', 300, 475);
                  
                  // simple axes for now
                  ctx.font = '14px sans-serif';
                  ctx.fillText(xmin, loffset-3, 470);
                  ctx.fillText(xmax.toPrecision(3), width-loffset-15, 470);
                  ctx.fillText(ymin, loffset-35, height-boffset-3);
                  ctx.fillText(ymax.toPrecision(3), loffset-35, 15);

                }
              }


              <!-- for making a table of 2-d data -->
              function makeTable( data ){
                newwindow = window.open("","plot window", "height=480,width=640", false);
                newdoc = newwindow.document;
                var dataList = data.split(" ");
                var length = dataList.length / 2;
                newdoc.write( "&lt;html&gt;&lt;head/&gt;" );
                newdoc.write( "&lt;body&gt;&lt;h1&gt;Raw data:&lt;/h1&gt;&lt;table&gt;" );
                for (var i=0; i&lt;length; i++) {
                  newdoc.write( "&lt;tr&gt;&lt;td&gt;" + dataList[2*i] + "&lt;/td&gt;" );
                  newdoc.write( "&lt;td&gt;" + dataList[2*i+1] + "&lt;/td&gt;&lt;/tr&gt;" );
                  }
                newdoc.write( "&lt;/table&gt;&lt;/body&gt;&lt;/html&gt;" );
                }

            </script>
    
        </head>
          
        <!-- ********** main page formatting *********** -->
        <!-- the javascript is for collapsing a list. Register which lists should be collapsed: -->
        <body onLoad=" makeCollapsible(document.getElementById('reactions'));
              makeCollapsible(document.getElementById('summedReactions'));
              makeCollapsible(document.getElementById('Docs'));
              makeCollapsible(document.getElementById('particles'));
              makeCollapsible(document.getElementById('levels'));
              makeCollapsible(document.getElementById('gammas')); ">
            
              <h1>A Generalized Nuclear Data reactionSuite for 
                    <xsl:value-of select="@projectile"/> +
                    <xsl:value-of select="@target"/>
                  </h1>
                            
            <p>
                XML files (such as GNDS) can be easily transformed into html for viewing online.
                This page contains an overview of the GNDS file that was automatically generated
                using an XML stylesheet.
            </p>
    
            <div class="evaluationMeta">
                <h3>About this file:</h3>
                <ul>
                    <li>
                        Incident channel:
                        <span class="inChannel"><xsl:value-of select="@projectile"/> + <xsl:value-of select="@target"/></span>
                    </li>
                    <li>Format: <i><xsl:value-of select="@format"/></i></li>
                    <!-- <li>Schema: <i>???? (need a URL here!)</i></li> -->
                    <li>Temperature:  <i><xsl:value-of select="@temperature"/></i></li>
                    <li>
                        Available styles for this reactionSuite:  <xsl:value-of select="@styles"/>
                        <ul>
                            <xsl:apply-templates select="./styles/style"/>
                        </ul>
                    </li>
                </ul>
            </div>
    
            <xsl:apply-templates select="./documentations"/>
            <xsl:apply-templates select="./particles"/>
            <xsl:apply-templates select="./resonances"/>
    
            <h2>List of reactions:</h2>
            <table>
                <tr>
                    <td>&#160;</td>
                    <td>
                        <h3>Inclusive (summed) channels:</h3>
                        <ul id="summedReactions">
                            <xsl:apply-templates select="./summedReaction"/>
                        </ul>
                    </td>
                </tr>
            </table>
            <table>
                <tr>
                    <td>&#160;</td>
                    <td>
                        <h3>Exclusive (regular) channels:</h3>
                        <ul id="reactions">
                            <xsl:apply-templates select="./reaction"/>
                        </ul>
                    </td>
                </tr>
            </table>
        </body>
    </html>
</xsl:template>



<!-- now add template elements to html -->


<!-- ********** evalaution meta data formatting *********** -->
<xsl:template match="styles/style">
    <li>
        Style: <i><xsl:value-of select="@name"/></i>; 
        Library: <i><xsl:value-of select="@library"/></i>; 
        Evaluation version <i><xsl:value-of select="@version"/></i>
    </li>
</xsl:template>


<!-- ********** documentation formatting *********** -->
<xsl:template match="documentations">
    <h2>Documentation:</h2>
    <table>
        <tr>
            <td>&#160;</td>
            <td>
              <xsl:for-each select="documentation">
                <ul id="Docs">
                    <li>
                      <i><xsl:value-of select="@name"/>: click to expand</i>
                        <ul>
                            <pre>
                                <li>
                                    <xsl:value-of select="."/>
                                </li>
                            </pre>
                        </ul>
                    </li>
                  </ul>
                </xsl:for-each>
            </td>
        </tr>
    </table>
</xsl:template>


<!-- ********** particle properties formatting *********** -->
<xsl:template match="particles">
    <h2>Particles used in this evaluation:</h2>
    <table>
        <tr>
            <td>&#160;</td>
            <td>
                <ul id="particles">
                    <li>
                        <i>click to expand</i>
                        <ul>
                            <xsl:apply-templates select="./particle"/>
                        </ul>
                    </li>
                </ul>
            </td>
        </tr>
    </table>
</xsl:template>

<xsl:template match="particle">
    <li> &#160; 
        <b><xsl:value-of select="@name"/></b>:  
        (<xsl:value-of select="@genre"/>, mass = <xsl:value-of select="@mass"/>)
        <xsl:if test="@transportable='true'">&#160;-- <i>transportable</i>
        </xsl:if>
        <ul id="levels">
            <xsl:apply-templates select="./level"/>
        </ul>
    </li>
</xsl:template>

<xsl:template match="level">
    <li> &#160; 
        Level <b><xsl:value-of select="@name"/></b>
        <!-- # <xsl:value-of select="@index"/>: --> 
        at E[<xsl:value-of select="@label"/>] = <xsl:value-of select="@energy"/>
        <ul id="gammas">
            <xsl:apply-templates select="./gamma"/>
        </ul>
    </li>
</xsl:template>

<xsl:template match="gamma">
    <li class="gammaList"> &#160; 
        Eg = <xsl:value-of select="@energy"/> to <xsl:value-of select="@finalLevel"/> w/ 
        prob =  <xsl:value-of select="@probability"/>
    </li>
</xsl:template>


<!-- ********** resonance formatting *********** -->
<xsl:template match="resonances">
    <h2>Resonance region:</h2>
    <table>
        <tr>
            <td>&#160;</td>
            <td>
                <ul>
                    <xsl:apply-templates select="./scatteringRadius"/>
                    <xsl:apply-templates select="./resolved"/>
                    <xsl:apply-templates select="./unresolved"/>
                </ul>
            </td>
        </tr>
    </table>
</xsl:template>

<xsl:template match="scatteringRadius">
    <li>Scattering radius 
        (<xsl:value-of select="@lowerBound"/> - 
        <xsl:value-of select="@upperBound"/>):  
        <xsl:value-of select="@value"/>
    </li>
</xsl:template>

<xsl:template match="resolved">
    <xsl:choose>
        <xsl:when test="@multipleRegions='true'">
            <xsl:for-each select="region">
                <li>
                    Resolved region # <xsl:value-of select="@index"/>       
                    (<xsl:value-of select="@lowerBound"/> - 
                    <xsl:value-of select="@upperBound"/>): 
                    <xsl:value-of select="@nativeData"/>
                </li>
            </xsl:for-each>
        </xsl:when>
        <xsl:otherwise>
            <li>
                Resolved region         
                (<xsl:value-of select="@lowerBound"/> - 
                <xsl:value-of select="@upperBound"/>): 
                  <xsl:value-of select="@nativeData"/> format with
                    <xsl:value-of select=".//table/@rows"/> resonances.
            </li>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template match="unresolved">
    <li>
        Unresolved region
        (<xsl:value-of select="@lowerBound"/> - 
        <xsl:value-of select="@upperBound"/>): 
        <xsl:value-of select="@nativeData"/>
    </li>
</xsl:template>


<!-- ********** reaction channel formatting *********** -->
<xsl:template match="summedReaction">
    <li class="summedReactionList"> &#160; <xsl:value-of select="@name"/>
        <ul>
            <xsl:apply-templates select="./crossSection"/>
        </ul>
    </li>
</xsl:template>

<xsl:template match="reaction">
    <li class="reactionList"> &#160; <xsl:value-of select="@outputChannel"/>
        <ul>
            <xsl:apply-templates select="./crossSection"/>
            <xsl:apply-templates select="./outputChannel/product"/>
        </ul>
    </li>
</xsl:template>

  
<!-- ********** observable formatting *********** -->
<xsl:template match="crossSection">
    <li style="color:blue"> 
        Native cross section is in <xsl:value-of select="@nativeData"/> format:
        <ul>
            <xsl:apply-templates select="./linear"/>
            <xsl:apply-templates select="./pointwise"/>
            <xsl:apply-templates select="./piecewise"/>
            <xsl:apply-templates select="./resonancesWithBackground"/>
        </ul>
    </li>
</xsl:template>

<xsl:template match="product">
    <li style="color:green">
        Product <xsl:value-of select="@name"/>:
        <ul>
            <xsl:apply-templates select="./distributions"/>
            <xsl:apply-templates select="./multiplicity"/>
        </ul>
    </li>
</xsl:template>

<xsl:template match="distributions">
    <li><xsl:value-of select="@nativeData"/> distribution </li>
</xsl:template>

<xsl:template match="multiplicity">
    <li><xsl:value-of select="@nativeData"/> multiplicity </li>
</xsl:template>



<!-- ********** low-level data formatting *********** -->
<xsl:template match="linear">
    <xsl:variable name="xx" select="@xData"/>
    <li> 
      Pointwise lin-lin cross section:
      <button>
        <xsl:attribute name="onclick">
          makeTable('<xsl:value-of select="normalize-space(./data)"/>')
        </xsl:attribute>
        Show data
      </button>
      <button>
        <xsl:attribute name="onclick">
          <!-- drawGroovyPlot() -->
          popupPlot('<xsl:value-of select="normalize-space(./data)"/>')
        </xsl:attribute>
        Plot
      </button>
      <!--
          <ul>
            <li>
                <xsl:apply-templates select="./axes"/>
            </li>
            <li>
                Pointwise data:
                <table class="dataTable">
                    <tr>
                        <td>
                            <b>Axis #<xsl:value-of select="./axes/axis[@index='0']/@index"/></b>
                        </td>
                        <td>
                            <b>Axis #<xsl:value-of select="./axes/axis[@index='1']/@index"/></b>
                        </td>
                    </tr>
                    <tr>
                        <td>
                            <b><xsl:value-of select="./axes/axis[@index='0']/@label"/></b>
                        </td>
                        <td>
                            <b><xsl:value-of select="./axes/axis[@index='1']/@label"/></b>
                        </td>
                    </tr>
                    <tr>
                        <td>
                            <b>(<xsl:value-of select="./axes/axis[@index='0']/@unit"/>)</b>
                        </td>
                        <td>
                            <b>(<xsl:value-of select="./axes/axis[@index='1']/@unit"/>)</b>
                        </td>
                      </tr>
                      <xsl:call-template name="makeTable" select="./data">
                        <xsl:with-param name="data" select="substring-after(./data,' ')"></xsl:with-param>
                        <xsl:with-param name="delimiter" select="' '"></xsl:with-param>
                      </xsl:call-template>
                </table>
            </li>
          </ul>  -->
    </li>
</xsl:template>

<xsl:template match="pointwise">
    <li> 
        Pointwise cross section:
        <ul>
            <li>
                <xsl:apply-templates select="./axes"/>
            </li>
            <li>
                Pointwise data:
                <table class="dataTable">
                    <tr>
                        <td>
                            <b>Axis #<xsl:value-of select="./axes/axis[@index='0']/@index"/></b>
                        </td>
                        <td>
                            <b>Axis #<xsl:value-of select="./axes/axis[@index='1']/@index"/></b>
                        </td>
                    </tr>
                    <tr>
                        <td>
                            <b><xsl:value-of select="./axes/axis[@index='0']/@label"/></b>
                        </td>
                        <td>
                            <b><xsl:value-of select="./axes/axis[@index='1']/@label"/></b>
                        </td>
                    </tr>
                    <tr>
                        <td>
                            <b>(<xsl:value-of select="./axes/axis[@index='0']/@unit"/>)</b>
                        </td>
                        <td>
                            <b>(<xsl:value-of select="./axes/axis[@index='1']/@unit"/>)</b>
                        </td>
                    </tr>
                    <xsl:call-template name="makeTable" select="./data">
                        <xsl:with-param name="data" select="substring-after(./data,' ')"></xsl:with-param>
                        <xsl:with-param name="delimiter" select="' '"></xsl:with-param>
                    </xsl:call-template>
                    <!--
                    <xsl:call-template name="str.tokenize.to.columns">
                        <xsl:with-param name="string" select="./data"></xsl:with-param>
                        <xsl:with-param name="numcolumns" select="2"></xsl:with-param>
                    </xsl:call-template> 
                    -->
                </table>
            </li>
        </ul>
    </li>
</xsl:template>

<xsl:template match="piecewise">
    <li> Piecewise cross section:
        <ul>
            <li>
                <xsl:apply-templates select="./axes"/>
            </li>
            <li>
                Piecewise data:
                <ul>
                    <xsl:for-each select="./region">
                        <li>
                            Region # <xsl:value-of select="@index"/> 
                            (interpolate using <xsl:value-of select="interpolationAxes/@interpolation"/> to accuracy <xsl:value-of select="@accuracy"/>):
                            <table class="dataTable">
                                <tr>
                                    <td>
                                        <b>Axis #<xsl:value-of select="../axes/axis[@index='0']/@index"/></b>
                                    </td>
                                    <td>
                                        <b>Axis #<xsl:value-of select="../axes/axis[@index='1']/@index"/></b>
                                    </td>
                                </tr>
                                <tr>
                                    <td>
                                        <b><xsl:value-of select="../axes/axis[@index='0']/@label"/></b>
                                    </td>
                                    <td>
                                        <b><xsl:value-of select="../axes/axis[@index='1']/@label"/></b>
                                    </td>
                                </tr>
                                <tr>
                                    <td>
                                        <b>(<xsl:value-of select="../axes/axis[@index='0']/@unit"/>)</b>
                                    </td>
                                    <td>
                                        <b>(<xsl:value-of select="../axes/axis[@index='1']/@unit"/>)</b>
                                    </td>
                                </tr>
                                <xsl:call-template name="makeTable" select="./data">
                                    <xsl:with-param name="data" select="substring-after(./data,' ')"></xsl:with-param>
                                    <xsl:with-param name="delimiter" select="' '"></xsl:with-param>
                                </xsl:call-template>
                                <!--
                                <xsl:call-template name="str.tokenize.to.columns">
                                    <xsl:with-param name="string" select="./data"></xsl:with-param>
                                    <xsl:with-param name="numcolumns" select="2"></xsl:with-param>
                                </xsl:call-template>
                                -->
                            </table>
                        </li>
                    </xsl:for-each>
                </ul>
            </li>
        </ul>
    </li>
</xsl:template>

<xsl:template match="resonancesWithBackground">
    <li> This should be resonancesWithBackground data </li>
</xsl:template>

<xsl:template match="axes">
    Axes:
    <table border="1" cellpadding="2" cellspacing="0" class="axisTable">
        <tr>
            <td><B>Axis</B></td>
            <td><B>Quantity</B></td>
            <td><B>Unit</B></td>
            <td><B>Interpolation Scheme</B></td>
            <td><B>Reference Frame</B></td>
        </tr>
        <xsl:apply-templates select="./axis"/>
    </table>
</xsl:template>

<xsl:template match="axis">
    <tr>
        <td>axis # <xsl:value-of select="@index"/></td>
        <td><xsl:value-of select="@label"/></td>
        <td><xsl:value-of select="@unit"/></td>
        <td><xsl:value-of select="@interpolation"/></td>
        <td><xsl:value-of select="@frame"/></td>
    </tr>
</xsl:template>


<!-- Tokenizes 'string' then arranges in 'numcolumns' columns.  Does not do any checking. -->
<xsl:template name="str.tokenize.to.columns">
    <xsl:param name="string"></xsl:param>
    <xsl:param name="numcolumns"></xsl:param>
    <tr>
        <td>
            <xsl:value-of select="$numcolumns"></xsl:value-of>
        </td>
        <td>
            <xsl:value-of select="$string"></xsl:value-of>
        </td>
    </tr>
    <!-- <xsl:for-each select="tokenize('a b c',' ')"> This needs XPath 2.0 </xsl:for-each> -->
</xsl:template>


<!-- turn list of data into 2-column table: -->
<xsl:template name="makeTable">
    <xsl:param name="data"/>
    <xsl:param name="delimiter"/>
    <xsl:variable name="sub_string" select="substring-after(normalize-space($data),$delimiter)"/>
    <!--<xsl:variable name="td_count" select='number(13)'/>
    <xsl:if test="$td_count mod 2 = 1">Yippie</xsl:if>-->
    <tr>
        <td><xsl:value-of select="substring-before(normalize-space($data),$delimiter)"/></td>
        <xsl:choose>
            <xsl:when test="contains($sub_string,$delimiter)">
                <td><xsl:value-of select="substring-before($sub_string,$delimiter)"/></td>
            </xsl:when>
            <xsl:otherwise>
                <!-- last value in the list -->
                <td><xsl:value-of select="$sub_string"/></td>
            </xsl:otherwise>
        </xsl:choose>
    </tr>

    <!-- keep going while string still contains data -->
    <xsl:if test="substring-after($sub_string,$delimiter) != ''">
        <xsl:call-template name="makeTable">
            <xsl:with-param name="data" select="substring-after($sub_string,$delimiter)"/>
            <xsl:with-param name="delimiter" select="' '"/>
        </xsl:call-template>
    </xsl:if>
</xsl:template>

    
</xsl:stylesheet>
