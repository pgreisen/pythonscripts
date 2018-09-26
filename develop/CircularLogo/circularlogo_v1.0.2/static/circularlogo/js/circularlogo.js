    
function cartesian2Polar(x, y){
    distance = Math.sqrt(x*x + y*y)
    radians = Math.atan2(y,x) //This takes y first
    polarCoor = { distance:distance, radians:radians }
    return polarCoor
};
  
/*
function connectLinks(cluster_nodes, motifGraph ) {
  var map = {},
      imports = [];

  // Compute a map from name to node.
  cluster_nodes.forEach(function(d) {
    map[d.name] = d;
  });

  var motifLinks = motifGraph.links;
  vislinks = [];
  for(i=0; i<motifLinks.length; i++){
    obj = motifLinks[i];
    newobj = {'source':map[obj.source], 'target':map[obj.target], 'value':obj.value};
    vislinks.push(newobj);
  };  
  return vislinks;
};
*/
  
function mouseovered(d, visParam) {
    var links = motifGraph.links;
    if(links) {
        links.sort(function(a, b) {
               //return a.value - b.value;
               return b.value - a.value;
        });
        var minval = d3.min(links, function(d){ return d.value; });
        var maxval = d3.max(links, function(d){ return d.value; });
        //var scale = d3.scale.linear()
        //                    .domain([minval, maxval])
        //                    .range([1, Math.min(20, arclength)]);  //map to "stroke-width"    
        var numColors = 9;
        //var heatmapColour = d3.scale.quantize()  
        var heatmapColour = d3.scale.linear()
          .domain([Math.log(minval), 0.5*(Math.log(minval) + Math.log(maxval)), Math.log(maxval)])
          .range(["red", "blue", "green"]);
          //.range(colorbrewer.Reds[numColors]);
          
        nid = d.data.index;
        var visnodes = d3.selectAll("#nodeG").selectAll("path");
                    visnodes.filter(function(n) { return n.data.index == nid;})
                            .style("stroke", "#d62728");
                    visnodes.filter(function(n) {return n.data.index != nid;})
                            .style("stroke", "steelblue");    
                            
        var vislinks = d3.select("#linkG").selectAll("path"); 
        vislinks
            .filter(function(l) { return l.target.id === nid || l.source.id === nid; })
            //.classed("link", function(l){return true;});
            //.style('stroke', function(l){return heatmapColour(Math.log(l.value))})
            .style('stroke', function(l){return "steelblue";})
            .style('stroke-opacity', 0.8);
            
        vislinks
            .filter(function(l) { return l.target.id != nid && l.source.id != nid; })
            //.classed("link--fade", function(l){return true;});
            .style('stroke-opacity', 0);
    }        
};    
    
function mouseouted(d) {
    var visnodes = d3.selectAll("#nodeG").selectAll("path");
    visnodes
         .style("stroke", "steelblue");
                        
    var vislinks = d3.select("#linkG").selectAll("path"); 
    if(vislinks){
        vislinks
            //.classed("link--highlight", false)
            //.classed("link--fade", false);
            .style("stroke", "steelblue")
            .style('stroke-opacity', 0.25);
    }        
};
        
    
function drawingMainSVG(motifGraph, visParam){ 
    var g_K = 4;  //the length of base table (A, T, C, G)          
    var g_width = visParam.width;
    var g_height = visParam.height;
    var g_padding = visParam.padding; //30
    var g_startAng = visParam.startAng;  // -3*Math.PI/4
    var g_endAng = visParam.endAng;   //3*Math.PI/4;
    var g_padAng = visParam.padAng;   //0.05
    
    
    var g_outerRadius = visParam.outerRadius;  //0.5*g_height - g_padding,
    var g_innerRadius = visParam.innerRadius;
    var g_offset = g_outerRadius - g_innerRadius;
    
    var g_oriX = 0.5*visParam.width;
    var g_oriY = 0.5*(visParam.height - visParam.padding);
    

    var arc = d3.svg.arc()
                    .innerRadius(g_innerRadius);
                    //.outerRadius(g_outerRadius);
    
    var pie = d3.layout.pie()
              .padAngle(g_padAng)
              .startAngle(g_startAng)
              .endAngle(g_endAng)
              .value(function(d, i) { return 1 - 0.000000000001*i; });
              //.value(function(d) { return 1; });
              //some bug exists in d3.pie when valuse are equal
              //here just for fix this bug by a temporary way now
    

    //Create SVG element
                
    var svg = d3.select("#svgcontainer")
                .append("svg")
                .attr("id", "root")
                .attr("width", g_width)
                .attr("height", g_height);      
                                   
     
    //Set up groups
    var arcs = svg.selectAll("g.arc")
                  .data(pie(motifGraph.nodes))
                  .enter()
                  .append("g")
                  .attr("id", "nodeG")
                  .attr("class", "node")
                  .attr("transform", "translate(" + g_oriX + "," + g_oriY + ")")
                  //.on("mouseover", mouseovered)
                  .on("mouseover", function(d){ mouseovered(d, visParam);})
                  .on("mouseout",  mouseouted);
                  
    
    svg.select("g").append('circle')
        .attr('cx', 0)
        .attr('cy', 0)
        .attr('r', 3)
        .style('fill', 'red');
            
   
    //Draw arc paths
    arcs.append("path")
        .each(function(d) { d.outerRadius = g_innerRadius + g_offset; })
        //.attr("stroke", "steelblue")  //#d62728
        //.attr("stroke-width", "3px")
        .attr("class", "npath")
        .attr("d", arc);
        
    
     
    //Locus
    var locus = arcs.append("g")
               .attr("transform", function(d, i) {
                                    arcO = arc.centroid(d);
                                    radians = cartesian2Polar(arcO[0], arcO[1]).radians;
                                    degree = radians * (180/Math.PI) + 90 ;
                                    myAction = "translate(" + arc.centroid(d) + ") rotate(" + degree + ")";
                                    return myAction;
                                }
                );

    function calcTransform(textnode, i){
        atom = textnode.parentNode.__data__;
        theta = Math.abs(atom.startAngle - atom.endAngle) - atom.padAngle;
        innerArcLen = Math.abs(theta) * g_innerRadius; 
        bbox = textnode.getBBox();
        ascent = bbox.height * (1 - 0.68 ) / 2;
        descent = ascent;
        charhgt = bbox.height * 0.68;
        xs = 0.75*innerArcLen/bbox.width;
        //ys = g_offset*atom.data.freq[i]/charhgt;
        ys = atom.data.bit/2*g_offset*atom.data.freq[i]/charhgt;
        //dy = g_offset*d3.sum(atom.data.freq.slice(0,i));
        dy = atom.data.bit/2*g_offset*d3.sum(atom.data.freq.slice(0,i));
        yt = 0.5*g_offset - dy;
        xt = 0;  
        //console.log(atom.data.index);
        //arcO = arc.centroid(atom);
        //console.log(arcO);
       return [xs, ys, xt, yt];
    }; 
    

    var bases = locus.selectAll("g")
                .append("text")
                .data(function(n) { 
                       return n.data.base; })         
                .enter().append("text")
                   .text(function(d) { return d; })      
                   .attr("text-anchor", "middle")                           
                   .attr("transform", function(d, i){
                         xyst = calcTransform(this, i);
                         return "translate(0," + xyst[3] + ") scale(" + xyst[0] + "," + xyst[1] + ")";
                     })
                   .attr("class", function(d, i) {
                      return "base" + d;  // The first smallest base  ##
                    });
           
    
    locus.each(function(d) {
        d3.select(this).append('circle')
            .attr('cx', 0)
            .attr('cy', 0)
            .attr('r', 1.5)
            .style('fill', 'red');
        
        
        var lbl = d3.select(this).append("text")
            .attr("x", 0 )
            .attr("y", - 0.5*g_offset - 5)
            .attr("text-anchor", "middle") 
            .text(d.data.label);   
        
        lbl.style("font-size", "0.8em");  
           
    });  

      
    //links section            
    //http://bl.ocks.org/mbostock/7607999  
    var nodeMap = {};
    locus.each(function(d, i) {
            ctm = this.getCTM();
            var point = document.getElementById('root').createSVGPoint();
            point.x = 0; 
            point.y = 0.5*g_offset + 5;
            var newPoint = point.matrixTransform(ctm);//new point after the transform
            x = newPoint.x;
            y = newPoint.y;
            //console.log(i, x, y);
            nodeMap[i] = {x:x,y:y, id:i};
        });   
    //console.log(nodeMap);
    
    var links = motifGraph.links;
    if(links) {
        links.sort(function(a, b) {
               //return a.value - b.value;
               return b.value - a.value;
            });
            
        
        vlinks = [];
        for(i=0; i<links.length; i++){
            obj = links[i];
            if(obj.value <= visParam.maxVal && obj.value >= visParam.minVal){
                newobj = {'source':nodeMap[obj.source], 'target':nodeMap[obj.target], 'value':obj.value};
                vlinks.push(newobj);
            };
        };    
      
        //Set up links group
        var arclinks = svg.append("g")
                            .attr("id", "linkG")
                            .attr("transform", "translate(" + 0 + "," + 0 + ")");
        
        //var minval = d3.min(links, function(d){ return d.value; });
        //var maxval = d3.max(links, function(d){ return d.value; });
        //var scale = d3.scale.linear()
        //                    .domain([minval, maxval])
        //                    .range([1, Math.min(20, arclength)]);  //map to "stroke-width"
        
        var scale = d3.scale.linear()
                            .domain([visParam.minVal, visParam.maxVal])
                            .range([visParam.minArcLen, visParam.maxArcLen]);  
                            
                            
        var bezierLine = d3.svg.line()
            .x(function(d) { return d[0]; })
            .y(function(d) { return d[1]; })
            .interpolate("basis"); 
            
         
        function customline(d){         
              var points = [
                            [ d.target.x, d.target.y ],
                            [ g_oriX, g_oriY ],
                            [ d.source.x, d.source.y ]
                          ];            
              return bezierLine(points);
        };

        var vislinks = arclinks.selectAll('path')
             .data(vlinks)
             //.data(links.filter(function(d){ return d.value > 6; }))
             .enter()
             .append('path')
             .attr("d", function(d) {return customline(d);})
             .attr("class", "link")
             .attr("stroke-width", function(d){return scale(d.value);});
    };

    /*
    var cluster_ang = (g_endAng - g_startAng)/Math.PI*180;
    //console.log(cluster_ang);
    var cluster = d3.layout.cluster()
                   .size([cluster_ang, g_innerRadius - 10]);
    
    var nodeMap = {};
    locus.each(function(d, i) {
            theta = 0.5*(d.startAngle + d.endAngle)/Math.PI*180 - g_startAng/Math.PI*180;
            //console.log(theta);
            nodeMap[i] = {x:theta};
            //ctm = this.getCTM();
            //var point = document.getElementById('root').createSVGPoint();
            //point.x = 0; 
            //point.y = 0.5*g_offset + 5;
            //var newPoint = point.matrixTransform(ctm);//new point after the transform
            //x = newPoint.x;
            //y = newPoint.y;
            //console.log(i, x, y);
            //nodeMap[i] = {x:x,y:y};   //need to transform to polar coordinates against the original point (g_oriX, g_oriY)
        });   
    //console.log(nodeMap);
    
    var bundle = d3.layout.bundle();

    var line = d3.svg.line.radial()
        .interpolate("bundle")
        .tension(.85)
        .radius(function(d) { return d.y; })
        .angle(function(d) { return d.x / 180 * Math.PI; }); 

    var cluster_nodes = cluster.nodes(d3tree);

    //we need to adjust each node position
    cluster_nodes.forEach(function(d){
        if(d.name in nodeMap){
            d.x = nodeMap[d.name].x;
        };
    });
    
    //var cluster_links = cluster.links(cluster_nodes);  
    var cluster_links = connectLinks(cluster_nodes, motifGraph );

    var cluster_rot = g_startAng/Math.PI*180;
    var vnode = svg.append("g")
                   .attr("transform", "translate(" + g_oriX + "," + g_oriY + "), rotate(" + cluster_rot + ")")
                   .selectAll(".node");
                   
    var vlink = svg.append("g")
                   .attr("transform","translate(" + g_oriX + "," + g_oriY + "), rotate(" + cluster_rot + ")")
                   .selectAll(".link");
    
    var minval = d3.min(motifGraph.links, function(d){ return d.value; });
    var maxval = d3.max(motifGraph.links, function(d){ return d.value; });
    var scale = d3.scale.linear()
                        .domain([minval, maxval])
                        .range([1, 10]);  //map to "stroke-width"

    
    var vlink = vlink
      .data(bundle(cluster_links))
    .enter().append("path")
      .each(function(d) { d.source = d[0], d.target = d[d.length - 1]; })
      .attr("class", "link")
      .attr("stroke-width", function(d) {
                               source = d[0].name; 
                               target = d[d.length - 1].name;
                               return 2;
                               //return scale(d.value);
       })
      .attr("d", line);   
      
      
    var vnode = vnode
      .data(cluster_nodes.filter(function(n) { return !n.children; }))
    .enter().append("text")
      .attr("class", "node")
      .attr("dy", ".31em")
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (d.y + 8) + ",0)" + (d.x < 180 ? "" : "rotate(180)"); })
      .style("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      //.text(function(d) { return d.key; })
      .text(function(d) { return d.name; })
      .on("mouseover", mouseovered)
      .on("mouseout", mouseouted);  
    */         

}



function formatSelect(){
    var format = $("input[type='radio'][name='format']:checked").val();
    if(format == "json") {
        document.getElementById("id_pseudoA").disabled = true;
        document.getElementById("id_pseudoT").disabled = true;
        document.getElementById("id_pseudoC").disabled = true;
        document.getElementById("id_pseudoG").disabled = true;
        document.getElementById("id_method").disabled = true;
        document.getElementById("id_gbkgA").disabled = true;
        document.getElementById("id_gbkgT").disabled = true;
        document.getElementById("id_gbkgC").disabled = true;
        document.getElementById("id_gbkgG").disabled = true;
        document.getElementById("id_pvalue").disabled = true;
    }else{
        document.getElementById("id_pseudoA").disabled = false;
        document.getElementById("id_pseudoT").disabled = false;
        document.getElementById("id_pseudoC").disabled = false;
        document.getElementById("id_pseudoG").disabled = false;
        document.getElementById("id_method").disabled = false;
        document.getElementById("id_gbkgA").disabled = false;
        document.getElementById("id_gbkgT").disabled = false;
        document.getElementById("id_gbkgC").disabled = false;
        document.getElementById("id_gbkgG").disabled = false;
        document.getElementById("id_pvalue").disabled = false;
    };
};

function loadJsonDemo(){
    var json_instances = ''
        +'{                                                                                                               \n'
        +'    "id": "Toy motif",                                                                                          \n'
        +'    "background":{"key":["A","T","C","G"],"val":[0.25,0.25,0.25,0.25]},                                         \n'
        +'    "pseudocounts":{"key":["A","T","C","G"],"val":[0.25,0.25,0.25,0.25]},                                       \n'        
        +'    "nodes": [                                                                                                  \n'
        +'        {"index": 0, "label": "1", "bit": 1.9, "base": ["A", "T", "C", "G"], "freq": [0.0, 0.0, 0.2, 0.8]  },   \n'
        +'        {"index": 1, "label": "2", "bit": 1.8, "base": ["T", "C", "A", "G"], "freq": [0.1, 0.2, 0.3, 0.4]  },   \n'
        +'        {"index": 2, "label": "3", "bit": 1.6, "base": ["C", "G", "A", "T"], "freq": [0.05, 0.05, 0.2, 0.7]},   \n'
        +'        {"index": 3, "label": "4", "bit": 0.8, "base": ["G", "C", "A", "T"], "freq": [0.1, 0.2, 0.3, 0.4]  },   \n'
        +'        {"index": 4, "label": "5", "bit": 1.5, "base": ["G", "T", "C", "A"], "freq": [0.04, 0.06, 0.2, 0.7]  }, \n'
        +'        {"index": 5, "label": "6", "bit": 0.8, "base": ["T", "C", "G", "A"], "freq": [0.1, 0.2, 0.3, 0.4]  },   \n'
        +'        {"index": 6, "label": "7", "bit": 0.9, "base": ["G", "T", "C", "A"], "freq": [0.1, 0.2, 0.2, 0.5]  },   \n'
        +'        {"index": 7, "label": "8", "bit": 1.8, "base": ["C", "G", "A", "T"], "freq": [0.02, 0.08, 0.2,0.7]  },  \n'
        +'        {"index": 8, "label": "9", "bit": 0.8, "base": ["T", "G", "A", "C"], "freq": [0.1, 0.2, 0.3, 0.4]  },   \n'
        +'        {"index": 9, "label": "10", "bit": 0.8, "base": ["T", "A", "G", "C"], "freq": [0.1, 0.2, 0.3, 0.4]  }    \n'
        +'      ],                                                                                                        \n'
        +'    "links": [                                                                                                  \n'
        +'        {"source": 0, "target": 1, "value": 2},                                                                 \n'
        +'        {"source": 2, "target": 4, "value": 8},                                                                 \n'
        +'        {"source": 2, "target": 7, "value": 10},                                                                \n'
        +'        {"source": 7, "target": 8, "value": 5},                                                                 \n'
        +'        {"source": 3, "target": 8, "value": 5}                                                                  \n'
        +'      ]                                                                                                         \n'
        +'}                                                                                                               \n';
    document.getElementById("id_inputText").value = json_instances;    
    document.getElementById("id_format_1").checked = true;    
    var hidElem = document.getElementById("id_hidden_config");
    hidElem.value = "minArcLen:5;maxArcLen:30";  //;minVal:3.0;maxVal:8.0
    formatSelect();
};                 


function loadFastaDemo(){
    var demo_instances = prepareFasta();
    document.getElementById("id_inputText").value = demo_instances;   
    document.getElementById("id_format_0").checked = true;      
    var hidElem = document.getElementById("id_hidden_config");
    hidElem.value = "case:'HNF6Demo'";
    //hidElem.value = "minVal:102.0";    //chi-square 
    //hidElem.value = "minVal:0.035";  //mutual information
    document.getElementById("id_title").value = "HNF6";
    formatSelect();   
};

