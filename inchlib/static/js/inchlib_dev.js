/**
* InCHlib is an interactive JavaScript library which facilitates data
* visualization and exploration by means of a cluster heatmap. InCHlib
* is a versatile tool, and its use is not limited only to chemical or
* biological data. Source code, tutorial, documentation, and example
* data are freely available from InCHlib website <a
* href="http://openscreen.cz/software/inchlib"
* target=blank>http://openscreen.cz/software/inchlib</a>. At the
* website, you can also find a Python script <a
* href="http://openscreen.cz/software/inchlib/inchlib_clust"
* target=blank>inchlib_clust</a> which performs data clustering and
* prepares <a href="http://openscreen.cz/software/inchlib/input_format"
* target=blank>input data for InCHlib</a>.
* 
* @author <a href="mailto:ctibor.skuta@img.cas.cz">Ctibor Škuta</a>
* @author <a href="mailto:petr.bartunek@img.cas.cz">Petr Bartůněk</a>
* @author <a href="mailto:svozild@vscht.cz">Daniel Svozil</a>
* @version 1.0.0
* @category 1
* @license InCHlib - Interactive Cluster Heatmap Library http://openscreen.cz/software/inchlib Copyright 2014, Ctibor Škuta, Petr Bartůněk, Daniel Svozil Licensed under the MIT license.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
* OTHER DEALINGS IN THE SOFTWARE.
* 
* @requires <a href='http://code.jquery.com/jquery-2.0.3.min.js'>jQuery Core 2.0.3</a>
* @dependency <script language="JavaScript" type="text/javascript" src="http://code.jquery.com/jquery-2.0.3.min.js"></script>
* 
* @requires <a href='http://kineticjs.com/'>KineticJS 5.0.0</a>
* @dependency <script language="JavaScript" type="text/javascript" src="http://openscreen.cz/software/inchlib/static/js/kinetic-v5.0.0.min.js"></script>
*
* @param {Object} options An object with the options for the InCHlib component.
*
* @option {string} target
*   identifier of the DIV tag where the component should be displayed

* @option {boolean} [column_dendrogram=false]
*   turn on/off the column dendrogram

* @option {boolean} [count_column=false]
*   turn on/off the count column

* @option {boolean} [dendrogram=true]
*   turn on/off the row dendrogram

* @option {string} [font="Trebuchet&nbsp;MS"]
*   font family

* @option {string} [heatmap_colors="Greens"]
*   the heatmap color scale

* @option {number} [heatmap_part_width=0.7]
*   define the heatmap part width from the width of the whole graph

* @option {string} [highlight_colors="Reds"]
*   color scale for highlighted rows

* @option {obejct} [highlighted_rows=[]]
*   array of row IDs to highlight

* @option {boolean} [independent_columns=true]
*   determines whether the color scale is based on the values from all columns together or for each column separately

* @option {string} [label_color=grey]
*   color of column label

* @option {number} [max_column_width=100]
*   maximum column width in pixels

* @option {number} [max_height=800]
*   maximum graph height in pixels

* @option {number} [max_row_height=25]
*   maximum row height in pixels

* @option {boolean} [metadata=false]
*   turn on/off the metadata

* @option {string} [metadata_colors="Oranges"]
*   the metadata color scale

* @option {number} [min_row_height=false]
*   minimum row height in pixels

* @option {number} [width="the width of target DIV"]
*   width of the graph in pixels

* @option {boolean} [heatmap=true]
*   turn on/off the heatmap

* @option {string} [heatmap_font_color="black"]
*   the color of the text values in the heatmap

* @option {string} [count_column_colors="Reds"]
*   the color scale of count column

* @option {string} [values_center="median"]
*   the center point of color scale when using scale with three colors (possible values: median, mean, zero)

* @option {boolean} [draw_row_ids=false]
*   draws the row IDs next to the heatmap when there is enough space to visualize them

* @option {boolean} [header_as_heatmap_row=false]
*   use the heatmap header to draw one special heatmap row at the top (useful as column metadata active/inactive)

* @option {string} [header_row_colors="YlOrB"]
*   color scale when using header as a heatmap row
* 
* @example
*       window.instance = new InCHlib({
*                target : "YourOwnDivId",
*                metadata: true, 
*                max_height: 800,
*                width: 700,
*                metadata_colors: "RdLrBu"
*            });
*       instance.read_data_from_file("../biojs/data/chembl_gr.json");
*       instance.draw();
*/

var InCHlib;

(function($){
  InCHlib = function(settings){
      // When measuring the rendering duration
      // this.start = new Date().getTime();
      var self = this;
      this.user_settings = settings;
      this.target_element = $("#" + settings.target);
      var target_width = this.target_element.width();
      this.target_element.css({"position": "relative"});

      /**
      * Default values for the settings
      * @name InCHlib#settings
      */
      
      this.settings = {
          "target" : "YourOwnDivId",
          "heatmap" : true,
          "heatmap_header": true,
          "dendrogram": true,
          "metadata": false,
          "column_metadata": false,
          "column_metadata_row_height": 8,
          "column_metadata_colors": "RdLrBu",
          "max_height" : 800,
          "width" : target_width,
          "heatmap_colors" : "Greens",
          "heatmap_font_color" : "black",
          "heatmap_part_width" : 0.7,
          "column_dendrogram" : false,
          "independent_columns" : true,
          "metadata_colors" : "Reds",
          "highlight_colors" : "Oranges",
          "highlighted_rows" : [],
          "label_color": "#9E9E9E",
          "count_column": false,
          "count_column_colors": "Reds",
          "min_row_height": false,
          "max_row_height": 25,
          "max_column_width": 150,
          "font": "Helvetica",
          "values_center": "median",
          "draw_row_ids": false,
          "show_export_button": true,
          "max_quantile": 100,
          "min_quantile": 0,
          "middle_quantile": 50,
      };

      $.extend(this.settings, settings);
      this.settings.width = (settings.max_width && settings.max_width < target_width)?settings.max_width:this.settings.width;
      this.settings.heatmap_part_width = (this.settings.heatmap_part_width>0.9)?0.9:this.settings.heatmap_part_width;

      this.header_height = 150;
      this.footer_height = 70;
      this.dendrogram_heatmap_distance = 5;
      this.heatmap_width = 0;
      this.highlighted_row = null;
      this.leaves_y_coordinates = {};

      /**
      * Default function definitions for the InCHlib events
      * @name InCHlib#events
      */
      this.events = {
          /**
            * @name InCHlib#row_onclick
            * @event
            * @param {function} function() callback function for click on the heatmap row event
            * @eventData {array} array array of object IDs represented by row
            * @eventData {object} event event object

            * @example 
            * instance.events.row_onclick = (
            *    function(object_ids, evt) {
            *       alert(object_ids);
            *    }
            * ); 
            * 
            */
          "row_onclick": function(object_ids, evt){
              return;
          },

          /**
            * @name InCHlib#row_onmouseover
            * @event
            * @param {function} function() callback function for mouse cursor over the heatmap row event
            * @eventData {array} array array of object IDs represented by row
            * @eventData {object} event event object

            * @example 
            * instance.events.row_onmouseover = (
            *    function(object_ids, evt) {
            *       alert(object_ids);
            *    }
            * ); 
            * 
            */
          "row_onmouseover": function(object_ids, evt){
              return;
          },

          /**
            * @name InCHlib#row_onmouseout
            * @event
            * @param {function} function() callback function for mouse cursor out of the heatmap row event
            * @eventData {object} event event object

            * @example 
            * instance.events.row_onmouseout = (
            *    function(evt) {
            *       alert("now");
            *    }
            * ); 
            * 
            */
          "row_onmouseout": function(evt){
              return;
          },

          /**
            * @name InCHlib#dendrogram_node_onclick
            * @event
            * @param {function} function() callback function for dendrogram node click event
            * @eventData {array} array array of object IDs represented by the node
            * @eventData {string} node_id Id of the dendrogram node
            * @eventData {object} event event object

            * @example 
            * instance.events.dendrogram_node_onclick = (
            *    function(object_ids, node_id, evt) {
            *    alert(node_id + ": " + object_ids.length+" rows");
            *    }
            * ); 
            * 
            */
          "dendrogram_node_onclick": function(object_ids, node_id, evt){
              return;
          },

          /**
            * @name InCHlib#column_dendrogram_node_onclick
            * @event
            * @param {function} function() callback function for column dendrogram click event
            * @eventData {array} array array of column indexes or headers represented by the node
            * @eventData {string} node_id Id of the dendrogram node
            * @eventData {object} event event object

            * @example 
            * instance.events.column_dendrogram_node_onclick = (
            *    function(column_ids, node_id, evt) {
            *    alert(node_id + ": " + column_ids.length+" columns");
            *    }
            * ); 
            * 
            */
          "column_dendrogram_node_onclick": function(object_ids, node_id, evt){
              return;
          },

          /**
            * @name InCHlib#dendrogram_node_highlight
            * @event
            * @param {function} function() callback function for the dendrogram node highlight event
            * @eventData {array} array array of object IDs represented by row
            * @eventData {string} node_id Id of the dendrogram node
            * @eventData {object} event event object

            * @example 
            * instance.events.dendrogram_node_highlight = (
            *    function(object_ids, node_id, evt) {
            *       alert(node_id + ": " + object_ids.length+" rows");
            *    }
            * ); 
            * 
            */
          "dendrogram_node_highlight": function(object_ids, node_id){
              return;
          },

          /**
            * @name InCHlib#dendrogram_node_unhighlight
            * @event
            * @param {function} function() callback function for the dendrogram node unhighlight event
            * @eventData {string} node_id Id of the dendrogram node

            * @example 
            * instance.events.dendrogram_node_unhighlight = (
            *    function(node_id) {
            *       alert(node_id);
            *    }
            * ); 
            * 
            */
          "dendrogram_node_unhighlight": function(node_id){
              return;
          },

          /**
            * @name InCHlib#heatmap_onmouseout
            * @event
            * @param {function} function() callback function for mouse cursor out of hte heatmap area
            * @eventData {object} event event object

            * @example 
            * instance.events.heatmap_onmouseout = (
            *    function(evt) {
            *       alert("now");
            *    }
            * ); 
            * 
            */
          "heatmap_onmouseout": function(evt){
              return;
          },

          /**
            * @name InCHlib#on_zoom
            * @event
            * @param {function} function() callback function for zoom icon click event
            * @eventData {string} node_id Id of the dendrogram node

            * @example 
            * instance.events.on_zoom = (
            *    function(node_id) {
            *       alert(node_id);
            *    }
            * ); 
            * 
            */
          "on_zoom": function(node_id){
              return;
          },

          /**
            * @name InCHlib#on_unzoom
            * @event
            * @param {function} function() callback function for unzoom icon click event
            * @eventData {string} node_id Id of the dendrogram node

            * @example 
            * instance.events.on_unzoom = (
            *    function(node_id) {
            *       alert(node_id);
            *    }
            * ); 
            * 
            */
          "on_unzoom": function(node_id){
              return;
          },

          /**
            * @name InCHlib#on_refresh
            * @event
            * @param {function} function() callback function for refresh icon click event
            * @eventData {object} event event object
            * @example 
            * instance.events.on_refresh = (
            *    function() {
            *       alert("now");
            *    }
            * ); 
            * 
            */
          "on_refresh": function(){
              return;
          },

          /**
            * @name InCHlib#empty_space_onclick
            * @event
            * @param {function} function() callback function for click on empty(inactive) space in the visualization (e.g., around the heatmap)
            * @eventData {object} event event object

            * @example 
            * instance.events.empty_space_onclick = (
            *    function(evt) {
            *       alert("now");
            *    }
            * ); 
            * 
            */
          "empty_space_onclick": function(evt){
              return;
          }

      }

      /**
      * Default color scales
      * @name InCHlib#colors
      */
      this.colors = {
              "YlGn": {"start": {"r":255, "g": 255, "b": 204}, "end": {"r": 35, "g": 132, "b": 67}},
              "GnBu": {"start": {"r":240, "g": 249, "b": 232}, "end": {"r": 43, "g": 140, "b": 190}},
              "BuGn": {"start": {"r":237, "g": 248, "b": 251}, "end": {"r": 35, "g": 139, "b": 69}},
              "PuBu": {"start": {"r":241, "g": 238, "b": 246}, "end": {"r": 5, "g": 112, "b": 176}},
              "BuPu": {"start": {"r":237, "g": 248, "b": 251}, "end": {"r": 136, "g": 65, "b": 157}},
              "RdPu": {"start": {"r":254, "g": 235, "b": 226}, "end": {"r": 174, "g": 1, "b": 126}},
              "PuRd": {"start": {"r":241, "g": 238, "b": 246}, "end": {"r": 206, "g": 18, "b": 86}},
              "OrRd": {"start": {"r":254, "g": 240, "b": 217}, "end": {"r": 215, "g": 48, "b": 31}},
              "Purples2": {"start": {"r":242, "g": 240, "b": 247}, "end": {"r": 106, "g": 81, "b": 163}},
              "Blues": {"start": {"r":239, "g": 243, "b": 255}, "end": {"r": 33, "g": 113, "b": 181}},
              "Greens": {"start": {"r":237, "g": 248, "b": 233}, "end": {"r": 35, "g": 139, "b": 69}},
              "Oranges": {"start": {"r":254, "g": 237, "b": 222}, "end": {"r": 217, "g": 71, "b": 1}},
              "Reds": {"start": {"r":254, "g": 229, "b": 217}, "end": {"r": 203, "g": 24, "b": 29}},
              "Greys": {"start": {"r":247, "g": 247, "b": 247}, "end": {"r": 82, "g": 82, "b": 82}},
              "PuOr": {"start": {"r":230, "g": 97, "b": 1}, "end": {"r": 94, "g": 60, "b": 153}},
              "BrBG": {"start": {"r":166, "g": 97, "b": 26}, "end": {"r": 1, "g": 133, "b": 113}},
              "RdBu": {"start": {"r":202, "g": 0, "b": 32}, "end": {"r": 5, "g": 113, "b": 176}},
              "RdGy": {"start": {"r":202, "g": 0, "b": 32}, "end": {"r": 64, "g": 64, "b": 64}},
              "BuYl": {"start": {"r": 5, "g": 113, "b": 176}, "end": {"r": 250, "g": 233, "b": 42}},
              "YlOrR": {"start": {"r":255, "g": 255, "b": 178}, "end": {"r": 227, "g": 26, "b": 28}, "middle": {"r": 204, "g": 76, "b": 2}},
              "YlOrB": {"start": {"r":255, "g": 255, "b": 212}, "end": {"r": 5, "g": 112, "b": 176}, "middle": {"r": 204, "g": 76, "b": 2}},
              "PRGn2": {"start": {"r":123, "g": 50, "b": 148}, "end": {"r": 0, "g": 136, "b": 55}, "middle": {"r":202, "g": 0, "b": 32}},
              "PiYG2": {"start": {"r":208, "g": 28, "b": 139}, "end": {"r": 77, "g": 172, "b": 38}, "middle": {"r":255, "g": 255, "b": 178},},
              "YlGnBu": {"start": {"r":255, "g": 255, "b": 204}, "end": {"r": 34, "g": 94, "b": 168}, "middle": {"r": 35, "g": 132, "b": 67}},
              "RdYlBu": {"start": {"r":215, "g": 25, "b": 28}, "end": {"r": 44, "g": 123, "b": 182}, "middle": {"r":255, "g": 255, "b": 178}},
              "RdYlGn": {"start": {"r":215, "g": 25, "b": 28}, "end": {"r": 26, "g": 150, "b": 65}, "middle": {"r":255, "g": 255, "b": 178}},
              "BuWhRd": {"start": {"r": 33, "g": 113, "b": 181}, "middle": {"r": 255, "g": 255, "b": 255}, "end": {"r":215, "g": 25, "b": 28}},
              "RdLrBu": {"start": {"r":215, "g": 25, "b": 28}, "middle": {"r":254, "g": 229, "b": 217}, "end": {"r": 44, "g": 123, "b": 182}},
              "RdBkGr": {"start": {"r":215, "g": 25, "b": 28}, "middle": {"r": 0, "g": 0, "b": 0}, "end": {"r": 35, "g": 139, "b": 69}},
              "RdLrGr": {"start": {"r":215, "g": 25, "b": 28}, "middle": {"r":254, "g": 229, "b": 217}, "end": {"r": 35, "g": 139, "b": 69}},
      };

      /**
      * Default kineticjs objects references
      * @name InCHlib#objects_ref
      */
      this.objects_ref = {
          "tooltip_label": new Kinetic.Label({
                              opacity: 1,
                              listening: false,
                           }),

          "tooltip_tag": new Kinetic.Tag({
                              fill: this.settings.label_color,
                              pointerWidth: 10,
                              pointerHeight: 10,
                              lineJoin: 'round',
                              listening: false,
                          }),
      
          "tooltip_text": new Kinetic.Text({
                              fontFamily: this.settings.font,
                              fontSize: 12,
                              padding: 8,
                              fill: 'white',
                              fontStyle: "bold",
                              listening: false,
                              align: "center",
                              lineHeight: 1.2,
                          }),

          "node": new Kinetic.Line({
                              stroke: "grey",
                              strokeWidth: 2,
                              lineCap: 'sqare',
                              lineJoin: 'round',
                              listening: false
                          }),

          "node_rect" : new Kinetic.Rect({
                              fill: "white",
                              opacity: 0,
                          }),

          "icon_overlay": new Kinetic.Rect({
                              width: 32,
                              height: 32,
                              opacity: 0,
                          }),

          "heatmap_value": new Kinetic.Text({
                              fontFamily: this.settings.font,
                              fill: this.settings.heatmap_font_color,
                              fontStyle: "bold",
                              listening: false,
                          }),

          "heatmap_line": new Kinetic.Line({
                             lineCap: 'butt',
                             value: false,
                          }),

          "column_header": new Kinetic.Text({
                              fontFamily: this.settings.font,
                              fontStyle: "bold",
                              fill: 'black',
                           }),

          "count": new Kinetic.Text({
                          fontSize: 20,
                          fontFamily: this.settings.font,
                          fontStyle: 'bold',
                          fill: 'grey',
                          listening: false,
                       }),

          "cluster_overlay": new Kinetic.Rect({
                                  fill: "white",
                                  opacity: 0.5,
                              }),

          "cluster_border": new Kinetic.Line({
                                  stroke: "black",
                                  strokeWidth: 1,
                              }),

          "icon": new Kinetic.Path({
                      fill: "grey",
                      scale: 1,
                  }),

          "rect_gradient": new Kinetic.Rect({
                              x: 0,
                              y: 80,
                              width: 100,
                              height: 20,
                              fillLinearGradientStartPoint: {x: 0, y: 80},
                              fillLinearGradientEndPoint: {x: 100, y: 80},
                              stroke: "#D2D2D2",
                              strokeWidth: "1px"
                          }),
      };
      
  }

  InCHlib.prototype._update_user_settings = function(settings){
    var updated_settings = {}, key;
    for(var i = 0, keys=Object.keys(settings), len = keys.length; i < len; i++){
      key = keys[i];
      if(this.user_settings[key] !== undefined && this.user_settings[key] !== settings[key] && this.user_settings[key] === true){
        updated_settings[key] = false;
      }
      else if(this.user_settings[key] === undefined){
        updated_settings[key] = settings[key];
      }
    }
    $.extend(this.settings, updated_settings);
  }

  /**
    * Read data from JSON variable.
    * 
    * @param {object} [variable] Clustering in proper JSON format.
    */
  InCHlib.prototype.read_data = function(json){
      this.json = json;
      this.data = this.json.data;
      var settings = {};
      if(this.json["metadata"] !== undefined){
        this.metadata = this.json.metadata;
        settings.metadata = true;
      }
      else{
        settings.metadata = false;
      }
      if(this.json["column_dendrogram"] !== undefined){
        this.column_dendrogram = this.json.column_dendrogram;
        settings.column_dendrogram = true;
      }
      else{
        settings.column_dendrogram = false;
      }
      if(this.json["column_metadata"] !== undefined){
        this.column_metadata = this.json.column_metadata;
        settings.column_metadata = true;
      }
      else{
        settings.column_metadata = false;
      }

      this._update_user_settings(settings);
      // this.data.nodes = this.data.data.nodes;
  }

  /**
    * Read data from JSON file.
    * 
    * @param {string} [filename] Path to the JSON data file.
    *
    */
  InCHlib.prototype.read_data_from_file = function(json){
      var self = this;
      
      $.ajax({
          type: 'GET',
          url: json,
          dataType: 'json',
          success: function(json){
            self.read_data(json);
          },
          async: false
      });
  }

  InCHlib.prototype._add_prefix = function(){
      var i, id;
      var data = {};
      for(var i = 0, keys = Object.keys(this.data.nodes), len = keys.length; i < len; i++){
          id = [this.settings.target, keys[i]].join("#");
          data[id] = this.data.nodes[keys[i]];
          
          if(data[id]["parent"] !== undefined){
              data[id].parent = [this.settings.target, data[id].parent].join("#");
          }

          if(data[id]["count"] != 1){
              data[id].left_child = [this.settings.target, data[id].left_child].join("#");
              data[id].right_child = [this.settings.target, data[id].right_child].join("#");
          }
      }

      this.data.nodes = data;

      if(this.settings.metadata){
        var metadata = {};
        for(var i = 0, keys = Object.keys(this.metadata.nodes), len = keys.length; i < len; i++){
            id = [this.settings.target, keys[i]].join("#");
            metadata[id] = this.metadata.nodes[keys[i]];
        }
        this.metadata.nodes = metadata;
      }

      return;
  }

  InCHlib.prototype._get_root_id = function(nodes){
      var root_id;
      for(var i = 0, keys = Object.keys(nodes), len = keys.length; i < len; i++){
          if(nodes[keys[i]]["parent"] === undefined){
              root_id = keys[i];
              break;
          }
      }
      return root_id;
  }

  InCHlib.prototype._get_dimensions = function(){
      var dimensions = {"data": 0, "metadata": 0, "overall": 0}, key, keys, i;

      for(i = 0, keys = Object.keys(this.data.nodes), len = keys.length; i < len; i++){
          key = keys[i];
          if(this.data.nodes[key].count == 1){
              dimensions["data"] = this.data.nodes[key].features.length;
              break;
          }
      }

      if(this.settings.metadata){
        key = Object.keys(this.metadata.nodes)[0];
        dimensions["metadata"] = this.metadata.nodes[key].length;
      }
      
      dimensions["overall"] = dimensions["data"] + dimensions["metadata"];
      return dimensions;
  }

  InCHlib.prototype._get_min_max_middle = function(data){
      var i, len;
      var min_max_middle = [];
      var all = [];

      for(i = 0, len = data.length; i<len; i++){
          all = all.concat(data[i].filter(function(x){return x !== null}));
      }

      var len = all.length;
      all.sort(function(a,b){return a - b});
      min_max_middle.push((this.settings.min_quantile > 0)?all[this._hack_round(len*this.settings.min_quantile/100)]:Math.min.apply(null, all));
      min_max_middle.push((this.settings.max_quantile < 100)?all[this._hack_round(len*this.settings.max_quantile/100)]:Math.max.apply(null, all));
      min_max_middle.push((this.settings.middle_quantile != 50)?all[this._hack_round(len*this.settings.middle_quantile/100)]:this._middle2fnc(all));
      return min_max_middle;
  }

  InCHlib.prototype._get_data_min_max_middle = function(data, axis){
      if(axis === undefined){
          axis = "column";
      }

      var self = this;
      var columns = data;
      var i, j, value, len;
      var data_length = data[0].length;

      if(axis == "column"){
          columns = [];

          for(i = 0; i<data_length; i++){
              columns.push([]);
          }

          for(i = 0; i<data.length; i++){
              for(j = 0; j < data_length; j++){
                  value = data[i][j];
                  if(value !== null && value !== undefined){
                      columns[j].push(value);
                  }
              }
          }
      }

      var data2descs = {}
      var data_min_max_middle = [], min, max, middle;

      for(i = 0; i<columns.length; i++){
          if(this._is_number(columns[i][0])){
              columns[i] = columns[i].map(parseFloat);
              columns[i].sort(function(a,b){return a - b});
              len = columns[i].length;
              max = (this.settings.max_quantile < 100)?columns[i][this._hack_round(len*this.settings.max_quantile/100)]:Math.max.apply(null, columns[i]);
              min = (this.settings.min_quantile > 0)?columns[i][this._hack_round(len*this.settings.min_quantile/100)]:Math.min.apply(null, columns[i]);
              middle = (this.settings.middle_quantile != 50)?columns[i][this._hack_round(len*this.settings.middle_quantile/100)]:this._middle2fnc(columns[i]);
              data2descs[i] = {"min": min, "max": max, "middle": middle};
          }
          else{
              var hash_object = this._get_hash_object(columns[i]);
              min = 0;
              max = this._hack_size(hash_object)-1;
              middle = max/2;
              data2descs[i] = {"min": min, "max": max, "middle": middle, "str2num": hash_object};
          }
      }

      return data2descs;
  }

  InCHlib.prototype._get_hash_object = function(array){
      var i, count=0, hash_object = {};

      for(i = 0; i<array.length; i++){
          if(hash_object[array[i]] === undefined){
              hash_object[array[i]] = count;
              count++;
          }
      }
      return hash_object;
  }

  InCHlib.prototype._get_max_length = function(items){
      var lengths = items.map(function(x){return (""+x).length});
      var max = Math.max.apply(Math, lengths);
      return max;
  }

  InCHlib.prototype._get_max_value_length = function(){
      var nodes = this.data.nodes;
      var max_length = 0;
      var node_data, key;

      for(var i = 0, keys = Object.keys(nodes), len = keys.length; i < len; i++){
          key = keys[i];
          if(nodes[key].count == 1){
              node_data = nodes[key].features;
              for(var j = 0, len_2 = node_data.length; j < len_2; j++){
                  if((""+node_data[j]).length > max_length){
                      max_length = (""+node_data[j]).length;
                  }
              }
          }
      }
      
      if(this.settings.metadata){
          nodes = this.metadata.nodes;
          for(var i = 0, keys = Object.keys(nodes), len = keys.length; i < len; i++){
              key = keys[i];
              node_data = nodes[key];
              for(var j = 0, len_2 = node_data.length; j < len_2; j++){
                  if((""+node_data[j]).length > max_length){
                      max_length = (""+node_data[j]).length;
                  }
              }
          }
      }
      return max_length;
  }

  InCHlib.prototype._preprocess_heatmap_data = function(){
      var heatmap_array = [], i, j = 0, keys, key, len, data, node;

      for(i = 0, keys = Object.keys(this.data.nodes), len = keys.length; i < len; i++){
          key = keys[i];
          node = this.data.nodes[key];
          if(node.count == 1){
              data = node.features;
              heatmap_array.push([key]);
              heatmap_array[j].push.apply(heatmap_array[j], data);
              if(this.settings.metadata){
                  heatmap_array[j].push.apply(heatmap_array[j], this.metadata.nodes[key]);
              }
              j++;
          }
      }
      return heatmap_array;
  }

  InCHlib.prototype._reorder_heatmap = function(column_index){
      this.leaves_y_coordinates = {};
      
      if(this.ordered_by_index == column_index){
          this.heatmap_array.reverse();
      }
      else{
          if(this._is_number(this.heatmap_array[0][column_index])){
            this.heatmap_array.sort(function(a,b){return (a[column_index] == null)?-1:(b[column_index] == null)?1:a[column_index] - b[column_index]});
          }
          else{
            this.heatmap_array.sort(function(a,b){return (a[column_index] == null)?-1:(b[column_index] == null)?1:(a[column_index] > b[column_index])?1:(a[column_index] < b[column_index])?-1:0});
          }
      }

      var pixels_for_leaf = this.pixels_for_leaf;
      var y = pixels_for_leaf/2 + this.header_height;

      for(var i = 0, len = this.heatmap_array.length; i<len; i++){
          this.leaves_y_coordinates[this.heatmap_array[i][0]] = y;
          y += pixels_for_leaf;
      }

      this.ordered_by_index = column_index;
      return;
  }

  /**
    * Draw already read data (from file/JSON variable).
    */
  InCHlib.prototype.draw = function(){
      this._add_prefix();
      this.zoomed_clusters = {"row": [], "column": []};
      this.heatmap_array = this._preprocess_heatmap_data();

      if(this.settings.heatmap){
        this.last_column = null;
        this.on_features = {"data":[], "metadata":[]}
        this.dimensions = this._get_dimensions();
        this.column_metadata_rows = (this.settings.column_metadata)?this.column_metadata.length:0;
        this.column_metadata_height = this.column_metadata_rows * this.settings.column_metadata_row_height;
        this._set_heatmap_settings();
        this._adjust_leaf_size(this.heatmap_array.length);
      }
      else{
          this.dimensions = {"data": 0, "metadata": 0, "overall": 0};
      }

      if(this.settings.column_dendrogram && this.heatmap_header){
          this.footer_height = 150;
      }

      this.stage = new Kinetic.Stage({
          container: this.settings.target,
      });

      this.settings.height = this.heatmap_array.length*this.pixels_for_leaf+this.header_height+this.footer_height;
      this.stage.setWidth(this.settings.width);
      this.stage.setHeight(this.settings.height);
      this._draw_stage_layer();

      if(this.settings.dendrogram){
          this.root_id = this._get_root_id(this.data.nodes);
          this.distance = this.settings.width - this.heatmap_width-this.right_margin;

          var time_1 = new Date().getTime();
          this._draw_row_dendrogram(this.root_id);
          var time_2 = new Date().getTime();
          console.log("Dendrogram:", time_2 - time_1);

          var time_3 = new Date().getTime();
          this._draw_heatmap();
          var time_4 = new Date().getTime();
          console.log("Heatmap:", time_4 - time_3);
      }
      else{
        this.settings.heatmap_part_width = 1;
        this.distance = this._hack_round((this.settings.width-this.heatmap_width)/2);
        this._reorder_heatmap(0);
        this.ordered_by_index = 0;
        this._draw_heatmap();
      }

      if(this.settings.column_dendrogram){
        this.column_root_id = this._get_root_id(this.column_dendrogram.nodes);
        this.nodes2columns = false;
        this.columns_start_index = 0;
        this._draw_column_dendrogram(this.column_root_id);
      }
      this._draw_navigation();
      this.highlight_rows(this.settings.highlighted_rows);
      // When measuring the rendering duration
      // this.end = new Date().getTime();
      // console.log(this.end - this.start);
  }

  InCHlib.prototype._draw_row_dendrogram = function(node_id){
      this.dendrogram_layer = new Kinetic.Layer();
      var node = this.data.nodes[node_id];
      var count = node.count;

      this.distance_step = this.distance/node.distance;
      this.last_highlighted_cluster = null;
      this.leaves_y_coordinates = {};
      this.objects2leaves = {};

      this._adjust_leaf_size(count);
      this.settings.height = count*this.pixels_for_leaf+this.header_height+this.footer_height+this.column_metadata_height;
      
      this.stage.setWidth(this.settings.width);
      this.stage.setHeight(this.settings.height);

      var current_left_count = 0;
      var current_right_count = 0;
      var y = this.header_height + this.column_metadata_height + this.pixels_for_leaf/2;
      
      if(node.count > 1){
          current_left_count = this.data.nodes[node.left_child].count;
          current_right_count = this.data.nodes[node.right_child].count;
      }
      this._draw_row_dendrogram_node(node_id, node, current_left_count, current_right_count, 0, y);
      this.middle_item_count = (this.min_item_count+this.max_item_count)/2;
      this._draw_distance_scale(node.distance);

      this.cluster_layer = new Kinetic.Layer();
      this.dendrogram_overlay_layer = new Kinetic.Layer();
      this.dendrogram_hover_layer = new Kinetic.Layer();
      this.stage.add(this.dendrogram_layer, this.cluster_layer, this.dendrogram_overlay_layer, this.dendrogram_hover_layer);

      var self = this;

      this._bind_dendrogram_hover_events(this.dendrogram_layer);
      this._bind_dendrogram_hover_events(this.dendrogram_overlay_layer);
      
      this.dendrogram_layer.on("click", function(evt){
          self._dendrogram_layers_click(this, evt);
      });

      this.dendrogram_overlay_layer.on("click", function(evt){
          self._dendrogram_layers_click(this, evt);
      });
      
      this.timer = 0;
      this.dendrogram_layer.on("mousedown", function(evt){
        self._dendrogram_layers_mousedown(this, evt);
      });

      this.dendrogram_layer.on("mouseup", function(evt){
        self._dendrogram_layers_mouseup(this, evt);
      });

      this.dendrogram_overlay_layer.on("mousedown", function(evt){
        self._dendrogram_layers_mousedown(this, evt);
      });

      this.dendrogram_overlay_layer.on("mouseup", function(evt){
        self._dendrogram_layers_mouseup(this, evt);
      });
      
  }

  InCHlib.prototype._draw_row_dendrogram_node = function(node_id, node, current_left_count, current_right_count, x, y){
      if(node.count != 1){
          var node_neighbourhood = this._get_node_neighbourhood(node, this.data.nodes);
          var right_child = this.data.nodes[node.right_child];
          var left_child = this.data.nodes[node.left_child];

          var y1 = this._get_y1(node_neighbourhood, current_left_count, current_right_count);
          var y2 = this._get_y2(node_neighbourhood, current_left_count, current_right_count);
          var x1 = this._hack_round(this.distance - this.distance_step*node.distance);
          x1 = (x1 == 0)? 2: x1;
          
          var x2 = x1;
          var left_distance = this.distance - this.distance_step*this.data.nodes[node.left_child].distance;
          var right_distance = this.distance - this.distance_step*this.data.nodes[node.right_child].distance;

          if(right_child.count == 1){
              y2 = y2 + this.pixels_for_leaf/2;
          }

          this.dendrogram_layer.add(this._draw_horizontal_path(node_id, x1, y1, x2, y2, left_distance, right_distance));
          this._draw_row_dendrogram_node(node.left_child, left_child, current_left_count - node_neighbourhood.left_node.right_count, current_right_count + node_neighbourhood.left_node.right_count, left_distance, y1);
          this._draw_row_dendrogram_node(node.right_child, right_child, current_left_count + node_neighbourhood.right_node.left_count, current_right_count - node_neighbourhood.right_node.left_count, right_distance, y2);
      }
      else{
          var objects = node.objects;
          this.leaves_y_coordinates[node_id] = y;

          for(var i = 0, len = objects.length; i<len; i++){
              this.objects2leaves[objects[i]] = node_id;
          }

          var count = node.objects.length;
          if(count<this.min_item_count){
              this.min_item_count = count;
          }
          if(count>this.max_item_count){
              this.max_item_count = count;
          }
      }
  }

  InCHlib.prototype._draw_stage_layer = function(){
      this.stage_layer = new Kinetic.Layer();
      var stage_rect = new Kinetic.Rect({
                                  x: 0,
                                  y: 0,
                                  width: this.settings.width,
                                  height: this.settings.height,
                                  opacity: 0,
                              });
      this.stage_layer.add(stage_rect);
      stage_rect.moveToBottom();
      this.stage.add(this.stage_layer);

      var self = this;

      this.stage_layer.on("click", function(evt){
          self.unhighlight_cluster();
          self.unhighlight_column_cluster();
          self.events.empty_space_onclick(evt);
      });
  }

  InCHlib.prototype._draw_column_dendrogram = function(node_id){
      this.column_dendrogram_layer = new Kinetic.Layer();
      this.column_x_coordinates = {};
      var node = this.column_dendrogram.nodes[node_id];
      this.current_column_count = node.count;
      this.current_column_ids = [];
      this.vertical_distance = this.header_height-5;
      this.vertical_distance_step = this.vertical_distance/node.distance;

      this.last_highlighted_column_cluster = null;
      var current_left_count = this.column_dendrogram.nodes[node.left_child].count;
      var current_right_count = this.column_dendrogram.nodes[node.right_child].count;
      this._draw_column_dendrogram_node(node_id, node, current_left_count, current_right_count, 0, 0);
      this.stage.add(this.column_dendrogram_layer);

      if(!this.nodes2columns){
        this.nodes2columns = this._get_nodes2columns();
      }
      
      this.column_dendrogram_overlay_layer = new Kinetic.Layer();
      this.stage.add(this.column_dendrogram_overlay_layer);

      this._bind_dendrogram_hover_events(this.column_dendrogram_layer);
      this._bind_dendrogram_hover_events(this.column_dendrogram_overlay_layer);
      
      var self = this;

      this.column_dendrogram_layer.on("click", function(evt){
          self._column_dendrogram_layers_click(this, evt);
      });

      this.column_dendrogram_overlay_layer.on("click", function(evt){
          self._column_dendrogram_layers_click(this, evt);
      });

      this.column_dendrogram_layer.on("mousedown", function(evt){
        self._column_dendrogram_layers_mousedown(this, evt);
      });

      this.column_dendrogram_layer.on("mouseup", function(evt){
        self._dendrogram_layers_mouseup(this, evt);
      });

      this.column_dendrogram_overlay_layer.on("mousedown", function(evt){
        self._column_dendrogram_layers_mousedown(this, evt);
      });

      this.column_dendrogram_overlay_layer.on("mouseup", function(evt){
        self._dendrogram_layers_mouseup(this, evt);
      });
  }

  InCHlib.prototype._get_nodes2columns = function(){
    var coordinates = [];
    var coordinates2nodes = {};
    var nodes2columns = {};
    var key, value, i;

    for(i = 0, keys = Object.keys(this.column_x_coordinates), len = keys.length; i < len; i++){
      key = keys[i];
      value = this.column_x_coordinates[key];
      coordinates2nodes[value] = key;
      coordinates.push(value);
    }
    coordinates.sort(function(a,b){return a - b});

    for(i = 0, len = coordinates.length; i<len; i++){
      nodes2columns[coordinates2nodes[coordinates[i]]] = i;
    }
    return nodes2columns;
  }

  InCHlib.prototype._bind_dendrogram_hover_events = function(layer){
      var self = this;

      layer.on("mouseover", function(evt){
          self._dendrogram_layers_mouseover(this, evt);
      });

      layer.on("mouseout", function(evt){
          self._dendrogram_layers_mouseout(this, evt);
      });    
  }

  InCHlib.prototype._delete_layers = function(to_destroy, to_remove_children){
      for(var i = 0, len = to_destroy.length; i < len; i++){
        if(to_destroy[i] !== undefined){
          to_destroy[i].destroy();
        }
      }

      if(to_remove_children !== undefined){
        for(var i = 0, len = to_remove_children.length; i < len; i++){
          to_remove_children[i].removeChildren();
          to_remove_children[i].draw();
        }
      }
  }

  InCHlib.prototype._delete_all_layers = function(){
      this.stage.destroyChildren();
  }

  InCHlib.prototype._adjust_leaf_size = function(leaves){
      this.pixels_for_leaf = (this.settings.max_height-this.header_height-this.footer_height-this.column_metadata_height)/leaves;
      if(this.pixels_for_leaf > 3*this.pixels_for_dimension){
        this.pixels_for_leaf = this.pixels_for_dimension;
      }

      if(this.pixels_for_leaf > this.settings.max_row_height){
          this.pixels_for_leaf = this.settings.max_row_height;
      }

      if(this.settings.min_row_height > this.pixels_for_leaf){
          this.pixels_for_leaf = this.settings.min_row_height;
      }
      return;
  }

  InCHlib.prototype._adjust_dimension_size = function(dimensions){
      this.heatmap_width = (this.settings.width - this.right_margin)*this.settings.heatmap_part_width;
      this.distance = this.settings.width - this.heatmap_width - this.right_margin;
      this.heatmap_distance = this.distance + this.dendrogram_heatmap_distance;
      this.pixels_for_dimension = dimensions?this.heatmap_width/dimensions:0;

      if(this.settings.max_column_width && this.settings.max_column_width < this.pixels_for_dimension){
        this.pixels_for_dimension = this.settings.max_column_width;
        this.heatmap_width = dimensions*this.pixels_for_dimension;
        this.distance = this.settings.width - this.heatmap_width - this.right_margin;
        this.heatmap_distance = this.distance + this.dendrogram_heatmap_distance;
      }
      return;
  }

  InCHlib.prototype._set_heatmap_settings = function(){
      var i, keys, key, len, node;

      this.header = [];
      for(i = 0; i<this.dimensions["overall"]; i++){
          this.header.push("");
      }

      this.heatmap_header = false;
      this.metadata_header = false;
      this.current_label = null;

      if(this.data.feature_names !== undefined){
          this.heatmap_header = this.data.feature_names;
          for(i=0; i<this.dimensions["data"]; i++){
              this.header[i] = this.heatmap_header[i];
          }
      }

      var data = [];
      for(i = 0, keys = Object.keys(this.data.nodes), len = keys.length; i < len; i++){
          node = this.data.nodes[keys[i]];
          if(node.count == 1){
              data.push(node.features);
          };
      }
      
      this.data_descs = {};
      if(this.settings.independent_columns){
          this.data_descs = this._get_data_min_max_middle(data);
      }
      else{
          var min_max_middle = this._get_min_max_middle(data);
          for(i = 0; i < this.dimensions["data"]; i++){
            this.data_descs[i] = {"min": min_max_middle[0], "max": min_max_middle[1], "middle": min_max_middle[2]};
          }
      }

      if(this.settings.metadata){

          if(this.metadata.feature_names){
              this.metadata_header = this.metadata.feature_names;

              for(i=0; i<this.dimensions["metadata"]; i++){
                  this.header[this.dimensions["data"]+i] = this.metadata_header[i];
              }
          }

          var metadata = [];

          for(i = 0, keys = Object.keys(this.metadata.nodes), len = keys.length; i < len; i++){
              metadata.push(this.metadata.nodes[keys[i]]);
          }
          this.metadata_descs = this._get_data_min_max_middle(metadata);
      }

      if(this.settings.count_column){
          this.max_item_count = 1;
          this.min_item_count = 1;
          this.dimensions["overall"]++;
          this.header.push("Count");
      }

      this.settings.heatmap_part_width = this.dimensions["overall"]?this.settings.heatmap_part_width:0;
      this.right_margin = 100;
      this.heatmap_width = this.dimensions["overall"]?this._hack_round(this.settings.heatmap_part_width*(this.settings.width-this.right_margin)):0;

      this.features = {};
      for(i=0; i<this.dimensions["overall"]; i++){
          this.features[i] = 1;
      }
      this._adjust_dimension_size(this.dimensions["overall"]);
      this._set_on_features();
      this.top_heatmap_distance = this.header_height + this.column_metadata_height + this.settings.column_metadata_row_height/2;
  }

  InCHlib.prototype._set_on_features = function(features){
    var key;
    if(features === undefined){
      var features = [];
      for(var i = 0, keys = Object.keys(this.features), len = keys.length; i < len; i++){
        key = keys[i];
        if(this.features[key] === 1){
          features.push(key);
        }
      }
    }

    this.on_features = {"data":[], "metadata":[]}

    for(var i = 0, len = features.length; i < len; i++){
      key = features[i];
      if(key < this.dimensions["data"]){
          this.on_features["data"].push(parseInt(key));
      }
      else{
          if(parseInt(key) <= this.dimensions["overall"]-1){
              this.on_features["metadata"].push(parseInt(key)-this.dimensions["data"]);
          }
      }
    }
  }

  InCHlib.prototype._draw_heatmap = function(){
      if(!this.settings.heatmap || this.dimensions["overall"]==0){
          return;
      }

      var heatmap_row, row_id, col_number, col_label, row_values, y;
      this.heatmap_layer = new Kinetic.Layer();
      this.heatmap_overlay = new Kinetic.Layer();
      this.current_draw_values = true;
      this.max_value_length = this._get_max_value_length();
      this.value_font_size = this._get_font_size(this.max_value_length, this.pixels_for_dimension, this.pixels_for_leaf, 12);

      if(this.value_font_size < 4){
          this.current_draw_values = false;
      }
      this.heatmap_distance = this._hack_round(this.distance + this.dendrogram_heatmap_distance);

      var x1 = this.heatmap_distance;
      var current_leaves_y = [];
      var self = this;

      var time_1 = new Date().getTime();
      for(var i = 0, keys = Object.keys(this.leaves_y_coordinates), len = keys.length; i < len; i++){
          key = keys[i];
          y = this.leaves_y_coordinates[key];
          heatmap_row = this._draw_heatmap_row(key, x1, y);
          this.heatmap_layer.add(heatmap_row);
          current_leaves_y.push([key, y]);
          this._bind_row_events(heatmap_row);
      }

      var time_2 = new Date().getTime();
      console.log("Heatmap rows:", time_2 - time_1);

      if(this.settings.column_metadata){
          var time_3 = new Date().getTime();
          this.column_metadata_descs = this._get_data_min_max_middle(this.column_metadata, "row");
          y1 = this.header_height + 0.5*this.settings.column_metadata_row_height;

          for(var i = 0, len = this.column_metadata.length; i < len; i++){
              heatmap_row = this._draw_column_metadata_row(this.column_metadata[i], i, x1, y1);
              this.heatmap_layer.add(heatmap_row);
              this._bind_row_events(heatmap_row);
              y1 = y1 + this.settings.column_metadata_row_height;
          }
          var time_4 = new Date().getTime();
          console.log("Column metadata:", time_4 - time_3);
          
      }

      if(this.settings.draw_row_ids){
          this._draw_row_ids(current_leaves_y);
      }

      this.highlighted_rows_layer = new Kinetic.Layer();
      this.stage.add(this.heatmap_layer, this.heatmap_overlay, this.highlighted_rows_layer);

      this.highlighted_rows_layer.moveToTop();
      this.row_overlay = this.objects_ref.heatmap_line.clone();
      this.column_overlay = this.objects_ref.heatmap_line.clone();

      var self = this;

      this.heatmap_layer.on("mouseleave", function(evt){
          self.last_header = null;
          self.heatmap_overlay.destroyChildren();
          self.heatmap_overlay.draw();
          self.events.heatmap_onmouseout(evt);
      });

      if(this.settings.heatmap_header){
          this._draw_heatmap_header();
      }
  }

  InCHlib.prototype._draw_heatmap_row = function(node_id, x1, y1){
      var node = this.data.nodes[node_id];
      var row = new Kinetic.Group({id:node_id});
      var x2, y2, color, line, value, text, text_value, col_index;
      
      for (var i = 0, len = this.on_features["data"].length; i < len; i++){
          col_index = this.on_features["data"][i];
          x2 = x1 + this.pixels_for_dimension;
          y2 = y1;
          value = node.features[col_index];

          if(value !== null){
            color = this._get_color_for_value(value, this.data_descs[col_index]["min"], this.data_descs[col_index]["max"], this.data_descs[col_index]["middle"], this.settings.heatmap_colors);

            line = this.objects_ref.heatmap_line.clone({
                stroke: color,
                points: [x1, y1, x2, y2],
                value: value,
                column: ["d", col_index].join("_"),
                strokeWidth: this.pixels_for_leaf,
            });
            row.add(line);

            if(this.current_draw_values){
                text = this.objects_ref.heatmap_value.clone({
                    x: this._hack_round((x1 + x2)/2-(""+value).length*(this.value_font_size/4)),
                    y: this._hack_round(y1-this.value_font_size/2),
                    fontSize: this.value_font_size,
                    text: value,
                });
                row.add(text);
            }
          }
          x1 = x2;
      }

      if(this.settings.metadata){
          var metadata = this.metadata.nodes[node_id];

          if(metadata !== undefined){
            for (var i = 0, len = this.on_features["metadata"].length; i < len; i++){
                col_index = this.on_features["metadata"][i];
                value = metadata[col_index];
                x2 = x1 + this.pixels_for_dimension;
                y2 = y1;

                if(value !== null && value !== undefined){
                  text_value = value;
                  
                  if(this.metadata_descs[col_index]["str2num"] !== undefined){
                      value = this.metadata_descs[col_index]["str2num"][value];
                  }
                  color = this._get_color_for_value(value, this.metadata_descs[col_index]["min"], this.metadata_descs[col_index]["max"], this.metadata_descs[col_index]["middle"], this.settings.metadata_colors);
                      
                  line = this.objects_ref.heatmap_line.clone({
                          stroke: color,
                          points: [x1, y1, x2, y2],
                          value: text_value,
                          column: ["m", col_index].join("_"),
                          strokeWidth: this.pixels_for_leaf,
                      });
                  row.add(line);

                  if(this.current_draw_values){
                      text = this.objects_ref.heatmap_value.clone({
                          text: text_value,
                          fontSize: this.value_font_size,
                      });

                      width = text.getWidth();
                      x = this._hack_round((x1+x2)/2-width/2);
                      y = this._hack_round(y1-this.value_font_size/2);
                      text.position({x:x, y:y});
                      row.add(text);
                  }
                }
                x1 = x2;
            }
          }
      }

      if(this.settings.count_column && this.features[this.dimensions["overall"]-1] == 1){
          x2 = x1 + this.pixels_for_dimension;
          var count = node.objects.length;
          color = this._get_color_for_value(count, this.min_item_count, this.max_item_count, this.middle_item_count, this.settings.count_column_colors);
          
          line = this.objects_ref.heatmap_line.clone({
                  stroke: color,
                  points: [x1, y1, x2, y2],
                  value: count,
                  column: "Count",
                  strokeWidth: this.pixels_for_leaf,
          });
          row.add(line);

          if(this.current_draw_values){
              text = this.objects_ref.heatmap_value.clone({
                  text: count,
              });

              width = text.getWidth();
              x = this._hack_round((x1+x2)/2-width/2);
              y = this._hack_round(y1-this.value_font_size/2);
              text.position({x:x, y:y});
              row.add(text);
          }
          x1 = x2;
      }
      return row;
  }

  InCHlib.prototype._draw_column_metadata_row = function(data, row_index, x1, y1){
      var row = new Kinetic.Group({"class": "column_metadata"});
      var x2, y2, color, line, value, text, text_value, width, col_index;
      var str2num = (this.column_metadata_descs[row_index]["str2num"] === undefined)?false:true;

      for (var i = 0, len = this.on_features["data"].length; i < len; i++){
          col_index = this.on_features["data"][i];
          value = data[col_index];
          text_value = value;
          
          if(str2num){
              value = this.column_metadata_descs[row_index]["str2num"][value];
          }

          color = this._get_color_for_value(value, this.column_metadata_descs[row_index]["min"], this.column_metadata_descs[row_index]["max"], this.column_metadata_descs[row_index]["middle"], this.settings.column_metadata_colors);
          x2 = x1 + this.pixels_for_dimension;
          y2 = y1;
              
          line = this.objects_ref.heatmap_line.clone({
                  strokeWidth: this.settings.column_metadata_row_height,
                  stroke: color,
                  value: text_value,
                  points: [x1, y1, x2, y2],
                  column: ["d", col_index].join("_"),
              });
          row.add(line);
          x1 = x2;
      }
      return row;
  }

  InCHlib.prototype._bind_row_events = function(row){
      var self = this;
      row.on("mouseenter", function(evt){
          self._row_mouseenter(evt);
      });

      row.on("mouseleave", function(evt){
          self._row_mouseleave(evt);
      });

      row.on("mouseover", function(evt){
          self._draw_col_label(evt);
      });

      row.on("mouseout", function(evt){
          self.heatmap_overlay.find("#col_label")[0].destroy();
      });

      row.on("click", function(evt){
          var row_id = evt.target.parent.attrs.id;
          if(evt.target.parent.attrs.class !== "column_metadata"){
              var items = self.data.nodes[row_id].objects;
              var item_ids = [];
              
              for(i = 0; i < items.length; i++){
                  item_ids.push(items[i]);
              }
              self.events.row_onclick(item_ids, evt);
          }
      });
  }

  InCHlib.prototype._draw_row_ids = function(leaves_y){
      if(this.pixels_for_leaf < 6){
          return;
      }
      var i, objects, object_y = [], leaf, values = [], text;
      
      for(i = 0; i < leaves_y.length; i++){
          leaf = leaves_y[i];
          objects = this.data.nodes[leaf[0]].objects;
          if(objects.length > 1){
              return;
          }
          values.push(objects[0]);
          object_y.push([objects[0], leaf[1]]);
      }
      var max_length = this._get_max_length(values);
      var font_size = this._get_font_size(max_length, 85, this.pixels_for_leaf, 10);
      var x = this.distance + (this.on_features["data"].length + this.on_features["metadata"].length)*this.pixels_for_dimension + 15;
      x = (this.settings.count_column)?x+this.pixels_for_dimension: x;
      
      if(font_size > 4){
          for(i = 0; i < object_y.length; i++){
              text = this.objects_ref.heatmap_value.clone({
                  x: x,
                  y: this._hack_round(object_y[i][1] - font_size/2),
                  fontSize: font_size,
                  text: object_y[i][0],
                  fontStyle: 'italic',
                  fill: "gray"
              });
              this.heatmap_layer.add(text);
          }
      }
      return;
  }

  InCHlib.prototype._draw_heatmap_header = function(){
      this.header_layer = new Kinetic.Layer();
      var count = this._hack_size(this.leaves_y_coordinates);
      var y = (this.settings.column_dendrogram && this.heatmap_header)? this.header_height+(this.pixels_for_leaf*count) + 10: this.header_height - 20;
      var rotation = (this.settings.column_dendrogram && this.heatmap_header) ? 45 : -45;
      var distance_step = 0;
      var x, i, column_header, key;
      var current_headers = [];
      
      for(i = 0, len = this.on_features["data"].length; i < len; i++){
        current_headers.push(this.header[this.on_features["data"][i]]);
      }

      for(i = 0, len = this.on_features["metadata"].length; i < len; i++){
        current_headers.push(this.header[this.on_features["metadata"][i] + this.dimensions["data"]]);
      }

      var max_text_length = this._get_max_length(current_headers);
      var font_size = this._get_font_size(max_text_length, this.header_height, this.pixels_for_dimension, 16);
      if(font_size < 8){
          return;
      }
      
      for(i = 0, len = current_headers.length; i<len; i++){
        x = this.heatmap_distance+distance_step*this.pixels_for_dimension+this.pixels_for_dimension/2;
        column_header = this.objects_ref.column_header.clone({
                x: x,
                y: y,
                text: current_headers[i],
                position_index: i+1,
                fontSize: font_size,
                rotationDeg: rotation,
        });
        this.header_layer.add(column_header);
        distance_step++;
      }

      this.stage.add(this.header_layer);

      if(!(this.settings.dendrogram)){
          var self = this;

          this.header_layer.on("click", function(evt){
              var column = evt.target;
              var order = column.attrs.order;
              var position_index = column.attrs.position_index;

              for(i = 0; i<self.header_layer.getChildren().length; i++){
                  self.header_layer.getChildren()[i].setFill("black");
              }

              evt.target.setAttrs({"order":order, "fill": "red"});

              self.heatmap_layer.destroyChildren();
              self.heatmap_layer.draw();
              self._reorder_heatmap(position_index, order);
              self._draw_heatmap();
          })

          this.header_layer.on("mouseover", function(evt){
              var label = evt.target;
              label.setOpacity(0.7);
              this.draw();
          })

          this.header_layer.on("mouseout", function(evt){
              var label = evt.target;
              label.setOpacity(1);
              this.draw();
          })
      }
  }

  InCHlib.prototype._draw_distance_scale = function(distance){
      var y1 = this.header_height + this.column_metadata_height + this.settings.column_metadata_row_height/2 -10;
      var y2 = y1;
      var x1 = 0;
      var x2 = this.distance;
      var path = new Kinetic.Line({
          points: [x1, y1, x2, y2],
          stroke: "black",
          listening: false,
      })

      var circle = new Kinetic.Circle({
          x: x2, 
          y: y2,
          radius: 3,
          fill: "black",
          listening: false,
      })

      var number = 0;
      var marker_tail = 3;
      var marker_distance = x2;
      var marker_number_distance = this._hack_round(30/this.distance_step*10)/10;
      var distance = Math.round(100*this.distance/this.distance_step)/100;
      var marker_distance_step = this._hack_round(this.distance_step*marker_number_distance);
      var marker_counter = 0;

      var distance_number = new Kinetic.Text({
              x: 0,
              y: y1-20,
              text: distance,
              fontSize: 12,
              fontFamily: this.settings.font,
              fontStyle: 'bold',
              fill: 'black',
              align: 'right',
              listening: false,
      });
      this.dendrogram_layer.add(path, circle, distance_number);

      if(marker_distance_step==0){
          marker_distance_step=0.5;
      }

      var path;
      if(marker_number_distance > 0.1){
          while(marker_distance > 0){
              path = new Kinetic.Line({
                  points: [marker_distance, (y1-marker_tail), marker_distance, (y2+marker_tail)],
                  stroke: "black",
                  listening: false,
              })
              this.dendrogram_layer.add(path);

              number = this._hack_round((number + marker_number_distance)*10)/10;
              if(number>10){
                  number = this._hack_round(number);
              }
              
              marker_distance = marker_distance - marker_distance_step;
              marker_counter++;
          }
      }
  }

  InCHlib.prototype._draw_navigation = function(){
      this.navigation_layer = new Kinetic.Layer();
      var self = this;
      var x = 0;
      var y = 10;

      this._draw_color_scale();

      if(!this.settings.column_dendrogram){
          var filter_icon = this.objects_ref.icon.clone({
                  data: "M26.834,6.958c0-2.094-4.852-3.791-10.834-3.791c-5.983,0-10.833,1.697-10.833,3.791c0,0.429,0.213,0.84,0.588,1.224l8.662,15.002v4.899c0,0.414,0.709,0.75,1.583,0.75c0.875,0,1.584-0.336,1.584-0.75v-4.816l8.715-15.093h-0.045C26.625,7.792,26.834,7.384,26.834,6.958zM16,9.75c-6.363,0-9.833-1.845-9.833-2.792S9.637,4.167,16,4.167c6.363,0,9.834,1.844,9.834,2.791S22.363,9.75,16,9.75z",
                  x: x,
                  y: y,
                  label: "Filter\ncolumns"
          });
      
          var filter_overlay = self._draw_icon_overlay(x, y);
          self.navigation_layer.add(filter_icon, filter_overlay);
          x = x + 40;    
      
          filter_overlay.on("click", function(){
              self._filter_icon_click(this);
          });
              
          filter_overlay.on("mouseover", function(){
              self._icon_mouseover(filter_icon, filter_overlay, self.navigation_layer);
          });
             
          filter_overlay.on("mouseout", function(){
              self._icon_mouseout(filter_icon, filter_overlay, self.navigation_layer);
          });
      }

      if(self.zoomed_clusters["row"].length > 0){
          var unzoom_icon = this.objects_ref.icon.clone({
              data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z",
              x: x,
              y: y,
              label: "Unzoom\nrows"
          });
          var unzoom_overlay = self._draw_icon_overlay(x, y);
          x = x + 40;

          self.navigation_layer.add(unzoom_icon, unzoom_overlay);

          unzoom_overlay.on("click", function(){
              self.events.on_unzoom(self._unprefix(self.zoomed_clusters["row"][self.zoomed_clusters["row"].length-1]));
              self._unzoom_icon_click(this);
          });

          unzoom_overlay.on("mouseover", function(){
              self._icon_mouseover(unzoom_icon, unzoom_overlay, self.navigation_layer);
          });

          unzoom_overlay.on("mouseout", function(){
              self._icon_mouseout(unzoom_icon, unzoom_overlay, self.navigation_layer);
          });
      }

      if(self.zoomed_clusters["row"].length > 0 || self.zoomed_clusters["column"].length > 0){
        var refresh_icon = this.objects_ref.icon.clone({
              data: "M24.083,15.5c-0.009,4.739-3.844,8.574-8.583,8.583c-4.741-0.009-8.577-3.844-8.585-8.583c0.008-4.741,3.844-8.577,8.585-8.585c1.913,0,3.665,0.629,5.09,1.686l-1.782,1.783l8.429,2.256l-2.26-8.427l-1.89,1.89c-2.072-1.677-4.717-2.688-7.587-2.688C8.826,3.418,3.418,8.826,3.416,15.5C3.418,22.175,8.826,27.583,15.5,27.583S27.583,22.175,27.583,15.5H24.083z",
              x: x,
              y: y,
              id: "refresh_icon",
              label: "Refresh"
        });
        var refresh_overlay = self._draw_icon_overlay(x, y);
        self.navigation_layer.add(refresh_icon, refresh_overlay);

        refresh_overlay.on("click", function(){
            self._refresh_icon_click();
            self.events.on_refresh();
        });

        refresh_overlay.on("mouseover", function(){
            self._icon_mouseover(refresh_icon, refresh_overlay, self.navigation_layer);
        });

        refresh_overlay.on("mouseout", function(){
            self._icon_mouseout(refresh_icon, refresh_overlay, self.navigation_layer);
        });
      }

      if(self.zoomed_clusters["column"].length > 0){
          x = self.settings.width - 60;
          y = 45;
          var column_unzoom_icon = this.objects_ref.icon.clone({
              data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z",
              x: x,
              y: y,
              label: "Unzoom\ncolumns"
          });
          var column_unzoom_overlay = self._draw_icon_overlay(x, y);

          self.navigation_layer.add(column_unzoom_icon, column_unzoom_overlay);

          column_unzoom_overlay.on("click", function(){
              // self.events.on_unzoom(self._unprefix(self.zoomed_clusters["row"][self.zoomed_clusters["row"].length-1]));
              self._column_unzoom_icon_click(this);
          });

          column_unzoom_overlay.on("mouseover", function(){
              self._icon_mouseover(column_unzoom_icon, column_unzoom_overlay, self.navigation_layer);
          });

          column_unzoom_overlay.on("mouseout", function(){
              self._icon_mouseout(column_unzoom_icon, column_unzoom_overlay, self.navigation_layer);
          });
      }

      if(self.settings.show_export_button){
        var export_icon = this.objects_ref.icon.clone({
              data: "M24.25,10.25H20.5v-1.5h-9.375v1.5h-3.75c-1.104,0-2,0.896-2,2v10.375c0,1.104,0.896,2,2,2H24.25c1.104,0,2-0.896,2-2V12.25C26.25,11.146,25.354,10.25,24.25,10.25zM15.812,23.499c-3.342,0-6.06-2.719-6.06-6.061c0-3.342,2.718-6.062,6.06-6.062s6.062,2.72,6.062,6.062C21.874,20.78,19.153,23.499,15.812,23.499zM15.812,13.375c-2.244,0-4.062,1.819-4.062,4.062c0,2.244,1.819,4.062,4.062,4.062c2.244,0,4.062-1.818,4.062-4.062C19.875,15.194,18.057,13.375,15.812,13.375z",
              x: self.settings.width - 62,
              y: 10,
              id: "export_icon",
              label: "Export\nin png format"
        });

        var export_overlay = self._draw_icon_overlay(self.settings.width - 62, 10);
        self.navigation_layer.add(export_icon, export_overlay);

        export_overlay.on("click", function(){
            self._export_icon_click(this);
        });
            
        export_overlay.on("mouseover", function(){
            self._icon_mouseover(export_icon, export_overlay, self.navigation_layer);
        });
           
        export_overlay.on("mouseout", function(){
            self._icon_mouseout(export_icon, export_overlay, self.navigation_layer);
        });
      }

      if(self.settings.show_settings_button){
        var settings_icon = this.objects_ref.icon.clone({
              data: "M26.834,14.693c1.816-2.088,2.181-4.938,1.193-7.334l-3.646,4.252l-3.594-0.699L19.596,7.45l3.637-4.242c-2.502-0.63-5.258,0.13-7.066,2.21c-1.907,2.193-2.219,5.229-1.039,7.693L5.624,24.04c-1.011,1.162-0.888,2.924,0.274,3.935c1.162,1.01,2.924,0.888,3.935-0.274l9.493-10.918C21.939,17.625,24.918,16.896,26.834,14.693z",
              x: self.settings.width - 62,
              y: 80,
              id: "settings_icon",
              label: "Settings"
        });

        var settings_overlay = self._draw_icon_overlay(self.settings.width - 62, 80);
        self.navigation_layer.add(settings_icon, settings_overlay);

        settings_overlay.on("click", function(evt){
            self._settings_icon_click(this, evt);
        });
            
        settings_overlay.on("mouseover", function(){
            self._icon_mouseover(settings_icon, settings_overlay, self.navigation_layer);
        });
           
        settings_overlay.on("mouseout", function(){
            self._icon_mouseout(settings_icon, settings_overlay, self.navigation_layer);
        });
      }

      self.stage.add(self.navigation_layer);
  }

  InCHlib.prototype._draw_color_scale = function(){
      var self = this;
      var color_steps = [this.settings.min_quantile/100, this._get_color_for_value(0, 0, 1, 0.5, this.settings.heatmap_colors), this.settings.middle_quantile/100, this._get_color_for_value(0.5, 0, 1, 0.5, this.settings.heatmap_colors), this.settings.max_quantile/100, this._get_color_for_value(1, 0, 1, 0.5, this.settings.heatmap_colors)];
      var color_scale = this.objects_ref.rect_gradient.clone({"label": "Color settings",
                                                              "fillLinearGradientColorStops": color_steps,
                                                              "id": "color_scale"});

      color_scale.on("mouseover", function(){
        self._color_scale_mouseover(color_scale, self.navigation_layer);
      });

      color_scale.on("mouseout", function(){
        self._color_scale_mouseout(color_scale, self.navigation_layer);
      });

      color_scale.on("click", function(){
        self._color_scale_click(color_scale, self.navigation_layer);
      });

      this.navigation_layer.add(color_scale);
  }

  InCHlib.prototype._update_color_scale = function(){
    var color_scale = this.navigation_layer.find("#color_scale");
    color_scale.fillLinearGradientColorStops([this.settings.min_quantile/100, this._get_color_for_value(0, 0, 1, 0.5, this.settings.heatmap_colors), this.settings.middle_quantile/100, this._get_color_for_value(0.5, 0, 1, 0.5, this.settings.heatmap_colors), this.settings.max_quantile/100, this._get_color_for_value(1, 0, 1, 0.5, this.settings.heatmap_colors)]);
    this.navigation_layer.draw();
  }

  InCHlib.prototype._draw_icon_overlay = function(x, y){
      var icon_overlay = this.objects_ref.icon_overlay.clone({x: x, y: y});
      return icon_overlay;
  }

  InCHlib.prototype._highlight_path = function(path_id){
      var node = this.data.nodes[path_id];
      if(node.count != 1){
          var path = this.dendrogram_layer.get("#"+path_id)[0].clone({"stroke": "#F5273C"});
          var rect = this.dendrogram_layer.get("#"+path_id+"_rect")[0].clone({path: path, path_id: path_id});            
          this.dendrogram_overlay_layer.add(path, rect);
          this._highlight_path(node.left_child);
          this._highlight_path(node.right_child);
      }
      else{
          this.highlighted_rows_y.push(this.leaves_y_coordinates[path_id]);
       }
   }

   InCHlib.prototype._highlight_column_path = function(path_id){
      var node = this.column_dendrogram.nodes[path_id];
      if(node.count != 1){
          var path = this.column_dendrogram_layer.get("#col"+path_id)[0].clone({"stroke": "#F5273C"});
          var rect = this.column_dendrogram_layer.get("#col_rect"+path_id)[0].clone({path: path, path_id: path_id});            
          this.column_dendrogram_overlay_layer.add(path, rect);
          this._highlight_column_path(node.left_child);
          this._highlight_column_path(node.right_child);
      }
   }

  /**
    * Unhighlight highlighted heatmap rows. 
    *
    * @example 
    * instance.unhighlight_rows();
    */
   InCHlib.prototype.unhighlight_rows = function(){
      this.highlight_rows([]);
      return;
   }

  /**
    * Highlight heatmap rows with color defined in instance.settings.highlight_colors. 
    * When the empty array is passed it unhighlights all highlighted rows.
    * 
    * @param {object} [row_ids] The array of heatmap row (object) IDs.
    *
    * @example 
    * instance.highlight_rows(["CHEMBL7781", "CHEMBL273658", "CHEMBL415309", "CHEMBL267231", "CHEMBL8007", "CHEMBL7987", "CHEMBL7988", "CHEMBL266282", "CHEMBL7655", "CHEMBL7817", "CHEMBL8637", "CHEMBL8639", "CHEMBL8055", "CHEMBL7843", "CHEMBL266488", "CHEMBL8329"]);
    */

   InCHlib.prototype.highlight_rows = function(row_ids){
      var i, row, row_id;
      if(!this.settings.heatmap){
        return;
      }

      this.settings.highlighted_rows = row_ids;
      this.highlighted_rows_layer.destroyChildren();

      var original_colors = this.settings.heatmap_colors;
      var original_metadata_colors = this.settings.metadata_colors;
      this.settings.heatmap_colors = this.settings.highlight_colors;
      this.settings.metadata_colors = this.settings.highlight_colors;

      var done_rows = {};
      var unique_row_ids = [];

      for(i = 0; i<row_ids.length; i++){
          if(this.objects2leaves[row_ids[i]] !== undefined){
              row_id = this.objects2leaves[row_ids[i]];
              if(done_rows[row_id] === undefined){
                  unique_row_ids.push(row_id);
                  done_rows[row_id] = null;
              }
          }
      }

      for(i = 0; i<unique_row_ids.length; i++){
          row = this._draw_heatmap_row(unique_row_ids[i], this.heatmap_distance, this.leaves_y_coordinates[unique_row_ids[i]]);
          this.highlighted_rows_layer.add(row);
          row.setAttr("listening", false);
      }


      this.highlighted_rows_layer.draw();
      this.heatmap_overlay.moveToTop();

      this.settings.heatmap_colors = original_colors;
      this.settings.metadata_colors = original_metadata_colors;

      var self = this;

      this.highlighted_rows_layer.on("click", function(evt){
          self.heatmap_layer.fire("click");
      });
      
   }

  InCHlib.prototype._highlight_cluster = function(path_id){
      this.current_object_ids = [];

      if(this.last_highlighted_cluster){
        if(this.settings.heatmap){
            this.row_cluster_group.destroy();
        }
        this.dendrogram_overlay_layer.destroyChildren();
      }

      this.highlighted_rows_y = [];
      this.last_highlighted_cluster = path_id;
      this._highlight_path(path_id);
      this.dendrogram_overlay_layer.draw();
      this._draw_cluster_layer(path_id);
      this._get_object_ids(path_id);
      this.events.dendrogram_node_highlight(this.current_object_ids, path_id);
      return;
  }

  InCHlib.prototype._highlight_column_cluster = function(path_id){
      if(this.last_highlighted_column_cluster){
        this.column_cluster_group.destroy();
        this.column_dendrogram_overlay_layer.destroyChildren();
      }
      this.last_highlighted_column_cluster = path_id;
      this._highlight_column_path(path_id);
      this.column_dendrogram_overlay_layer.draw();
      this._get_column_ids(path_id);
      this.current_column_ids.sort(function(a,b){return a - b});
      this._draw_column_cluster_layer(path_id);
      // this.events.dendrogram_node_highlight(this.current_object_ids, path_id);
      return;
  }

  InCHlib.prototype.unhighlight_column_cluster = function(){
      if(this.last_highlighted_column_cluster){
        this.current_column_ids = [];
        this.column_cluster_group.destroy();
        this.cluster_layer.draw(); 
        this.column_dendrogram_overlay_layer.destroyChildren();
        this.column_dendrogram_overlay_layer.draw();
        this.last_highlighted_column_cluster = null;
      }
      return [];
  }

  /**
    * Highlight cluster defined by the dendrogram node ID.
    * 
    * @param {string} node_id The ID of particular node in dendrogram.
    *
    * @example 
    * instance.highlight_cluster("node@715");
    */

  InCHlib.prototype.highlight_cluster = function(node_id){
      return this._highlight_cluster(this._prefix(node_id));
  }

  /**
    * Unhighlight highlighted dendrogram node (cluster).
    *
    * @example 
    * instance.unhighlight_cluster();
    */
  InCHlib.prototype.unhighlight_cluster = function(){
      if(this.last_highlighted_cluster){
        if(this.settings.heatmap){
            this.row_cluster_group.destroy();
            this.cluster_layer.draw();
          }
          this.dendrogram_overlay_layer.removeChildren();
          this.dendrogram_overlay_layer.draw();
          this.events.dendrogram_node_unhighlight(this.last_highlighted_cluster);
          this.last_highlighted_cluster = null;
      }
      return [];
  }

  InCHlib.prototype._neutralize_path = function(path_id){
      var node = this.data.nodes[path_id];

      if(node.count != 1){
          var path = this.dendrogram_layer.get("#"+path_id)[0];
          if(path){
              path.setStroke("grey");
              this._neutralize_path(node.right_child);
              this._neutralize_path(node.left_child);
          }
      }
  }

  InCHlib.prototype._draw_cluster_layer = function(path_id){
      this.row_cluster_group = new Kinetic.Group();
      var visible = this._get_visible_count();
      var count = this.data.nodes[path_id].count;
      var x = (this.settings.column_dendrogram)?0:40;
      var y = 10;

      var rows_desc = this.objects_ref.count.clone({x: x + 25,
                                                        y: y - 5,
                                                        text: count,
                                                        fontSize: 12,
                                                        fill: "#6d6b6a"
                                                    });

      var zoom_icon = this.objects_ref.icon.clone({
                      data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM15.687,9.051h-4v2.833H8.854v4.001h2.833v2.833h4v-2.834h2.832v-3.999h-2.833V9.051z",
                      x: x,
                      y: y,
                      label: "Zoom\nrows",
                  });

      var zoom_overlay = this._draw_icon_overlay(x, y);

      x = this.heatmap_distance;
      var width = visible*this.pixels_for_dimension;
      var upper_y = this.highlighted_rows_y[0]-this.pixels_for_leaf/2;
      var lower_y = this.highlighted_rows_y[this.highlighted_rows_y.length-1]+this.pixels_for_leaf/2;
      
      var cluster_overlay_1 = this.objects_ref.cluster_overlay.clone({
          x: x,
          y: this.header_height + this.column_metadata_height,
          width: width,
          height: upper_y-this.header_height - this.column_metadata_height,
      });

      var cluster_border_1 = this.objects_ref.cluster_border.clone({
          points: [x, upper_y, x+width, upper_y],
      });

      var cluster_overlay_2 = this.objects_ref.cluster_overlay.clone({
          x: x,
          y: lower_y,
          width: width,
          height: this.settings.height-lower_y-this.footer_height + 0.5*this.pixels_for_leaf,
      });

      var cluster_border_2 = this.objects_ref.cluster_border.clone({
          points: [x, lower_y, x+width, lower_y],
      });

      this.row_cluster_group.add(rows_desc, cluster_overlay_1, cluster_overlay_2, zoom_icon, zoom_overlay, cluster_border_1, cluster_border_2);
      this.cluster_layer.add(this.row_cluster_group);
      this.stage.add(this.cluster_layer);
      rows_desc.moveToTop();

      this.cluster_layer.draw();

      var self = this;

      zoom_overlay.on("mouseover", function(){
          self._icon_mouseover(zoom_icon, zoom_overlay, self.cluster_layer);
      });

      zoom_overlay.on("mouseout", function(){
          self._icon_mouseout(zoom_icon, zoom_overlay, self.cluster_layer);
      });

      zoom_overlay.on("click", function(){
          self._zoom_cluster(self.last_highlighted_cluster);
      });

      this.cluster_layer.on("click", function(evt){
          self.unhighlight_cluster();
          self.unhighlight_column_cluster();
          self.events.empty_space_onclick(evt);
      });
  }

  InCHlib.prototype._draw_column_cluster_layer = function(){
      this.column_cluster_group = new Kinetic.Group();
      var x = this.settings.width - 60;
      var y = 45;

      var cols_desc = this.objects_ref.count.clone({x: x + 25,
                                                        y: y,
                                                        text: this.current_column_ids.length,
                                                        fontSize: 12,
                                                        fill: "#6d6b6a"
                                                    });

      var zoom_icon = this.objects_ref.icon.clone({
                      data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM15.687,9.051h-4v2.833H8.854v4.001h2.833v2.833h4v-2.834h2.832v-3.999h-2.833V9.051z",
                      x: x,
                      y: y,
                      label: "Zoom\ncolumns",
                  });

      var zoom_overlay = this._draw_icon_overlay(x, y);

      var x1 = this._hack_round((this.current_column_ids[0] - this.columns_start_index)*this.pixels_for_dimension);
      var x2 = this._hack_round((this.current_column_ids[0] + this.current_column_ids.length - this.columns_start_index)*this.pixels_for_dimension);
      var y1 = this.header_height;
      var y2 = this.settings.height-this.footer_height+0.5*this.pixels_for_leaf;
      var height = this.settings.height-this.footer_height-this.header_height+this.settings.column_metadata_row_height;    
      
      var cluster_border_1 = this.objects_ref.cluster_border.clone({
          points: [this.heatmap_distance + x1, y1, this.heatmap_distance + x1, y2],
      });

      var cluster_overlay_1 = this.objects_ref.cluster_overlay.clone({
          x: this.heatmap_distance,
          y: this.header_height,
          width: x1,
          height: height,
      });

      var cluster_border_2 = this.objects_ref.cluster_border.clone({
          points: [this.heatmap_distance + x2, y1, this.heatmap_distance + x2, y2],
      });

      var cluster_overlay_2 = this.objects_ref.cluster_overlay.clone({
          x: x2+this.heatmap_distance,
          y: this.header_height,
          width: this.heatmap_width - x2,
          height: height,
      });
    

      this.column_cluster_group.add(cluster_overlay_1, cluster_overlay_2, zoom_icon, zoom_overlay, cols_desc, cluster_border_1, cluster_border_2);
      this.cluster_layer.add(this.column_cluster_group);
      this.stage.add(this.cluster_layer);
      this.cluster_layer.draw();

      var self = this;

      zoom_overlay.on("mouseover", function(){
          self._icon_mouseover(zoom_icon, zoom_overlay, self.cluster_layer);
      });

      zoom_overlay.on("mouseout", function(){
          self._icon_mouseout(zoom_icon, zoom_overlay, self.cluster_layer);
      });

      zoom_overlay.on("click", function(){
          self._zoom_column_cluster(self.last_highlighted_column_cluster);
      });

      this.cluster_layer.on("click", function(evt){
          self.unhighlight_cluster();
          self.unhighlight_column_cluster();
          self.events.empty_space_onclick(evt);
      });
  }

  InCHlib.prototype._zoom_column_cluster = function(node_id){
    if(node_id != this.zoomed_clusters["column"][this.zoomed_clusters["column"].length-1]){
      if(node_id !== this.column_root_id){
        this.zoomed_clusters["column"].push(node_id);
      }
      
      this._get_column_ids(node_id);
      this.current_column_ids.sort(function(a,b){return a - b});
      this.columns_start_index = this.current_column_ids[0];
      this.on_features["data"] = this.current_column_ids;
      var distance = this.distance;
      this._adjust_dimension_size(this.on_features["data"].length + this.on_features["metadata"].length);
      this._delete_layers([this.column_dendrogram_layer, this.column_dendrogram_overlay_layer, this.heatmap_layer, this.cluster_layer, this.navigation_layer], [this.dendrogram_hover_layer]);
      if(this.settings.heatmap_header){
        this._delete_layers([this.header_layer]);
      }
      this._draw_column_dendrogram(node_id);
      this._draw_heatmap();
      this._draw_navigation();

      var highlighted_row_node = this.last_highlighted_cluster;

      if(distance !== this.distance){
        this._delete_layers([this.dendrogram_layer, this.dendrogram_overlay_layer]);
        var row_node = (this.zoomed_clusters["row"].length > 0)?this.zoomed_clusters["row"][this.zoomed_clusters["row"].length - 1]:this.root_id;
        this._draw_row_dendrogram(row_node);
      }
      // this.highlight_rows(this.settings.highlighted_rows);
      // this.events.on_zoom(this._unprefix(node_id));

      if(highlighted_row_node !== null){
        this._highlight_cluster(highlighted_row_node);
      }
    }
  }

  InCHlib.prototype._zoom_cluster = function(node_id){
      if(node_id != this.zoomed_clusters["row"][this.zoomed_clusters["row"].length-1]){
          if(node_id !== this.root_id){
            this.zoomed_clusters["row"].push(node_id);
          }
          this._delete_layers([this.dendrogram_layer, this.heatmap_layer, this.dendrogram_overlay_layer, this.cluster_layer, this.navigation_layer, this.header_layer], [this.dendrogram_hover_layer]);
          this._draw_row_dendrogram(node_id);
          this._draw_heatmap();
          this._draw_navigation();
          if(this.settings.column_dendrogram && this.last_highlighted_column_cluster !== null){
            this._draw_column_cluster_layer(this.last_highlighted_column_cluster);
          }
          this.highlight_rows(this.settings.highlighted_rows);
          this.events.on_zoom(this._unprefix(node_id));
      }
  }

  InCHlib.prototype._get_node_neighbourhood = function(node, nodes){
      var self = this;
      var node_neighbourhood = {"left_node": {"left_node": {"left_count" : 0,
                                                            "right_count": 0}, 
                                              "right_node": {"left_count" : 0,
                                                             "right_count": 0},
                                              "left_count" : 0.5,
                                              "right_count": 0.5
                                             },
                                "right_node": {"left_node": {"left_count" : 0,
                                                            "right_count": 0}, 
                                              "right_node": {"left_count" : 0,
                                                             "right_count": 0},
                                              "left_count" : 0.5,
                                              "right_count": 0.5
                                             },
                                "left_count": nodes[node.left_child].count,
                                "right_count": nodes[node.right_child].count,
      };

      var left_child = nodes[node.left_child];
      var right_child = nodes[node.right_child];

      var left_child_left_child = nodes[left_child.left_child];
      var left_child_right_child = nodes[left_child.right_child];

      var right_child_left_child = nodes[right_child.left_child];
      var right_child_right_child = nodes[right_child.right_child];

      if(left_child.count != 1){
              node_neighbourhood.left_node.left_count = nodes[left_child.left_child].count;
              node_neighbourhood.left_node.right_count = nodes[left_child.right_child].count;

          if(left_child_left_child.count != 1){
              node_neighbourhood.left_node.left_node.left_count = nodes[left_child_left_child.left_child].count;
              node_neighbourhood.left_node.left_node.right_count = nodes[left_child_left_child.right_child].count;
          }
          else{
              node_neighbourhood.left_node.left_node.left_count = 0.5;
              node_neighbourhood.left_node.left_node.right_count = 0.5;
          }

          if(left_child_right_child.count != 1){
              node_neighbourhood.left_node.right_node.left_count = nodes[left_child_right_child.left_child].count;
              node_neighbourhood.left_node.right_node.right_count = nodes[left_child_right_child.right_child].count;
          }
          else{
              node_neighbourhood.left_node.right_node.left_count = 0.5;
              node_neighbourhood.left_node.right_node.right_count = 0.5;
          }
      }

      if(right_child.count != 1){
          node_neighbourhood.right_node.left_count = nodes[right_child.left_child].count;
          node_neighbourhood.right_node.right_count = nodes[right_child.right_child].count;

          if(right_child_left_child.count != 1){
              node_neighbourhood.right_node.left_node.left_count = nodes[right_child_left_child.left_child].count;
              node_neighbourhood.right_node.left_node.right_count = nodes[right_child_left_child.right_child].count;
          }
          else{
              node_neighbourhood.right_node.left_node.left_count = 0.5;
              node_neighbourhood.right_node.left_node.right_count = 0.5;
          }

          if(right_child_right_child.count != 1){
              node_neighbourhood.right_node.right_node.left_count = nodes[right_child_right_child.left_child].count;
              node_neighbourhood.right_node.right_node.right_count = nodes[right_child_right_child.right_child].count;
          }
          else{
              node_neighbourhood.right_node.right_node.left_count = 0.5;
              node_neighbourhood.right_node.right_node.right_count = 0.5;
          }
      }
      return node_neighbourhood;
  }

  InCHlib.prototype._draw_column_dendrogram_node = function(node_id, node, current_left_count, current_right_count, x, y){
      
      if(node.count > 1){
          var node_neighbourhood = this._get_node_neighbourhood(node, this.column_dendrogram.nodes);
          var right_child = this.column_dendrogram.nodes[node.right_child];
          var left_child = this.column_dendrogram.nodes[node.left_child];
          var x1 = this._get_x1(node_neighbourhood, current_left_count, current_right_count);
          var x2 = this._get_x2(node_neighbourhood, current_left_count, current_right_count);
          var y1 = this._hack_round(this.vertical_distance - this.vertical_distance_step*node.distance);
          y1 = (y1 == 0)? 2: y1;
          var y2 = y1;

          if(right_child.count == 1){
              x2 = x2 - this.pixels_for_dimension/2;
          }

          var left_distance = this.vertical_distance - this.vertical_distance_step*this.column_dendrogram.nodes[node.left_child].distance;
          var right_distance = this.vertical_distance - this.vertical_distance_step*this.column_dendrogram.nodes[node.right_child].distance;

          this.column_dendrogram_layer.add(this._draw_vertical_path(node_id, x1, y1, x2, y2, left_distance, right_distance));
          this._draw_column_dendrogram_node(node.left_child, left_child, current_left_count - node_neighbourhood.left_node.right_count, current_right_count + node_neighbourhood.left_node.right_count, left_distance, y1);
          this._draw_column_dendrogram_node(node.right_child, right_child, current_left_count + node_neighbourhood.right_node.left_count, current_right_count - node_neighbourhood.right_node.left_count, right_distance, y2);
      }
      else{
        this.column_x_coordinates[node_id] = current_right_count*this.pixels_for_dimension;
      }
  }

  InCHlib.prototype._get_y1 = function(node_neighbourhood, current_left_count, current_right_count){
      current_left_count = current_left_count-node_neighbourhood.left_node.right_count-node_neighbourhood.left_node.left_node.right_count;
      var y = (current_left_count+(node_neighbourhood.left_node.left_node.right_count+node_neighbourhood.left_node.right_node.left_count)/2)*this.pixels_for_leaf;
      return y + this.top_heatmap_distance;
  }

  InCHlib.prototype._get_y2 = function(node_neighbourhood, current_left_count, current_right_count){
      current_left_count = current_left_count+node_neighbourhood.right_node.left_node.left_count;
      var y = (current_left_count+(node_neighbourhood.right_node.left_node.right_count+node_neighbourhood.right_node.right_node.left_count)/2)*this.pixels_for_leaf;
      return y + this.top_heatmap_distance;
  }

  InCHlib.prototype._get_x1 = function(node_neighbourhood, current_left_count, current_right_count){
      current_left_count = current_left_count-node_neighbourhood.left_node.right_count-node_neighbourhood.left_node.left_node.right_count;
      var x = (current_left_count+(node_neighbourhood.left_node.left_node.right_count+node_neighbourhood.left_node.right_node.left_count)/2)*this.pixels_for_dimension;
      return (this.heatmap_distance+this.on_features["data"].length * this.pixels_for_dimension)-x;
  }

  InCHlib.prototype._get_x2 = function(node_neighbourhood, current_left_count, current_right_count){
      current_left_count = current_left_count+node_neighbourhood.right_node.left_node.left_count;
      var x = (current_left_count+(node_neighbourhood.right_node.left_node.right_count+node_neighbourhood.right_node.right_node.left_count)/2)*this.pixels_for_dimension;;
      return (this.heatmap_distance+this.on_features["data"].length * this.pixels_for_dimension)-x;
  }

  InCHlib.prototype._draw_vertical_path = function(path_id, x1, y1, x2, y2, left_distance, right_distance){
      var path_group = new Kinetic.Group({});
      var path = this.objects_ref.node.clone({points: [x1, left_distance, x1, y1, x2, y2, x2, right_distance], id: "col" + path_id,})
      var path_rect = this.objects_ref.node_rect.clone({x: x2-1,
                                                            y: y1-1,
                                                            width: x1 - x2 + 2,
                                                            height: this.header_height - y1,
                                                            id: "col_rect" + path_id,
                                                            path: path,
                                                            path_id: path_id,
                                                          });

      path_group.add(path, path_rect);
      return path_group;
  }

  InCHlib.prototype._draw_horizontal_path = function(path_id, x1, y1, x2, y2, left_distance, right_distance){
      var path_group = new Kinetic.Group({});
      var path = this.objects_ref.node.clone({points: [left_distance, y1, x1, y1, x2, y2, right_distance, y2],
                                                  id: path_id});

      var path_rect = this.objects_ref.node_rect.clone({x: x1-1,
                                                            y: y1-1,
                                                            width: this.distance - x1,
                                                            height: y2 - y1,
                                                            id: [path_id, "rect"].join("_"),
                                                            path: path,
                                                            path_id: path_id,
                                                          });
      path_group.add(path, path_rect);
      return path_group;
  }

  InCHlib.prototype._filter_icon_click = function(filter_button){
      var self = this;
      var filter_features_element = this.target_element.find(".filter_features");
      var symbol = "✖";

      if(filter_features_element.length){
          filter_features_element.fadeIn("fast");
          var overlay = self._draw_target_overlay();
      }
      else{
          filter_list = "";
          
          for(var attr in self.header){
              if(self.features[attr] == 1){
                  symbol = "✔";
              }
              if(attr < self.dimensions){
                  var text = self.header[attr];
                  if(text == ""){
                      text =  parseInt(attr) + 1 + ". column";
                  }
                  filter_list = filter_list + "<li class='feature_switch' data-num='" + attr + "'><span class='symbol'>" + symbol + "</span>  " + text +"</li>";
              }
          }
          
          self.target_element.append("<div class='filter_features'><ul>" + filter_list + "</ul><hr /><div><span class='cancel_filter_list'>Cancel</span>&nbsp;&nbsp;&nbsp;<span class='update_filter_list'>Update</span></div></div>");
          filter_features_element = self.target_element.find(".filter_features");
          
          filter_features_element.css({"display":"none",
              "top": 45,
              "left": 0,
              "border-radius":"5px",
              "text-align":"center",
              "position":"absolute",
              "background-color":"#ffffff",
              "border":"solid 2px #DEDEDE",
              "padding-top":"5px",
              "padding-left":"15px",
              "padding-bottom":"10px",
              "padding-right":"15px",
              "font-weight":"bold",
              "font-size": "14px",
              "z-index": 1000,
              "font-family": self.settings.font
          });

          filter_features_element.find("ul").css({
              "list-style-type":"none",
              "margin-left":"0",
              "padding-left":"0",
              "text-align":"left",
          });

          filter_features_element.find("li").css({
              "color":"green",
              "margin-top":"5px",
          });

          filter_features_element.find("div").css({
              "cursor":"pointer",
              "opacity":"0.7",
          });

          var overlay = self._draw_target_overlay();
          filter_features_element.fadeIn("fast");

          self.target_element.find(".feature_switch").click(function(){
              var num = parseInt($(this).attr("data-num"));
              var symbol_element = $(this).find("span");
              self.features[num] = -self.features[num];

              if(self.features[num] == 1){
                  symbol_element.text("✔");
                  $(this).css("color", "green");
              }
              else{
                  symbol_element.text("✖");
                  $(this).css("color", "red");
              }

              self._set_on_features();
          });

          $(function(){
              filter_features_element.click(function(){
                  return false;
              });

              filter_features_element.mousedown(function(){
                  return false;
              });

             $("#" + self.settings.target + " .filter_features ul li," + "#" + self.settings.target + " .filter_features div span").hover(
             function(){
                $(this).css({
                      "cursor": "pointer",
                      "opacity": "0.7",
                  });
             },
             function(){
                $(this).css({
                      "cursor": "default",
                      "opacity": "1",
                  });
             });
          });

          self.target_element.find(".cancel_filter_list").click(function(){
              filter_features_element.fadeOut("fast");
              overlay.fadeOut("fast");
          });

          overlay.click(function(){
              filter_features_element.fadeOut("fast");
              overlay.fadeOut("fast");
          });

          self.target_element.find(".update_filter_list").click(function(){
              filter_features_element.fadeOut("slow");
              overlay.fadeOut("slow");

              var node_id = (self.zoomed_clusters["row"].length > 0)?self.zoomed_clusters["row"][self.zoomed_clusters["row"].length-1]:self.root_id;
              var highlighted_cluster = self.last_highlighted_cluster;
              self.last_highlighted_cluster = null;
              self._set_on_features();
              var current_dimensions_count = self.on_features["data"].length+self.on_features["metadata"].length;
              self._adjust_dimension_size(current_dimensions_count);
              var curren_width = self.distance + self.heatmap_width;
              self.distance = curren_width - current_dimensions_count*self.pixels_for_dimension;
              self.heatmap_width = current_dimensions_count*self.pixels_for_dimension;

              
              self._delete_all_layers();
              self._draw_stage_layer();
              if(self.settings.dendrogram){
                  self._draw_row_dendrogram(node_id);
              }

              if(self._visible_features_equal_column_dendrogram_count() && self.settings.column_dendrogram){
                  self._draw_column_dendrogram(self.column_root_id);
              }

              self._draw_navigation();
              self._draw_heatmap();

              if(highlighted_cluster != null){
                  self._highlight_cluster(highlighted_cluster);
              }
          });
      }
  }
  InCHlib.prototype._draw_target_overlay = function(){
    var overlay = this.target_element.find(".target_overlay");

    if(overlay.length){
      overlay.fadeIn("fast");
    }
    else{
      overlay = $("<div class='target_overlay'></div>");
      overlay.css({"background-color": "white", 
                      "position": "absolute",
                      "top": 0,
                      "left": 0,
                      "right": 0,
                      "bottom": 0,
                      "opacity": 0.5
          });
      this.target_element.append(overlay);
    }

    return overlay;
  }

  InCHlib.prototype._refresh_icon_click = function(){
      var node_id = this.root_id;
      this.zoomed_clusters = {"row": [], "column": []};
      this._adjust_dimension_size(this.dimensions["overall"]);
      this._delete_all_layers();
      this._draw_stage_layer();
      this._set_on_features();

      this._draw_row_dendrogram(node_id);
      if(this.settings.column_dendrogram){
          this.columns_start_index = 0;
          this._get_column_ids(this.column_root_id);
          this._draw_column_dendrogram(this.column_root_id);
      }
      this._draw_navigation();
      this._draw_heatmap();
      this.highlight_rows(this.settings.highlighted_rows);
  }

  InCHlib.prototype._export_icon_click = function(){
    var self = this;
    var export_menu = this.target_element.find(".export_menu");
    var overlay = this._draw_target_overlay();

    if(export_menu.length){
      export_menu.fadeIn("fast");
    }
    else{
      export_menu = $("<div class='export_menu'><div><button type='submit' data-action='open'>Show image</button></div><div><button type='submit' data-action='save'>Save image</button></div></div>");
      this.target_element.append(export_menu);
      export_menu.css({"position": "absolute",
                      "top": 45,
                      "left": self.settings.width - 125,
                      "font-size": "12px",
                      "border": "solid #D2D2D2 1px",
                      "border-radius": "5px",
                      "padding": "2px",
                      "background-color": "white"});

      var buttons = export_menu.find("button");
      buttons.css({"padding-top": "7px", "padding-bottom": "5px", "padding-right": "8px", "padding-left": "8px", "color": "white", "border": "solid #D2D2D2 1px", "width": "100%", "background-color": "#2171b5", "font-weight": "bold"});  

      buttons.hover(
        function(){$(this).css({"cursor": "pointer", "opacity": 0.7})},
        function(){$(this).css({"opacity": 1})}
      );

      overlay.click(function(){
        export_menu.fadeOut("fast");
        overlay.fadeOut("fast");
      });

      buttons.click(function(){
        var action = $(this).attr("data-action");
        var zoom = 3;
        var width = self.stage.width();
        var height = self.stage.height();
        var loading_div = $("<h3 style='margin-top: 100px; margin-left: 100px; width: " + width + "px; height: " + height + "px;'>Loading...</h3>");
        self.target_element.after(loading_div);
        self.target_element.hide();
        self.stage.width(width*zoom);
        self.stage.height(height*zoom);
        self.stage.scale({x: zoom, y:zoom});
        self.stage.draw();
        self.navigation_layer.hide();
        self.stage.toDataURL({
          quality: 1,
          callback: function(dataUrl){
            if(action === "open"){
              open_image(dataUrl);
            }
            else{
              download_image(dataUrl);
            }
            self.stage.width(width);
            self.stage.height(height);
            self.stage.scale({x: 1, y:1});
            self.stage.draw();
            loading_div.remove();
            self.target_element.show();
            self.navigation_layer.show();
            self.navigation_layer.draw();
            overlay.trigger("click");
          }
        });
      });
    }

    function download_image(dataUrl){
      $('<a download="inchlib" href="'+ dataUrl + '"></a>')[0].click();
    };
    
    function open_image(dataUrl){
      window.open(dataUrl, '_blank');
    };
  };

  InCHlib.prototype._color_scale_click = function(icon, evt){
    var self = this;
    var i, option, key, value;
    var color_options = {"heatmap_colors": "Heatmap data colors"};

    var value_options = {"max_quantile": "Max quantile value",
                        "min_quantile": "Min quantile value",
                        "middle_quantile": "Middle quantile value",
                      };

    if(this.settings.metadata){
      color_options["metadata_colors"] = "Metadata colors";
    }

    if(this.settings.column_metadata){
      color_options["column_metadata_colors"] = "Column metadata colors";
    }

    var form_id = "settings_form_" + this.settings.target;
    var settings_form = $("#" + form_id);
    var overlay = this._draw_target_overlay();

    if(settings_form.length){
      settings_form.fadeIn("fast");
    }
    else{
      settings_form = $("<form class='settings_form' id='" + form_id + "'></form>");
      var options = "", color_1, color_2, color_3;

      for(i = 0, keys = Object.keys(color_options), len = keys.length; i < len; i++){
        key = keys[i];
        color_1 = this._get_color_for_value(0,0,1,0.5,this.settings[key]);
        color_2 = this._get_color_for_value(0.5,0,1,0.5,this.settings[key]);
        color_3 = this._get_color_for_value(1,0,1,0.5,this.settings[key]);

        option = "<div><div class='form_label'>" + color_options[key] + "</div><input type='text' name='" + key +"' value='"+ this.settings[key] + "'/> <div class='color_button' style='background: linear-gradient(to right, " + color_1 + "," + color_2 + "," + color_3 + ")'></div></div>";
        options += option;
      }

      for(i = 0, keys = Object.keys(value_options), len = keys.length; i < len; i++){
        key = keys[i];
        option = "<div><div class='form_label'>" + value_options[key] + "</div><input type='text' name='" + key +"' value='"+ this.settings[key] + "'/></div>";
        options += option;
      }
      option = "<div><div class='form_label'>Heatmap coloring</div>\
                <select name='independent_columns'>"
      
      if(this.settings.independent_columns){
        option += "<option value='true' selected>By columns</option>\
                  <option value='false'>Entire heatmap</option>"
      }
      else{
        option += "<option value='true'>By columns</option>\
                  <option value='false' selected>Entire heatmap</option>" 
      }
      option += "</select></div>";
      options += option;

      options = options + '<button type="submit">Redraw</button>'
      settings_form.html(options);

      this.target_element.append(settings_form);
      settings_form.css({"z-index": 1000, "position": "absolute", "top": 110, "left": 0, "padding": "10px", "border": "solid #D2D2D2 2px", "border-radius": "5px", "background-color": "white"});
      $("#" + form_id + " .color_button").css({"border": "solid #D2D2D2 1px", "height": "15px", "width": "30px", "display": "inline-block"});  
      $("#" + form_id + " > div").css({"font-size": "12px", "margin-bottom": "10px"});  
      $("#" + form_id + " input").css({"border-radius": "5px", "width": "100px"});  
      $("#" + form_id + " .form_label").css({"color": "gray", "margin-bottom": "5px", "font-style": "italic"});  
      $("#" + form_id + " button").css({"padding-top": "7px", "padding-bottom": "5px", "padding-right": "5px", "padding-left": "5px", "color": "white", "border": "solid #D2D2D2 1px", "border-radius": "5px", "width": "100%", "background-color": "#2171b5", "font-weight": "bold"});  
    }

    overlay.click(function(){
      settings_form.fadeOut("fast");
      overlay.fadeOut("fast");
    });

    var color_buttons = $("#" + form_id + " .color_button");
    
    color_buttons.hover(
      function(){$(this).css({"cursor": "pointer", "opacity": 0.7})},
      function(){$(this).css({"opacity": 1})}
    );

    color_buttons.click(function(evt){
      self._draw_color_scales_select(this, evt);
    });

    settings_form.submit(function(evt){
      var settings = {};
      var settings_fieldset = $(this).find("input, select");

      settings_fieldset.each(function(){
          option = $(this);
          key = option.attr("name");
          value = option.val();
          if(value != ""){
              if(value === "true"){
                  value = true;
              }
              else if(value === "false"){
                value = false;
              }
              settings[key] = value;
          }
      });
      self.update_settings(settings);
      self.redraw_heatmap();
      self._update_color_scale();
      overlay.trigger('click');
      evt.preventDefault();
      evt.stopPropagation();
    })
  }

  InCHlib.prototype._draw_color_scales_select = function(element, evt){
    var scales_div = this.target_element.find(".color_scales");
    var scale_divs;
    var self = this;

    if(scales_div.length){
      scales_div.fadeIn("fast");
      scale_divs = scales_div.find(".color_scale");
    }
    else{
      scales_div = $("<div class='color_scales'></div>");
      var scale, color_1, color_2, color_3, key;

      for(var i = 0, keys = Object.keys(this.colors), len = keys.length; i < len; i++){
        key = keys[i];
        color_1 = this._get_color_for_value(0,0,1,0.5,key);
        color_2 = this._get_color_for_value(0.5,0,1,0.5,key);
        color_3 = this._get_color_for_value(1,0,1,0.5,key);
        scale = "<div class='color_scale' data-scale_acronym='" + key + "' style='background: linear-gradient(to right, " + color_1 + "," + color_2 + "," + color_3 + ")'></div>";
        scales_div.append(scale);
      }
      this.target_element.append(scales_div);
      scales_div.css({"border": "solid #D2D2D2 2px",
                     "border-radius": "5px",
                     "padding": "5px",
                     "position": "absolute",
                     "top": 110,
                     "left": 170,
                     "background-color": "white"});

      scale_divs = this.target_element.find(".color_scale");
      scale_divs.css({"margin-top":"3px",
                      "width": "80px",
                      "height": "20px",
                      "border": "solid #D2D2D2 1px",});

      scale_divs.hover(
        function(){$(this).css({"cursor": "pointer", "opacity": 0.7})},
        function(){$(this).css({"opacity": 1})}
      );

      this.target_element.find(".target_overlay").click(function(){
        scales_div.fadeOut("fast");
      });
    }

    scale_divs.on("click", function(){
      var color = $(this).attr("data-scale_acronym");
      var input = $(element).prev("input:first").val(color);
      $(element).css({"background": "linear-gradient(to right, " + self._get_color_for_value(0,0,1,0.5,color) + "," + self._get_color_for_value(0.5,0,1,0.5,color) + "," + self._get_color_for_value(1,0,1,0.5,color) + ")"})
      scales_div.fadeOut("fast");
      scale_divs.off("click");
    });

  };

  InCHlib.prototype._color_scale_mouseover = function(color_scale, layer){
      var label = color_scale.getAttr("label");
      var x = color_scale.getAttr("x");
      var y = color_scale.getAttr("y");

      this.icon_tooltip = this.objects_ref.tooltip_label.clone({x: x,
          y: y + 25
      });

      this.icon_tooltip.add(this.objects_ref.tooltip_tag.clone());
      this.icon_tooltip.add(this.objects_ref.tooltip_text.clone({text: label}));

      layer.add(this.icon_tooltip);
      this.icon_tooltip.moveToTop();
      color_scale.setOpacity(0.7);
      layer.draw();
  }

  InCHlib.prototype._color_scale_mouseout = function(color_scale, layer){
      this.icon_tooltip.destroy();
      color_scale.setOpacity(1);
      layer.draw();
  }

  InCHlib.prototype._unzoom_icon_click = function(){
      var current_node_id = this.zoomed_clusters["row"].pop();
      var node_id = (this.zoomed_clusters["row"].length === 0)?this.root_id:this.zoomed_clusters["row"].pop();
      this._zoom_cluster(node_id);
      this._highlight_cluster(current_node_id);
  };

  InCHlib.prototype._column_unzoom_icon_click = function(){
      var current_node_id = this.zoomed_clusters["column"].pop();
      var node_id = (this.zoomed_clusters["column"].length === 0)?this.column_root_id:this.zoomed_clusters["column"].pop();
      this._zoom_column_cluster(node_id);
      this._highlight_column_cluster(current_node_id);
  };

  InCHlib.prototype._icon_mouseover = function(icon, icon_overlay, layer){
      var label = icon.getAttr("label");
      var x = icon_overlay.getAttr("x");
      var y = icon_overlay.getAttr("y");
      var width = icon_overlay.getWidth();
      var height = icon_overlay.getHeight();

      if(icon.getAttr("id") === "export_icon"){
        x = x - 50;
      }

      this.icon_tooltip = this.objects_ref.tooltip_label.clone({x: x,
          y: y+1.2*height
      });

      this.icon_tooltip.add(this.objects_ref.tooltip_tag.clone());
      this.icon_tooltip.add(this.objects_ref.tooltip_text.clone({text: label}));

      layer.add(this.icon_tooltip);    
      icon.setFill("black");
      layer.draw();
  }

  InCHlib.prototype._icon_mouseout = function(icon, icon_overlay, layer){
      this.icon_tooltip.destroy();
      icon.setFill("grey");
      layer.draw();
  }

  InCHlib.prototype._dendrogram_layers_click=function(layer, evt){
      var path_id = evt.target.attrs.path_id;
      layer.fire("mouseout", layer, evt);
      
      if(path_id != this.last_highlighted_cluster){
          this._highlight_cluster(path_id);
      }
      else{
          this.unhighlight_cluster();
      }

      layer.draw();
      this.events.dendrogram_node_onclick(this.current_object_ids, path_id, evt);
  }

  InCHlib.prototype._column_dendrogram_layers_click=function(layer, evt){
      var object_ids = [];
      var path_id = evt.target.attrs.path_id;
      layer.fire("mouseout", layer, evt);
      
      if(path_id != this.last_highlighted_column_cluster){
          object_ids = this._highlight_column_cluster(path_id);
      }
      else{
          this.unhighlight_column_cluster();
      }

      layer.draw();
      this.events.column_dendrogram_node_onclick(this.current_column_ids, path_id, evt);
  }

  InCHlib.prototype._dendrogram_layers_mousedown = function(layer, evt){
    var self = this;
    clearTimeout(this.timer);
    this.timer = setTimeout(function() {
        self._zoom_cluster(evt.target.attrs.path_id);
    }, 700);
  }

  InCHlib.prototype._column_dendrogram_layers_mousedown = function(layer, evt){
    var self = this;
    clearTimeout(this.timer);
    this.timer = setTimeout(function() {
        self._zoom_column_cluster(evt.target.attrs.path_id);
    }, 700);
  }

  InCHlib.prototype._dendrogram_layers_mouseup = function(layer, evt){
    clearTimeout(this.timer);
  }

  InCHlib.prototype._dendrogram_layers_mouseout = function(layer, evt){
    this.path_overlay.destroy();
    this.dendrogram_hover_layer.draw();
  }

  InCHlib.prototype._dendrogram_layers_mouseover = function(layer, evt){
    this.path_overlay = evt.target.attrs.path.clone({"strokeWidth": 4});
    
    this.dendrogram_hover_layer.add(this.path_overlay);
    this.dendrogram_hover_layer.draw();
  }

  InCHlib.prototype._visible_features_equal_column_dendrogram_count = function(){
      if((this.on_features["data"].length + this.on_features["metadata"].length) == this.current_column_count){
          return true;
      }
      return false;
  }

  InCHlib.prototype._get_color_for_value = function(value, min, max, middle, color_scale){
      var color = this.colors[color_scale];
      var c1 = color["start"];
      var c2 = color["end"];

      if(value > max){
        return 'rgb('+c2.r+','+c2.g+','+c2.b+')';
      }

      if(min == max || value < min){
        return 'rgb('+c1.r+','+c1.g+','+c1.b+')';
      }

      if(color["middle"] !== undefined){
          if(value >= middle){
              min = middle;
              c1 = color["middle"];
              c2 = color["end"];
          }
          else{
              max = middle;
              c1 = color["start"];
              c2 = color["middle"];
          }
      }

      var position = (value-min)/(max-min);
      var r = this._hack_round(c1.r+(position*(c2.r-c1.r)));
      var g = this._hack_round(c1.g+(position*(c2.g-c1.g)));
      var b = this._hack_round(c1.b+(position*(c2.b-c1.b)));
      return 'rgb('+r+','+g+','+b+')';
  }

  InCHlib.prototype._get_font_size = function(text_length, width, height, max_font_size){
      var max_possible_size = height - 2;
      var font_size = max_possible_size;

      if(font_size/2*text_length > width-10){
          font_size = font_size/(font_size/2*text_length/(width-10));
      };
      font_size = (font_size > max_possible_size)?max_possible_size:font_size;
      font_size = (font_size > max_font_size)?max_font_size:font_size;
      return font_size;
  }

  InCHlib.prototype._get_object_ids = function(node_id){
      if(this.data.nodes[node_id]["left_child"] !== undefined){
        this._get_object_ids(this.data.nodes[node_id]["left_child"]);
        this._get_object_ids(this.data.nodes[node_id]["right_child"]);
      }
      else{
        this.current_object_ids.push.apply(this.current_object_ids, this.data.nodes[node_id]["objects"])
      }
  }

  InCHlib.prototype._get_column_ids = function(node_id){
      this.current_column_ids = [];
      this._collect_column_ids(node_id);
  }

  InCHlib.prototype._collect_column_ids = function(node_id){
    if(this.column_dendrogram.nodes[node_id]["left_child"] !== undefined){
      this._collect_column_ids(this.column_dendrogram.nodes[node_id]["left_child"]);
      this._collect_column_ids(this.column_dendrogram.nodes[node_id]["right_child"]);
    }
    else{
      this.current_column_ids.push(this.nodes2columns[node_id]);
    }
  }

  InCHlib.prototype._hack_size = function(obj) {
      return Object.keys(obj).length;
  };

  InCHlib.prototype._hack_round = function(value){
      return (0.5 + value) >> 0;
  }

  InCHlib.prototype._middle2fnc = function(values){
      var self = this;
      var fncs = {"zero": function(values){return 0;},
        "median": function(values){
          values.sort(function(a,b){return a - b});
          var median_pos = self._hack_round(values.length/2);
          return values[median_pos];
          }, 
        "mean": function(values){
          var sum = values.reduce(function(a, b){return parseFloat(a) + b;});
          return sum/values.length;}
      };
      return fncs[this.settings.values_center](values);
  };

  InCHlib.prototype._is_number = function(n){
      return !isNaN(parseFloat(n)) && isFinite(n);
  }

  InCHlib.prototype._row_mouseenter = function(evt){
      var row_id = evt.target.parent.getAttr("id");
      var visible = this.on_features["data"].length + this.on_features["metadata"].length;
      
      if(this.settings.count_column){
          visible++;
      }

      if(evt.target.parent.attrs.class !== "column_metadata"){
          this.highlighted_row = row_id;
          var y = this.leaves_y_coordinates[row_id];
          var x = this.heatmap_distance;
          
          this.row_overlay = this.objects_ref.heatmap_line.clone({points: [x, y, x + this.heatmap_width, y],
            strokeWidth: this.pixels_for_leaf,
            stroke: "#FFFFFF",
            opacity: 0.3,
            listening: false});

          this.heatmap_overlay.add(this.row_overlay);
          this.heatmap_overlay.draw();
          this.events.row_onmouseover(this.data.nodes[row_id].objects, evt);
      }
  }

  InCHlib.prototype._row_mouseleave = function(evt){
      this.row_overlay.destroy();
      this.events.row_onmouseout(evt);
  };

  InCHlib.prototype._draw_col_label = function(evt){
      var i, line;
      var attrs = evt.target.attrs;
      var points = attrs.points;
      var x = this._hack_round((points[0] + points[2])/2);
      var y = points[1]-0.5*this.pixels_for_leaf;
      var column = attrs.column.split("_");
      var header_type2value = {"d": this.heatmap_header[column[1]],
                               "m": this.metadata_header[column[1]],
                               "Count": "Count"};
      
      var value = attrs.value;
      var header = header_type2value[column[0]];

      if(header !== this.last_column){
        this.column_overlay.destroy();
        this.last_column = attrs.column;
        this.column_overlay = this.objects_ref.heatmap_line.clone({points: [x, this.header_height, x, this.header_height + this.column_metadata_height + (this.heatmap_array.length+0.5)*this.pixels_for_leaf],
          strokeWidth: this.pixels_for_dimension,
          stroke: "#FFFFFF",
          opacity: 0.3,
          listening: false});

        this.heatmap_overlay.add(this.column_overlay);
      }
      
      if(header !== undefined){
          value = [header, value].join("\n");
      }

      var tooltip = this.objects_ref.tooltip_label.clone({x: x, y:y, id: "col_label",});
      tooltip.add(this.objects_ref.tooltip_tag.clone({pointerDirection: 'down'}), this.objects_ref.tooltip_text.clone({text: value}));
      
      this.heatmap_overlay.add(tooltip);
      this.heatmap_overlay.moveToTop();
      this.heatmap_overlay.draw();
      return;
  }

  InCHlib.prototype._unprefix = function(prefixed){
      return prefixed.split(this.settings.target+"#")[1];
  }

  InCHlib.prototype._prefix = function(nonprefixed){
      return this.settings.target + "#" + nonprefixed;
  }

  /**
    * Returns array of features for leaf by its ID. When sent leaf ID is not present, false is returned
    */
  InCHlib.prototype.get_features_for_leaf = function(leaf_id){
      if(this.objects2leaves[leaf_id] !== undefined){
        var row_id = this.objects2leaves[leaf_id];
        return this.data.nodes[row_id].features;
      }
      return false;
  }

  /**
    * Adds a user defined color scale defined by its name start color, end color and optionaly middle color
    */

  InCHlib.prototype.add_color_scale = function(color_scale_name, color_scale){
      this.colors[color_scale_name] = color_scale;
  }

  InCHlib.prototype._get_visible_count = function(){
      var visible = this.on_features["data"].length + this.on_features["metadata"].length;
      if(this.settings.count_column){
          visible++;
      }
      return visible;
  }

  /**
    * Update cluster heatmap settings
    */
  InCHlib.prototype.update_settings = function(settings_object){
    $.extend(this.settings, settings_object);
    return;
  }

  /**
    * Redraw cluster heatmap
    */
  InCHlib.prototype.redraw = function(){
    this._delete_all_layers();
    this.draw();
    return;
  }

  InCHlib.prototype.redraw_heatmap = function(){
    this._delete_layers([this.heatmap_layer]);
    this._draw_heatmap();
    this.heatmap_layer.moveToBottom();
    this.heatmap_layer.moveUp();
    return;
  }
}(jQuery));