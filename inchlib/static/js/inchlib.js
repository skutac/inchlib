function InCHlib(settings){
    var target_width = $("#" + settings.target).width();
    // When measuring the rendering duration
    // this.start = new Date().getTime();

    this.settings = {
        "target" : settings.target,
        "heatmap" : true,
        "dendrogram": true,
        "metadata": false,
        "max_height" : 800,
        "width" : target_width,
        "heatmap_colors" : "Greens",
        "heatmap_font_color" : "black",
        "heatmap_part_width" : 0.7,
        "column_dendrogram" : false,
        "independent_columns" : true,
        "metadata_colors" : "Oranges",
        "highlight_colors" : "Reds",
        "highlighted_rows" : [],
        "label_color": "#9E9E9E",
        "count_column": false,
        "count_column_colors": "Reds",
        "min_row_height": false,
        "max_row_height": 25,
        "max_column_width": 100,
        "font": "Arial",
        "values_center": "median",
        "draw_row_ids": false,
        "header_as_heatmap_row": false,
        "header_row_colors": "YlOrB",
        "current_row_ids_callback": function(row_ids){
            return;
        },
    };

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
            "YlOrR": {"start": {"r":255, "g": 255, "b": 178}, "end": {"r": 227, "g": 26, "b": 28}, "middle": {"r": 204, "g": 76, "b": 2}},
            "YlOrB": {"start": {"r":255, "g": 255, "b": 212}, "end": {"r": 5, "g": 112, "b": 176}, "middle": {"r": 204, "g": 76, "b": 2}},
            "PRGn2": {"start": {"r":123, "g": 50, "b": 148}, "end": {"r": 0, "g": 136, "b": 55}, "middle": {"r":202, "g": 0, "b": 32}},
            "PiYG2": {"start": {"r":208, "g": 28, "b": 139}, "end": {"r": 77, "g": 172, "b": 38}, "middle": {"r":255, "g": 255, "b": 178},},
            "YlGnBu": {"start": {"r":255, "g": 255, "b": 204}, "end": {"r": 34, "g": 94, "b": 168}, "middle": {"r": 35, "g": 132, "b": 67}},
            "RdYlBu": {"start": {"r":215, "g": 25, "b": 28}, "end": {"r": 44, "g": 123, "b": 182}, "middle": {"r":255, "g": 255, "b": 178}},
            "RdYlGn": {"start": {"r":215, "g": 25, "b": 28}, "end": {"r": 26, "g": 150, "b": 65}, "middle": {"r":255, "g": 255, "b": 178}},
            "BuWhRd": {"start": {"r": 33, "g": 113, "b": 181}, "middle": {"r": 255, "g": 255, "b": 255}, "end": {"r":215, "g": 25, "b": 28}},
            "RdBkGr": {"start": {"r":215, "g": 25, "b": 28}, "middle": {"r": 0, "g": 0, "b": 0}, "end": {"r": 35, "g": 139, "b": 69}},
            "RdGyBu": {"start": {"r":215, "g": 25, "b": 28}, "middle": {"r":254, "g": 229, "b": 217}, "end": {"r": 44, "g": 123, "b": 182}},
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
    this.rows = {};
}

InCHlib.prototype.read_data = function(json){
    this.data = json;
    this.data.header = this.data.data.header;
    this.data.nodes = this.data.data.nodes;
}

InCHlib.prototype.read_data_from_file = function(json){
    var self = this;

    $.ajax({
        type: 'GET',
        url: json,
        dataType: 'json',
        success: function(data){
            self.data = data;
            self.data.header = self.data.data.feature_names;
            self.data.nodes = self.data.data.nodes;
        },
        async: false
    });
}

InCHlib.prototype._add_prefix = function(){
    var i, id;
    for(i in this.data.nodes){
        id = [this.settings.target, i].join("#");
        this.data.nodes[id] = this.data.nodes[i];
        delete this.data.nodes[i];
        if("parent" in this.data.nodes[id]){
            this.data.nodes[id].parent = [this.settings.target, this.data.nodes[id].parent].join("#");
        }
        if(this.data.nodes[id]["count"] != 1){
            this.data.nodes[id].left_child = [this.settings.target, this.data.nodes[id].left_child].join("#");
            this.data.nodes[id].right_child = [this.settings.target, this.data.nodes[id].right_child].join("#");
        }
    }

    if(this.settings.metadata){
        for(i in this.data.metadata.nodes){
            id = [this.settings.target, i].join("#");
            this.data.metadata.nodes[id] = this.data.metadata.nodes[i];
            delete this.data.metadata.nodes[i];
        }
    }
    return;
}

InCHlib.prototype._get_root_id = function(nodes){
    var i, root_id;
    for(i in nodes){
        if(!("parent" in nodes[i])){
            root_id = i;
            break;
        }
    }
    return root_id;
}

InCHlib.prototype._get_dimensions = function(nodes){
    var dimensions, i;

    for(i in nodes){
        if(nodes[i].count == 1){
            dimensions = nodes[i].features.length;
            break;
        }
        else if(Object.prototype.toString.call(nodes[i]) == "[object Array]"){
            dimensions = nodes[i].length;
            break;
        };
    }
    return dimensions;
}

InCHlib.prototype._get_min_max_middle = function(data){
    var i;
    var min_max_middle = [];
    var all = [];

    for(i = 0; i<data.length; i++){
        all = all.concat(data[i]);
    }

    min_max_middle.push(Math.min.apply(null, all));
    min_max_middle.push(Math.max.apply(null, all));
    min_max_middle.push(this._middle2fnc(all));
    return min_max_middle;

}

InCHlib.prototype._get_columns_min_max_middle = function(data){
    var self = this;
    var columns = [];
    var i, j, value;
    var data_length = data[0].length;

    for(i = 0; i<data_length; i++){
        columns.push([]);
    }

    for(i = 0; i<data.length; i++){
        for(j = 0; j < data_length; j++){
            value = data[i][j];
            if(typeof value != "undefined"){
                columns[j].push(value);
            }
        }
    }

    var columns_min_max = [], min, max, middle;

    for(i = 0; i<columns.length; i++){
        if(this._is_number(columns[i][0])){
            columns[i] = columns[i].map(parseFloat)
            min = Math.min.apply(null, columns[i]);
            max = Math.max.apply(null, columns[i]);
            middle = this._middle2fnc(columns[i]);
        }
        else{
            var hash_object = this._get_hash_object(columns[i]);
            min = 0;
            max = this._hack_size(hash_object)-1;
            middle = max/2;
            this.categories2numbers[i] = hash_object;
        }
        columns_min_max.push([min,max,middle]);
    }
    return columns_min_max;
}

InCHlib.prototype._get_hash_object = function(array){
    var i, count=0, hash_object = {};

    for(i = 0; i<array.length; i++){
        if(typeof(hash_object[array[i]]) == "undefined"){
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
    var i, j, node_data;

    for(i in nodes){
        if(nodes[i].count == 1){
            node_data = nodes[i].features;
            for(j = 0; j < node_data.length; j++){
                if((""+node_data[j]).length > max_length){
                    max_length = (""+node_data[j]).length;
                }
            }
        }
    }
    
    if(this.settings.metadata){
        var nodes = this.data.metadata.nodes;
        for(i in nodes){
            node_data = nodes[i];
            for(j = 0; j < node_data.length; j++){
                if((""+node_data[j]).length > max_length){
                    max_length = (""+node_data[j]).length;
                }
            }
        }
    }

    return max_length;
}

InCHlib.prototype._preprocess_heatmap_data = function(){
    var heatmap_array = [], i, j = 0, data, node;

    for(i in this.data.nodes){
        node = this.data.nodes[i];
        if(node.count == 1){
            data = node.features;
            heatmap_array.push([i]);
            heatmap_array[j].push.apply(heatmap_array[j], data);
            if(this.settings.metadata){
                heatmap_array[j].push.apply(heatmap_array[j], this.data.metadata.nodes[i]);
            }
            j++;
        }
    }
    return heatmap_array;
}

InCHlib.prototype._reorder_heatmap = function(column_index){
    var i;
    this.leaves_y_coordinates = {};
    
    if(this.ordered_by_index == column_index){
        this.heatmap_array.reverse();
    }
    else{
        if(typeof(this.heatmap_array[0][column_index]) == "string"){
            this.heatmap_array.sort(function(a,b){return (a[column_index] == null)?-1:(b[column_index] == null)?1:(a[column_index] > b[column_index])?1:(a[column_index] < b[column_index])?-1:0});
        }
        else{
            this.heatmap_array.sort(function(a,b){return (a[column_index] == null)?-1:(b[column_index] == null)?1:a[column_index] - b[column_index]});
        }
    }

    var pixels_for_leaf = this.pixels_for_leaf;
    var y = pixels_for_leaf/2 + this.header_height;

    for(i = 0; i<this.heatmap_array.length; i++){
        this.leaves_y_coordinates[this.heatmap_array[i][0]] = y;
        y += pixels_for_leaf;
    }

    this.ordered_by_index = column_index;
    return;
}

/* *************************************************
Draw functions
************************************************* */

InCHlib.prototype.redraw = function(){
    return false;
}

InCHlib.prototype.draw = function(){
    this._add_prefix();
    this.dimensions = this._get_dimensions(this.data.nodes);
    this.zoomed_clusters = [];
    this.heatmap_array = this._preprocess_heatmap_data();
    this.data_dimensions = this.dimensions;

    if(this.settings.metadata){
        this.metadata_dimensions = this._get_dimensions(this.data.metadata.nodes);
        this.dimensions = this.data_dimensions+this.metadata_dimensions;
    }
    else{
        this.metadata_dimensions = 0
    }

    if(this.settings.heatmap){
        this._set_heatmap_settings();
    }
    else{
        this.dimensions = 0;
        this.visible_features = 0;
    }

    if(this.settings.column_dendrogram && this.heatmap_header){
        this.footer_height = 150;
    }

    this.stage = new Kinetic.Stage({
        container: this.settings.target,
    });

    this.pixels_for_leaf = (this.settings.max_height-this.header_height-this.footer_height)/this.heatmap_array.length;
    if(this.pixels_for_leaf > this.settings.max_row_height){
        this.pixels_for_leaf = this.settings.max_row_height;
    }

    if(this.settings.min_row_height > this.pixels_for_leaf){
        this.pixels_for_leaf = this.settings.min_row_height;
    }

    this.settings.height = this.heatmap_array.length*this.pixels_for_leaf+this.header_height+this.footer_height;
    this.stage.setWidth(this.settings.width);
    this.stage.setHeight(this.settings.height);
    this._draw_stage_layer();
    this._draw_navigation();

    if(this.settings.dendrogram){
        this.root_id = this._get_root_id(this.data.nodes);
        this.distance = this.settings.width - this.heatmap_width-100;
        this._draw_row_dendrogram(this.root_id);
        this._draw_heatmap_header();
        this._draw_heatmap();

        if(this.settings.column_dendrogram){
            this.column_root_id = this._get_root_id(this.data.column_dendrogram.nodes);
            this._draw_column_dendrogram(this.column_root_id);
        }
    }
    else{
        this.settings.heatmap_part_width = 1;
        this.settings.column_dendrogram = false;
        this.distance = this._hack_round((this.settings.width-this.heatmap_width)/2);
        this._reorder_heatmap(0);
        this.ordered_by_index = 0;
        this._draw_heatmap_header();
        this._draw_heatmap();
    }

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
    this.nodes_y_coordinates = {};
    this.leaves_y_coordinates = {};

    this.pixels_for_leaf = (this.settings.max_height-this.header_height-this.footer_height)/count;
    if(this.pixels_for_leaf > this.settings.max_row_height){
        this.pixels_for_leaf = this.settings.max_row_height;
    }

    if(this.settings.min_row_height > this.pixels_for_leaf){
        this.pixels_for_leaf = this.settings.min_row_height;
    }

    this.settings.height = count*this.pixels_for_leaf+this.header_height+this.footer_height;
    this.stage.setWidth(this.settings.width);
    this.stage.setHeight(this.settings.height);

    var current_left_count = 0;
    var current_right_count = 0;
    var y = this.header_height + this.pixels_for_leaf/2;
    if(node.count > 1){
        current_left_count = this.data.nodes[node.left_child].count;
        current_right_count = this.data.nodes[node.right_child].count;
    }
    this._draw_row_dendrogram_node(node_id, node, current_left_count, current_right_count, 0, y);
    this.middle_item_count = (this.min_item_count+this.max_item_count)/2;
    this._draw_distance_scale(node.distance);
    this.stage.add(this.dendrogram_layer);

    this.dendrogram_overlay_layer = new Kinetic.Layer();
    this.stage.add(this.dendrogram_overlay_layer);

    this.dendrogram_hover_layer = new Kinetic.Layer();
    this.stage.add(this.dendrogram_hover_layer);

    var self = this;

    this.dendrogram_layer.on("mouseover", function(evt){
        self._dendrogram_layers_mouseover(evt);
    });

    this.dendrogram_overlay_layer.on("mouseover", function(evt){
        self._dendrogram_layers_mouseover(evt);
    });

    this.dendrogram_layer.on("mouseout", function(evt){
        self._dendrogram_layers_mouseout(evt);
    });    

    this.dendrogram_overlay_layer.on("mouseout", function(evt){
        self._dendrogram_layers_mouseout(evt);
    });
    
    this.dendrogram_layer.on("click", function(evt){
        self._dendrogram_layers_click(evt, this);
    });

    this.dendrogram_overlay_layer.on("click", function(evt){
        self._dendrogram_layers_click(evt, this);
    });
    
    this.dendrogram_layer.on("dblclick", function(evt){
        if(evt.targetNode.attrs.path_id!=self.root_id){
            self._zoom_cluster(evt.targetNode.attrs.path_id);
        }
    });
}

InCHlib.prototype._draw_row_dendrogram_node = function(node_id, node, current_left_count, current_right_count, x, y){
    if(node.count != 1){
        var nodes = this.data.nodes;

        this.nodes_y_coordinates[node_id] = y;
        var node_neighbourhood = this._get_node_neighbourhood(node, nodes);
        var right_child = nodes[node.right_child];
        var left_child = nodes[node.left_child];

        var y1 = this._get_y1(node_neighbourhood, current_left_count, current_right_count);
        var y2 = this._get_y2(node_neighbourhood, current_left_count, current_right_count);
        var x1 = this._hack_round(this.distance - this.distance_step*node.distance);
        x1 = (x1 == 0)? 2: x1;
        
        var x2 = x1;
        var left_distance = this.distance - this.distance_step*nodes[node.left_child].distance;
        var right_distance = this.distance - this.distance_step*nodes[node.right_child].distance;

        if(right_child.count == 1){
            y2 = y2 + this.pixels_for_leaf/2;
        }

        this.dendrogram_layer.add(this._draw_horizontal_path(node_id, x1, y1, x2, y2, left_distance, right_distance));
        this._draw_row_dendrogram_node(node.left_child, left_child, current_left_count - node_neighbourhood.left_node.right_count, current_right_count + node_neighbourhood.left_node.right_count, left_distance, y1);
        this._draw_row_dendrogram_node(node.right_child, right_child, current_left_count + node_neighbourhood.right_node.left_count, current_right_count - node_neighbourhood.right_node.left_count, right_distance, y2);
    }
    else{
        this.leaves_y_coordinates[node_id] = y;

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

    this.stage_layer.on("click", function(){
        self._neutralize_cluster();
    });

    this.stage_layer.on("mouseenter", function(evt){
        self._row_mouseout();
    });
}

InCHlib.prototype._draw_column_dendrogram = function(node_id){
    this.column_dendrogram_layer = new Kinetic.Layer();
    var node = this.data.column_dendrogram.nodes[node_id];
    this.vertical_distance = (this.settings.header_as_heatmap_row)?this.header_height-18:this.header_height-5;
    this.vertical_distance_step = this.vertical_distance/node.distance;

    var current_left_count = this.data.column_dendrogram.nodes[node.left_child].count;
    var current_right_count = this.data.column_dendrogram.nodes[node.right_child].count;
    this._draw_column_dendrogram_node(node_id, node, current_left_count, current_right_count, 0, 0);
    this.stage.add(this.column_dendrogram_layer);
}

InCHlib.prototype._delete_all_layers = function(){
    this.stage.destroyChildren();
}

InCHlib.prototype._set_heatmap_settings = function(){
    var i, j, current_array, node;

    this.header = [];
    for(i = 0; i<this.dimensions; i++){
        this.header.push("");
    }

    this.heatmap_header = false;
    this.metadata_header = false;
    this.current_label = null;
    this.categories2numbers = {};

    if(this.data.header && this.data.header.length > 0){
        this.heatmap_header = this.data.header;
        for(i=0; i<this.data_dimensions; i++){
            this.header[i] = this.heatmap_header[i];
        }

        if(this.heatmap_header.length > 0 && this.settings.header_as_heatmap_row){
            var header_column = [];
            for(i = 0; i<this.heatmap_header.length; i++){
                header_column.push([this.heatmap_header[i]]);
            }
            this.header_min_max = this._get_columns_min_max_middle(header_column)[0];
            if(this.categories2numbers[0]){
                this.categories2numbers["header"] = this.categories2numbers[0];
            };
            var empty_arr = [];
            for(var i = 0; i < this.heatmap_header.length; i++) {
                empty_arr.push("");
            }
            this.data.nodes["header_row"] = {"features": empty_arr};
        }
    }
    else{
        this.settings.header_as_heatmap_row = false;
    }

    var data = [];
    for(i in this.data.nodes){
        node = this.data.nodes[i];
        if(node.count == 1){
            data.push(node.features);
        };
    }
        
    if(this.settings.independent_columns){
        this.columns_min_max = this._get_columns_min_max_middle(data);
    }
    else{
        var min_max_middle = this._get_min_max_middle(data);
        this.max_value = min_max_middle[1];
        this.min_value = min_max_middle[0];
        this.middle_value = min_max_middle[2];
    }

    if(this.settings.metadata){

        if(this.data.metadata.feature_names){
            this.metadata_header = this.data.metadata.feature_names;

            for(i=0; i<this.metadata_dimensions; i++){
                this.header[this.data_dimensions+i] = this.metadata_header[i];
            }
        }

        var metadata = [];

        for(i in this.data.metadata.nodes){
            metadata.push(this.data.metadata.nodes[i]);
        }
        this.metadata_columns_min_max = this._get_columns_min_max_middle(metadata);
    }

    if(this.settings.count_column){
        this.max_item_count = 1;
        this.min_item_count = 1;
        this.dimensions++;
        this.header.push("Count");
    }

    if(!(this.heatmap_header || this.metadata_header) && !(this.settings.count_column)){
        this.header = false;
    }

    this.settings.heatmap_part_width = this.dimensions?this.settings.heatmap_part_width:0;
    this.heatmap_width = this.dimensions?this._hack_round(this.settings.heatmap_part_width*(this.settings.width-100)):0;
    this.pixels_for_dimension = this.dimensions?this.heatmap_width/this.dimensions:0;

    if(this.settings.max_column_width && this.settings.max_column_width < this.pixels_for_dimension){
        this.pixels_for_dimension = this.settings.max_column_width;
        this.heatmap_width = this.dimensions*this.pixels_for_dimension;
    }

    this.features = {};
    for(i=0; i<this.dimensions; i++){
        this.features[i] = 1;
    }

    this.visible_features = 0
    for(i in this.features){
        if(this.features[i] == 1){
            this.visible_features++;
        }
    }
}

InCHlib.prototype._set_on_features = function(){
    this.on_features = [];
    this.on_metadata_features = [];

    for(var i in this.features){
        if(this.features[i] == 1){
            if(i < this.data_dimensions){
                this.on_features.push(parseInt(i));
            }
            else{
                if(!this.settings.count_column || parseInt(i) < this.dimensions-1){
                    this.on_metadata_features.push(parseInt(i)-this.data_dimensions);
                }
            }
        }
    }
}

InCHlib.prototype._draw_heatmap = function(){
    if(!this.settings.heatmap || this.dimensions==0){
        return;
    }

    var leaf_id, heatmap_row, row_id, col_number, col_label, i;
    this.heatmap_layer = new Kinetic.Layer();
    this.heatmap_overlay = new Kinetic.Layer();
    this.current_draw_values = true;

    this.max_value_length = this._get_max_value_length();
    this.value_font_size = this._get_font_size(this.max_value_length, this.pixels_for_dimension, this.pixels_for_leaf, 12);

    if(this.value_font_size < 4){
        this.current_draw_values = false;
    }
    

    this.value_text_ref = new Kinetic.Text({
        fontSize: this.value_font_size,
        fontFamily: this.settings.font,
        fill: this.settings.heatmap_font_color,
        listening: false,
    });


    this.line_ref = new Kinetic.Line({
            strokeWidth: this.pixels_for_leaf,
            lineCap: 'butt',
    });

    this._set_on_features();

    var x1 = this._hack_round(this.distance+this.dendrogram_heatmap_distance);
    var current_leaves_y = [];

    for(leaf_id in this.leaves_y_coordinates){
        var y = this.leaves_y_coordinates[leaf_id];
        heatmap_row = this._draw_heatmap_row(leaf_id, x1, y);
        this.rows[leaf_id] = heatmap_row;
        this.heatmap_layer.add(heatmap_row);
        current_leaves_y.push([leaf_id, y]);
    }

    if(this.settings.header_as_heatmap_row){
        heatmap_row = this._draw_header_row(x1);
        this.rows["header_row"] = heatmap_row;
        this.heatmap_layer.add(heatmap_row);
    }

    if(this.settings.draw_row_ids){
        this._draw_row_ids(current_leaves_y);
    }

    this.stage.add(this.heatmap_layer)
    this.stage.add(this.heatmap_overlay);

    this.highlighted_rows_layer = new Kinetic.Layer();
    this.stage.add(this.highlighted_rows_layer);
    this.highlighted_rows_layer.moveToTop();
    this.last_highlighted_row = null;
    this.last_col_number = null;

    var self = this;

    this.heatmap_layer.on("click", function(evt){
        if(self.highlighted_row != "header_row"){
            var items = self.data.nodes[self.highlighted_row].objects;
            var item_ids = [];
            
            for(i = 0; i < items.length; i++){
                item_ids.push(items[i]);
            }

            self.settings.current_row_ids_callback(item_ids);

        }
    });

    this.heatmap_layer.on("mouseover", function(evt){
        var row_id = evt.targetNode.parent.getAttr("id");
        var target_class = evt.targetNode.className;
        if(target_class == "Text"){
            return;
        }

        if(target_class == "Line"){
            col_number = self._hack_round((self.stage.getPointerPosition().x-self.distance-self.dendrogram_heatmap_distance-0.5*self.pixels_for_dimension)/self.pixels_for_dimension);
            row_id = evt.targetNode.parent.getAttr("id");
            var row_values = [];
            for(i = 0; i < self.data.nodes[row_id].features.length; i++){
                row_values.push(self.data.nodes[row_id].features[i]);
            }

            if(self.settings.metadata && col_number >= row_values.length){
                row_values = row_values.concat(self.data.metadata.nodes[row_id]);
            }

            if(self.settings.count_column){
                row_values.push(self.data.nodes[row_id].objects.length);
            }
            
            if(row_id != self.last_highlighted_row){
                self._row_mouseout();
                if(row_id != "header_row"){
                    self._row_mouseover(evt.targetNode.parent, col_number);
                }
                self._draw_col_label(evt.targetNode.parent, row_values, col_number);
                self.last_highlighted_row = row_id;
                self.last_col_number = col_number;
            }
            else if(col_number != self.last_col_number){
                col_label = self.stage.get('#col_label');
                col_label.destroy();
                self._draw_col_label(evt.targetNode.parent, row_values, col_number);
                self.last_col_number = col_number;
            }
        }
    });

    this.heatmap_overlay.on("mouseenter", function(){
        self._row_mouseout();
    });

}

InCHlib.prototype._draw_heatmap_row = function(node_id, x1, y1){
    var node = this.data.nodes[node_id];
    var row = new Kinetic.Group({id:node_id});
    var x, y, x2, y2, color, line, i, min_value, max_value, value, text, text_value, width, col_index;
    
    for (i = 0; i < this.on_features.length; i++){
        col_index = this.on_features[i];
        value = node.features[col_index];

        if(this.settings.independent_columns){
            this.min_value = this.columns_min_max[col_index][0];
            this.max_value = this.columns_min_max[col_index][1];
            this.middle_value = this.columns_min_max[col_index][2];
        }

        color = this._get_color_for_value(value, this.min_value, this.max_value, this.middle_value, this.settings.heatmap_colors);
        x2 = x1 + this.pixels_for_dimension;
        y2 = y1;

        line = this.line_ref.clone({
            stroke: color,
            points: [x1, y1, x2, y2]
        });
        row.add(line);

        if(this.current_draw_values){
            text = this.value_text_ref.clone({
                x: this._hack_round((x1 + x2)/2-(""+value).length*(this.value_font_size/4)),
                y: this._hack_round(y1-this.value_font_size/2),
                fontSize: this.value_font_size,
                text: value,
            });
            row.add(text);
        }
        x1 = x2;
    }

    if(this.settings.metadata){
        var metadata = this.data.metadata.nodes[node_id];
        for (i = 0; i < this.on_metadata_features.length; i++){
            col_index = this.on_metadata_features[i];
            value = metadata[col_index];
            text_value = value;
            this.metadata_min_value = this.metadata_columns_min_max[col_index][0];
            this.metadata_max_value = this.metadata_columns_min_max[col_index][1];
            this.metadata_middle_value = this.metadata_columns_min_max[col_index][2];

            if(!this._is_number(text_value)){
                value = this.categories2numbers[col_index][value];
            }

            color = this._get_color_for_value(value, this.metadata_min_value, this.metadata_max_value, this.metadata_middle_value, this.settings.metadata_colors);
            x2 = x1 + this.pixels_for_dimension;
            y2 = y1;
                
            line = this.line_ref.clone({
                    stroke: color,
                    points: [x1, y1, x2, y2]
                });
            row.add(line);

            if(this.current_draw_values){
                text = this.value_text_ref.clone({
                    text: text_value,
                })
                width = text.getWidth();
                x = this._hack_round((x1+x2)/2-width/2);
                y = this._hack_round(y1-this.value_font_size/2);
                text.position({x:x, y:y});
                row.add(text);
            }
            x1 = x2;
        }
    }

    if(this.settings.count_column && this.features[this.dimensions-1] == 1){
        x2 = x1 + this.pixels_for_dimension;
        var count = node.objects.length;
        color = this._get_color_for_value(count, this.min_item_count, this.max_item_count, this.middle_item_count, this.settings.count_column_colors);
        
        line = this.line_ref.clone({
                stroke: color,
                points: [x1, y1, x2, y2]
        });
        row.add(line);

        if(this.current_draw_values){
            text = this.value_text_ref.clone({
                text: count,
            })
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
    var font_size = this._get_font_size(max_length, 85, this.pixels_for_leaf, 12);
    var x = this.settings.width - 85;
    
    if(font_size > 4){
        for(i = 0; i < object_y.length; i++){
            text = this.value_text_ref.clone({
                x: x,
                y: this._hack_round(object_y[i][1] - font_size/2),
                fontSize: font_size,
                text: object_y[i][0],
            });
            this.heatmap_layer.add(text);
        }
    }
    return;
}

InCHlib.prototype._draw_header_row = function(x1){
    var row = new Kinetic.Group({id:"header_row"});
    var row_height = 15;
    var y1 = this.header_height - 0.5*row_height;
    var x, y, x2, y2, color, line, i, value, text, text_value, width, col_index;
    // var max_length = this._get_max_length(this.heatmap_header);
    // var font_size = 11;
    // var draw_values = true;
    
    // if(font_size/2*max_length > this.pixels_for_dimension-10){
    //     font_size = font_size/(font_size/2*max_length/(this.pixels_for_dimension-10));
    // };
    // font_size = (font_size > 12)?12:font_size;
    
    // if(font_size < 4){
    //     var draw_values = false;
    // }

    
    for (i = 0; i < this.on_features.length; i++){
        col_index = this.on_features[i];
        value = this.heatmap_header[col_index];
        text_value = value;

        if(!this._is_number(text_value)){
            value = this.categories2numbers["header"][value];
        }

        color = this._get_color_for_value(value, this.header_min_max[0], this.header_min_max[1], this.header_min_max[2], this.settings.header_row_colors);
        x2 = x1 + this.pixels_for_dimension;
        y2 = y1;
            
        line = this.line_ref.clone({
                strokeWidth: row_height,
                stroke: color,
                points: [x1, y1, x2, y2]
            });
        row.add(line);

        // if(draw_values){
        //     text = this.value_text_ref.clone({
        //         text: text_value,
        //         fontSize: font_size,

        //     })
        //     width = text.getWidth();
        //     x = this._hack_round((x1+x2)/2-width/2);
        //     y = this._hack_round(y1-font_size/2+2);
        //     text.position({x:x, y:y});
        //     row.add(text);
        // }
        x1 = x2;
    }
    return row;
}

InCHlib.prototype._draw_heatmap_header = function(){
    this.header_layer = new Kinetic.Layer();
    var row_count = this._hack_size(this.leaves_y_coordinates);
    var y = (this.settings.column_dendrogram && this.heatmap_header)? this.header_height+(this.pixels_for_leaf*row_count) + 10: this.header_height - 20;
    var rotation = (this.settings.column_dendrogram && this.heatmap_header) ? 45 : -45;
    y = (rotation == -45 && this.settings.header_as_heatmap_row)? y - 10: y;
    var distance = 0;
    var x, i, column_header;

    var items = [];
    for(i = 0; i<this.header.length; i++){
        items[i] = this.header[i];
    }
    var max_text_length = this._get_max_length(items);
    var font_size = this._hack_round(this.header_height*1.5/max_text_length);
    font_size = (font_size > 16)?16:font_size;
    font_size = (font_size > this.pixels_for_dimension)? this.pixels_for_dimension - 3: font_size;

    var column_header_ref = new Kinetic.Text({
                                fontSize: font_size,
                                fontFamily: this.settings.font,
                                fontStyle: "bold",
                                rotationDeg: rotation,
                                fill: 'black',
                            });

    if(font_size < 8){
        return;
    }
    
    for(i = 0; i<this.header.length; i++){
        if(this.features[i] == 1 || i>this.dimensions-1){
            x = this.distance+distance*this.pixels_for_dimension+this.pixels_for_dimension/2;
            column_header = column_header_ref.clone({
                    x: x,
                    y: y,
                    text: this.header[i],
                    position_index: i+1,
            });
            this.header_layer.add(column_header);
            distance++;
        }        
    }
    this.stage.add(this.header_layer);

    if(!(this.settings.dendrogram)){
        var self = this;
        this.header_layer.on("click", function(evt){
            var column = evt.targetNode;
            var order = column.attrs.order;
            var position_index = column.attrs.position_index;

            for(i = 0; i<self.header_layer.getChildren().length; i++){
                self.header_layer.getChildren()[i].setFill("black");
            }

            evt.targetNode.setAttrs({"order":order, "fill": "red"});

            self.heatmap_layer.destroyChildren();
            self.heatmap_layer.draw();
            self._reorder_heatmap(position_index, order);
            self._draw_heatmap();
        })

        this.header_layer.on("mouseover", function(evt){
            var label = evt.targetNode;
            label.setOpacity(0.7);
            this.draw();
        })

        this.header_layer.on("mouseout", function(evt){
            var label = evt.targetNode;
            label.setOpacity(1);
            this.draw();
        })
    }
}

InCHlib.prototype._draw_distance_scale = function(distance){
    this.distance_scale_layer = new Kinetic.Layer();
    var y1 = this.header_height-15;
    var y2 = y1;
    var x1 = 0;
    var x2 = this.distance;
    var path = new Kinetic.Line({
        points: [x1, y1, x2, y2],
        stroke: "black",
        listening: false,
    })

    this.distance_scale_layer.add(path);

    var circle = new Kinetic.Circle({
        x: x2, 
        y: y2,
        radius: 3,
        fill: "black",
        listening: false,
    })
    this.distance_scale_layer.add(circle);

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
            fontSize: 14,
            fontFamily: this.settings.font,
            fontStyle: 'bold',
            fill: 'black',
            align: 'right',
            listening: false,
    });
    this.distance_scale_layer.add(distance_number);

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
            this.distance_scale_layer.add(path);

            number = this._hack_round((number + marker_number_distance)*10)/10;
            if(number>10){
                number = this._hack_round(number);
            }
            
            marker_distance = marker_distance - marker_distance_step;
            marker_counter++;
        }
    }
    this.stage.add(this.distance_scale_layer);
}

InCHlib.prototype._draw_navigation = function(){
    this.navigation_layer = new Kinetic.Layer();
    var self = this;

    if(self.zoomed_clusters.length>0){
        var refresh_icon = new Kinetic.Path({
            data: "M24.083,15.5c-0.009,4.739-3.844,8.574-8.583,8.583c-4.741-0.009-8.577-3.844-8.585-8.583c0.008-4.741,3.844-8.577,8.585-8.585c1.913,0,3.665,0.629,5.09,1.686l-1.782,1.783l8.429,2.256l-2.26-8.427l-1.89,1.89c-2.072-1.677-4.717-2.688-7.587-2.688C8.826,3.418,3.418,8.826,3.416,15.5C3.418,22.175,8.826,27.583,15.5,27.583S27.583,22.175,27.583,15.5H24.083z",
            x: 50,
            y: 10,
            fill: "grey",
            scale: 1,
            id: "refresh_icon",
            label: "Refresh"
        });
        self.navigation_layer.add(refresh_icon);

        var back_icon = new Kinetic.Path({
            data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z",
            x: 90,
            y: 10,
            fill: "grey",
            scale: 1,
            label: "Zoom out"
        });

        var refresh_overlay = self._draw_icon_overlay(50, 10);
        var back_overlay = self._draw_icon_overlay(90, 10);

        self.navigation_layer.add(back_icon);
        self.navigation_layer.add(refresh_overlay);
        self.navigation_layer.add(back_overlay);

        back_overlay.on("click", function(){
            self._back_icon_click(this);
        });

        back_overlay.on("mouseover", function(){
            self._icon_mouseover(back_icon, back_overlay, self.navigation_layer);
        });

        back_overlay.on("mouseout", function(){
            self._icon_mouseout(back_icon, back_overlay, self.navigation_layer);
        });

        refresh_overlay.on("click", function(){
            self._refresh_icon_click(this);
        });

        refresh_overlay.on("mouseover", function(){
            self._icon_mouseover(refresh_icon, refresh_overlay, self.navigation_layer);
        });

        refresh_overlay.on("mouseout", function(){
            self._icon_mouseout(refresh_icon, refresh_overlay, self.navigation_layer);
        });
    }

    if(self.header){
        var filter_icon = new Kinetic.Path({
                data: "M26.834,6.958c0-2.094-4.852-3.791-10.834-3.791c-5.983,0-10.833,1.697-10.833,3.791c0,0.429,0.213,0.84,0.588,1.224l8.662,15.002v4.899c0,0.414,0.709,0.75,1.583,0.75c0.875,0,1.584-0.336,1.584-0.75v-4.816l8.715-15.093h-0.045C26.625,7.792,26.834,7.384,26.834,6.958zM16,9.75c-6.363,0-9.833-1.845-9.833-2.792S9.637,4.167,16,4.167c6.363,0,9.834,1.844,9.834,2.791S22.363,9.75,16,9.75z",
                x: 10,
                y: 10,
                fill: "grey",
                scale: 1,
                label: "Filter\ncolumns"
        });
        self.navigation_layer.add(filter_icon);

        var filter_overlay = self._draw_icon_overlay(10, 10);
        self.navigation_layer.add(filter_overlay);


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
    self.stage.add(self.navigation_layer);
}

InCHlib.prototype._draw_icon_overlay = function(x, y){
    var icon_overlay = new Kinetic.Rect({
            x: x,
            y: y,
            width: 32,
            height: 32,
            opacity: 0,
        });
    return icon_overlay;
}

InCHlib.prototype._highlight_path = function(path_id){
    var node = this.data.nodes[path_id];
    if(node.count != 1){
        var path = this._clone_path(this.dendrogram_layer.get("#"+path_id)[0]);
        // var path = this.dendrogram_layer.get("#"+path_id)[0].clone();
        
        if(path){
            path.setStroke("#F5273C");
            var rect = this._clone_rect(this.dendrogram_layer.get("#"+path_id+"_rect")[0], path);
            // var rect = this.dendrogram_layer.get("#"+path_id+"_rect")[0].clone();
            rect.setAttrs({path: path, path_id: path_id});
            
            this.dendrogram_overlay_layer.add(path);
            this.dendrogram_overlay_layer.add(rect);
            this._highlight_path(node.left_child);
            this._highlight_path(node.right_child);
        }
    }
    else{
        if(this.settings.heatmap){
            var y_coordinate = this.leaves_y_coordinates[path_id];
            this.row_y_coordinates.push(y_coordinate);
         }
     }
 }

 InCHlib.prototype.highlight_rows = function(row_ids, action){
    var i, row;

    if(action == "append"){
        this.settings.highlighted_rows.push.apply(this.settings.highlighted_rows, row_ids);
    }
    else{
        this.settings.highlighted_rows = row_ids;
        this.highlighted_rows_layer.removeChildren();
    }
    
    var prefixed_row_ids = [];
    for(i = 0; i<row_ids.length; i++){
        prefixed_row_ids.push([this.settings.target, row_ids[i]].join("#"));
    }

    var original_colors = this.settings.heatmap_colors;
    var original_metadata_colors = this.settings.metadata_colors;
    this.settings.heatmap_colors = this.settings.highlight_colors;
    this.settings.metadata_colors = this.settings.highlight_colors;

    for(i = 0; i<prefixed_row_ids.length; i++){
        if(prefixed_row_ids[i] in this.leaves_y_coordinates){
            row = this._draw_heatmap_row(prefixed_row_ids[i], this.distance+this.dendrogram_heatmap_distance, this.leaves_y_coordinates[prefixed_row_ids[i]]);
            this.highlighted_rows_layer.add(row);
        }
    }
    this.highlighted_rows_layer.draw();
    this.heatmap_overlay.moveToTop();

    this.settings.heatmap_colors = original_colors;
    this.settings.metadata_colors = original_metadata_colors;

    var self = this;

    this.highlighted_rows_layer.on("mouseover", function(evt){
        var rows = self.stage.getIntersection({x:evt.layerX, y:evt.layerY});
        self.heatmap_layer.fire("mouseover", evt);
    });

    this.highlighted_rows_layer.on("click", function(evt){
        self.settings.current_row_ids_callback([self.highlighted_row.split("#")[1]]);
    });
    
 }

 InCHlib.prototype._highlight_cluster = function(path_id){

    if(path_id != this.last_highlighted_cluster){
        if(this.last_highlighted_cluster){
            this.cluster_layer.removeChildren();
            this.cluster_layer.draw();
            this.dendrogram_overlay_layer.removeChildren();
        }
        this.row_y_coordinates = [];
        this.last_highlighted_cluster = path_id;
        this._highlight_path(path_id);
        this.dendrogram_overlay_layer.draw();
        this._draw_cluster_layer();
    }
    else{
        this._neutralize_cluster();
    }
}

InCHlib.prototype._neutralize_cluster = function(){
    if(this.last_highlighted_cluster){
        this.cluster_layer.removeChildren();
        this.cluster_layer.draw();
        this.dendrogram_overlay_layer.removeChildren();
        this.dendrogram_overlay_layer.draw();
        this.last_highlighted_cluster = null;
    }
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

InCHlib.prototype._draw_cluster_layer = function(){
    this.cluster_layer = new Kinetic.Layer();
    var row_count = this.data.nodes[this.last_highlighted_cluster].count;
    var width = this.distance+this.visible_features*this.pixels_for_dimension+100;
    var x = this._hack_round(this.distance+this.visible_features*this.pixels_for_dimension+20);
    
    if(this.settings.heatmap){
        var y1 = this._hack_round(this.row_y_coordinates[0]-this.pixels_for_leaf/2);
        var y2 = this._hack_round(this.row_y_coordinates[this.row_y_coordinates.length-1]+0.5*this.pixels_for_leaf);
    }
    else{
        var y = this._hack_round(this.nodes_y_coordinates[this.last_highlighted_cluster]-this.pixels_for_leaf/2-20);
    }
    
    var rows_desc = new Kinetic.Text({
                    x: x,
                    y: y1-20,
                    text: row_count + " rows",
                    fontSize: 14,
                    fontFamily: this.settings.font,
                    fontStyle: 'bold',
                    fill: 'grey',
                    listening: false,
                });

    var zoom_icon = new Kinetic.Path({
                    data: "M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM15.687,9.051h-4v2.833H8.854v4.001h2.833v2.833h4v-2.834h2.832v-3.999h-2.833V9.051z",
                    x: x-5,
                    y: y2 + 8,
                    fill: "grey",
                    scale: 1,
                    label: "Zoom in",
                });

    var upper_border = new Kinetic.Line({
        points: [0, y1, width, y1],
        stroke: 'black',
        strokeWidth: 2,
        lineJoin: 'round',
        dash: [10, 5],
        listening: false,
    });

    var lower_border = new Kinetic.Line({
        points: [0, y2, width, y2],
        stroke: 'black',
        strokeWidth: 2,
        lineJoin: 'round',
        dash: [10, 5],
        listening: false,
    });

    this._collect_row_ids(y1, y2);
    
    var x1 = this.distance+this.dendrogram_heatmap_distance;
    var y1 = this.header_height;
    var width = this.visible_features*this.pixels_for_dimension+90;
    var height = this.settings.height-this.header_height;    
    var upper_y = this.row_y_coordinates[0]-this.pixels_for_leaf/2;
    var lower_y = this.row_y_coordinates[this.row_y_coordinates.length-1]+this.pixels_for_leaf/2;
    
    var upper_overlay = new Kinetic.Rect({
                                x: x1,
                                y: y1,
                                width: width,
                                height: upper_y-y1-1,
                                fill: "white",
                                opacity: 0.6,
                            });

    var lower_overlay = new Kinetic.Rect({
                                x: x1,
                                y: lower_y+1,
                                width: width,
                                height: this.settings.height-lower_y-this.footer_height,
                                fill: "white",
                                opacity: 0.6,
                            });

    this.cluster_layer.add(upper_overlay);
    this.cluster_layer.add(lower_overlay);
    this.cluster_layer.add(upper_border);
    this.cluster_layer.add(lower_border);

    var zoom_overlay = this._draw_icon_overlay(x-5, y2);

    this.cluster_layer.add(rows_desc);
    this.cluster_layer.add(zoom_icon);
    this.cluster_layer.add(zoom_overlay);
    this.stage.add(this.cluster_layer);

    var self = this;

    this.cluster_layer.on("mouseenter", function(evt){
        self._row_mouseout();
    });

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
        self._neutralize_cluster();
    });

}

InCHlib.prototype._zoom_cluster = function(node_id){
    if(node_id != this.zoomed_clusters[this.zoomed_clusters.length-1] && node_id != this.root_id){
        this.zoomed_clusters.push(node_id);
        this._delete_all_layers();
        this._draw_stage_layer();
        this._draw_row_dendrogram(node_id);

        if(this.settings.column_dendrogram && this._visible_features_equal_column_dendrogram_count()){
            this._draw_column_dendrogram(this.column_root_id);
        }

        this._draw_navigation();
        this._draw_heatmap();
        this._draw_heatmap_header();
        this.highlight_rows(this.settings.highlighted_rows);
    }
    else{
        return false;
    }
}

InCHlib.prototype._get_node_neighbourhood = function(node, nodes){
    var node_neighbourhood = {};
    node_neighbourhood.left_node = {};
    node_neighbourhood.right_node = {}
    node_neighbourhood.left_node.left_node = {};
    node_neighbourhood.left_node.right_node = {};
    node_neighbourhood.right_node.left_node = {};
    node_neighbourhood.right_node.right_node = {};

    node_neighbourhood.left_node.left_count = 0.5;
    node_neighbourhood.left_node.right_count = 0.5;
    node_neighbourhood.left_node.right_node.left_count = 0;
    node_neighbourhood.left_node.right_node.right_count = 0;
    node_neighbourhood.left_node.left_node.left_count = 0;
    node_neighbourhood.left_node.left_node.right_count = 0;

    node_neighbourhood.right_node.left_count = 0.5;
    node_neighbourhood.right_node.right_count = 0.5;
    node_neighbourhood.right_node.right_node.left_count = 0;
    node_neighbourhood.right_node.right_node.right_count = 0;
    node_neighbourhood.right_node.left_node.left_count = 0;
    node_neighbourhood.right_node.left_node.right_count = 0;

    node_neighbourhood.left_count = nodes[node.left_child].count;
    node_neighbourhood.right_count = nodes[node.right_child].count;

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
    
    if(node.count != 1){
        var nodes = this.data.column_dendrogram.nodes;
        var node_neighbourhood = this._get_node_neighbourhood(node, nodes);
        var right_child = nodes[node.right_child];
        var left_child = nodes[node.left_child];
        var x1 = this._get_x1(node_neighbourhood, current_left_count, current_right_count);
        var x2 = this._get_x2(node_neighbourhood, current_left_count, current_right_count);
        var y1 = this._hack_round(this.vertical_distance - this.vertical_distance_step*node.distance);
        y1 = (y1 == 0)? 2: y1;
        var y2 = y1;

        if(right_child.count == 1){
            x2 = x2 - this.pixels_for_dimension/2;
        }

        var left_distance = this.vertical_distance - this.vertical_distance_step*nodes[node.left_child].distance;
        var right_distance = this.vertical_distance - this.vertical_distance_step*nodes[node.right_child].distance;

        this.column_dendrogram_layer.add(this._draw_vertical_path(node_id, x1, y1, x2, y2, left_distance, right_distance));
        this._draw_column_dendrogram_node(node.left_child, left_child, current_left_count - node_neighbourhood.left_node.right_count, current_right_count + node_neighbourhood.left_node.right_count, left_distance, y1);
        this._draw_column_dendrogram_node(node.right_child, right_child, current_left_count + node_neighbourhood.right_node.left_count, current_right_count - node_neighbourhood.right_node.left_count, right_distance, y2);
    }

}

InCHlib.prototype._get_y1 = function(node_neighbourhood, current_left_count, current_right_count){
    current_left_count = current_left_count-node_neighbourhood.left_node.right_count-node_neighbourhood.left_node.left_node.right_count;
    var y = (current_left_count+(node_neighbourhood.left_node.left_node.right_count+node_neighbourhood.left_node.right_node.left_count)/2)*this.pixels_for_leaf;
    return y + this.header_height;
}

InCHlib.prototype._get_y2 = function(node_neighbourhood, current_left_count, current_right_count){
    current_left_count = current_left_count+node_neighbourhood.right_node.left_node.left_count;
    var y = (current_left_count+(node_neighbourhood.right_node.left_node.right_count+node_neighbourhood.right_node.right_node.left_count)/2)*this.pixels_for_leaf;
    return y + this.header_height;
}

InCHlib.prototype._get_x1 = function(node_neighbourhood, current_left_count, current_right_count){
    var column_dendrogram_count = this.data.column_dendrogram.nodes[this.column_root_id].count;
    current_left_count = current_left_count-node_neighbourhood.left_node.right_count-node_neighbourhood.left_node.left_node.right_count;
    var x = (current_left_count+(node_neighbourhood.left_node.left_node.right_count+node_neighbourhood.left_node.right_node.left_count)/2)*this.pixels_for_dimension;
    return (this.distance+this.dendrogram_heatmap_distance+column_dendrogram_count*this.pixels_for_dimension)-x;
}

InCHlib.prototype._get_x2 = function(node_neighbourhood, current_left_count, current_right_count){
    var column_dendrogram_count = this.data.column_dendrogram.nodes[this.column_root_id].count;
    current_left_count = current_left_count+node_neighbourhood.right_node.left_node.left_count;
    var x = (current_left_count+(node_neighbourhood.right_node.left_node.right_count+node_neighbourhood.right_node.right_node.left_count)/2)*this.pixels_for_dimension;;
    return (this.distance+this.dendrogram_heatmap_distance+column_dendrogram_count*this.pixels_for_dimension)-x;
}

InCHlib.prototype._draw_vertical_path = function(path_id, x1, y1, x2, y2, left_distance, right_distance){
    var path = new Kinetic.Line({
                points: [x1, left_distance, x1, y1, x2, y2, x2, right_distance],
                stroke: "grey",
                strokeWidth: "2",
                lineCap: 'sqare',
                lineJoin: 'round',
                listening: false
            });

    return path;
}

InCHlib.prototype._draw_horizontal_path = function(path_id, x1, y1, x2, y2, left_distance, right_distance){
    var path_group = new Kinetic.Group({});

    var path = new Kinetic.Line({
                points: [left_distance, y1, x1, y1, x2, y2, right_distance, y2],
                stroke: "grey",
                strokeWidth: "2",
                lineCap: 'sqare',
                lineJoin: 'round',
                id: path_id,
                listening: false
            });

    var path_rect = new Kinetic.Rect({
        x: x1-1,
        y: y1-1,
        width: this.distance - x1,
        height: y2 - y1,
        fill: "white",
        opacity: 0,
        id: [path_id, "rect"].join("_"),
        path: path,
        path_id: path_id,
    })

    path_group.add(path);
    path_group.add(path_rect);

    return path_group;

}

InCHlib.prototype._filter_icon_click = function(filter_button){
    var filter_list = $("#dendrogram_filter_features").text();
    var position = $("#" + this.settings.target).offset();
    var position_top = position.top + 10;
    var position_left = position.left + 50;
    var symbol = "✖";
    var self = this;

    if(filter_list.length > 0){
        $("#dendrogram_filter_features").fadeIn("slow");
        var target_element = $("#" + self.settings.target);
        $("#dendrogram_overlay").css({"width": target_element.outerWidth(true), "height": target_element.outerHeight(true)});
        $("#dendrogram_overlay").fadeIn("slow");
    }
    else{
        filter_list = "";
        
        for(var attr in self.header){
            if(self.features[attr] == 1){
                symbol = "✔";
            }
            if(attr < this.dimensions){
                var text = self.header[attr];
                if(text == ""){
                    text =  parseInt(attr) + 1 + ". column";
                }

                filter_list = filter_list + "<li class='feature_switch' data-num='" + attr + "'><span class='symbol'>" + symbol + "</span>  " + text +"</li>";
            }
        }
        
        $("#" + self.settings.target).append("<div id='dendrogram_filter_features'><ul>" + filter_list + "</ul><hr /><div><span id='cancel_filter_list'>Cancel</span>&nbsp;&nbsp;&nbsp;<span id='update_filter_list'>Update</span></div></div>");
        var filter_features_element = $("#dendrogram_filter_features");
        
        filter_features_element.css({"display":"none",
            "top":position_top+"px",
            "left":position_left+"px",
            "border-radius":"5px",
            "text-align":"center",
            "position":"absolute",
            "background-color":"#ffffff",
            "border":"solid 3px #DEDEDE",
            "padding-top":"5px",
            "padding-left":"15px",
            "padding-bottom":"10px",
            "padding-right":"15px",
            "font-weight":"bold",
            "z-index": 1000
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

        draw_element_overlay();
        filter_features_element.toTo
        filter_features_element.fadeIn("slow");

        function draw_element_overlay(){
            var target_element = $("#" + self.settings.target);
            var outer_width = target_element.outerWidth(true);
            var width = target_element.outerWidth();
            var outer_height = target_element.outerHeight(true);
            var height = target_element.outerHeight();
            var position = $(target_element).offset();
            var overlay = $("<div id='dendrogram_overlay'></div>");

            overlay.css({"background-color": "white", 
                        "width": width, 
                        "height":height,
                        "position": "absolute",
                        "top": position.top-((outer_height-height)/2),
                        "left": position.left-((outer_width-width)/2),
                        "opacity": 0.5});
            
            target_element.append(overlay);
        }

        $(".feature_switch").click(function(){
            var num = parseInt($(this).attr("data-num"));
            var symbol_element = $(this).find("span");
            self.features[num] = -self.features[num];

            if(self.features[num] == 1){
                symbol_element.text("✔");
                $(this).css("color", "green");
                self.visible_features++;
            }
            else{
                symbol_element.text("✖");
                $(this).css("color", "red");
                self.visible_features--;
            }
        });

        $(function(){
           $("#dendrogram_filter_features ul li, #dendrogram_filter_features div span").hover(
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


        $("#cancel_filter_list").click(function(){
            $("#dendrogram_filter_features").fadeOut("slow");
            $("#dendrogram_overlay").fadeOut("slow");
        });

        $("#dendrogram_overlay").click(function(){
            $("#dendrogram_filter_features").fadeOut("slow");
            $("#dendrogram_overlay").fadeOut("slow");
        });

        $("#update_filter_list").click(function(){
            $("#dendrogram_filter_features").fadeOut("slow");
            $("#dendrogram_overlay").fadeOut("slow");

            var node_id = (self.zoomed_clusters.length > 0)?self.zoomed_clusters[self.zoomed_clusters.length-1]:self.root_id;

            var highlighted_cluster = self.last_highlighted_cluster;
            self.last_highlighted_cluster = null;
            
            self._delete_all_layers();
            if(self.settings.dendrogram){
                self._draw_row_dendrogram(node_id);
            }

            if(self._visible_features_equal_column_dendrogram_count() && self.settings.column_dendrogram){
                self._draw_column_dendrogram(self.column_root_id);
            }

            self._draw_navigation();
            self._draw_heatmap();
            self._draw_heatmap_header();

            if(highlighted_cluster != null){
                self._highlight_cluster(highlighted_cluster);
            }
        });
    }
}

InCHlib.prototype._refresh_icon_click = function(){
    var node_id = this.root_id;
    this.zoomed_clusters = [];
    this._delete_all_layers();
    this._draw_stage_layer();
    this._draw_row_dendrogram(node_id);
    if(this.settings.column_dendrogram && this._visible_features_equal_column_dendrogram_count()){
        this._draw_column_dendrogram(this.column_root_id);
    }
    this._draw_navigation();
    this._draw_heatmap();
    this._draw_heatmap_header();
    this.highlight_rows(this.settings.highlighted_rows);
}

InCHlib.prototype._back_icon_click = function(){
    var current_node_id = this.zoomed_clusters[this.zoomed_clusters.length-1];
    this.zoomed_clusters.pop();

    if(this.zoomed_clusters.length >= 1){
        var node_id = this.zoomed_clusters[this.zoomed_clusters.length-1];
        this.zoomed_clusters.pop();
        this._zoom_cluster(node_id);
    }
    else{
        this._refresh_icon_click();
    }
    this._highlight_cluster(current_node_id);
};

InCHlib.prototype._icon_mouseover = function(icon, icon_overlay, layer){
    var label = icon.getAttr("label");
    var x = icon_overlay.getAttr("x");
    var y = icon_overlay.getAttr("y");
    var width = icon_overlay.getWidth();
    var height = icon_overlay.getHeight();

    this.icon_tooltip = new Kinetic.Label({
        x: x,
        y: y+1.2*height,
        listening: false,
    });

    this.icon_tooltip.add(new Kinetic.Tag({
        fill: 'black',
        pointerWidth: 10,
        pointerHeight: 10,
        lineJoin: 'round',
        listening: false,
    }));

    this.icon_tooltip.add(new Kinetic.Text({
        text: label,
        fontFamily: this.settings.font,
        fontSize: 14,
        padding: 5,
        fill: 'white',
        listening: false,
    }));

    layer.add(this.icon_tooltip);
    
    icon.setFill("black");
    layer.draw();
}

InCHlib.prototype._icon_mouseout = function(icon, icon_overlay, layer){
    this.icon_tooltip.destroy();
    icon.setFill("grey");
    layer.draw();
}

InCHlib.prototype._dendrogram_layers_click=function(evt, layer){
    var path_id = evt.targetNode.attrs.path_id;
    layer.fire("mouseout", evt);
    this._highlight_cluster(path_id);
    layer.draw();
}

InCHlib.prototype._dendrogram_layers_mouseout = function(evt){
    var shape_type = evt.targetNode.className;
    if(shape_type == "Rect"){
        var overlay = this.dendrogram_hover_layer.get("#" + [this.settings.target, "path_hover"].join("_"))[0];
        if(overlay){
            this.dendrogram_hover_layer.get("#" + [this.settings.target, "path_hover"].join("_"))[0].destroy();
            this.dendrogram_hover_layer.draw();
        }
    }
}

InCHlib.prototype._dendrogram_layers_mouseover = function(evt){
    var shape_type = evt.targetNode.className;
    
    if(shape_type == "Rect"){
        var path = evt.targetNode.attrs.path;
        var path_overlay = this._clone_path(path);
        // var path_overlay = path.clone();
        path_overlay.setStrokeWidth(4);
        path_overlay.setId([this.settings.target, "path_hover"].join("_"))
        this.dendrogram_hover_layer.add(path_overlay);
        this.dendrogram_hover_layer.draw();
    }
}

InCHlib.prototype._row_mouseout = function(){
    this.heatmap_overlay.destroyChildren();
    this.heatmap_overlay.draw();
    this.highlighted_row = null;
    this.last_highlighted_row = null;
    this.last_col_number = null;
}

InCHlib.prototype._row_mouseover = function(row_group, col_number){
    var i, line;
    this.highlighted_row = row_group.attrs.id;

    for(i = 0; i < row_group.children.length; i++){
        if(row_group.children[i].className == "Line"){
            line = row_group.children[i];
            break;
        }
    }

    var y = line.attrs.points[1]-0.5*this.pixels_for_leaf;
    var x = this.distance+this.dendrogram_heatmap_distance;
    
    var line_ref = new Kinetic.Line({
                 points: [x, y, x+this.visible_features*this.pixels_for_dimension, y],
                 stroke: "black",
                 strokeWidth: 0.5,
                 lineCap: "round",
                 shadowOffset: 1,
                 dash: [4, 2],
                 listening: false,
             });
    this.heatmap_overlay.add(line_ref);

    y = y+this.pixels_for_leaf;
    this.heatmap_overlay.add(line_ref.clone({points:[x, y, x+this.visible_features*this.pixels_for_dimension, y]}));
}

InCHlib.prototype._draw_col_label = function(row_group, row_values, col_number){
    var i, line;
    this.highlighted_row = row_group.attrs.id;

    for(i = 0; i < row_group.children.length; i++){
        if(row_group.children[i].className == "Line"){
            line = row_group.children[i];
            break;
        }
    }

    var y = line.attrs.points[1]-0.5*this.pixels_for_leaf;
    var x = this.distance+this.dendrogram_heatmap_distance;

    var tooltip = new Kinetic.Label({
        x: x+(col_number+0.5)*this.pixels_for_dimension,
        y: y,
        opacity: 1,
        id: "col_label",
        listening: false,
    });

    tooltip.add(new Kinetic.Tag({
        fill: this.settings.label_color,
        pointerDirection: 'down',
        pointerWidth: 10,
        pointerHeight: 10,
        lineJoin: 'round',
        listening: false,
    }));
    
    tooltip.add(new Kinetic.Text({
        text: [this.header[col_number], row_values[col_number]].join(" "),
        fontFamily: this.settings.font,
        fontSize: 12,
        padding: 8,
        fill: 'white',
        fontStyle: "bold",
        listening: false,
    }));

    this.heatmap_overlay.add(tooltip);
    this.heatmap_overlay.moveToTop();
    this.heatmap_overlay.draw();
    return;
}

InCHlib.prototype._visible_features_equal_column_dendrogram_count = function(){
    var i;
    for(i = 0; i<this.data_dimensions; i++){
        if(this.features[i] == -1){
            return false;
        }
    }
    return true;
}

InCHlib.prototype._clone_path = function(path){
    var path_clone = new Kinetic.Line({
                points: path.attrs.points,
                stroke: path.attrs.stroke,
                strokeWidth: path.attrs.strokeWidth,
                lineCap: path.attrs.lineCap,
                lineJoin: path.attrs.lineJoin,
                listening: path.attrs.listening,
    });
    return path_clone;
}

InCHlib.prototype._clone_rect = function(rect){
    var rect_clone = new Kinetic.Rect({
                x: rect.attrs.x,
                y: rect.attrs.y,
                width: rect.attrs.width,
                height: rect.attrs.height,
                opacity: rect.attrs.opacity,
                listening: rect.attrs.listening,
    });
    return rect_clone;

}

InCHlib.prototype._get_color_for_value = function(value, min, max, middle, color_scale){
    var color = this.colors[color_scale];
    var c1 = color["start"];
    var c2 = color["end"];

    if(min == max){
        return 'rgb('+c1.r+','+c1.g+','+c1.b+')';
    }

    if("middle" in color){
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
    
    var calc_color = 'rgb('+r+','+g+','+b+')';
    return calc_color;
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

InCHlib.prototype._collect_row_ids = function(y1,y2){
    var current_node_ids = [];
    this.current_row_ids = [];
    var current_node, i, j;

    for(i in this.leaves_y_coordinates){
        if(this.leaves_y_coordinates[i] > y1 && this.leaves_y_coordinates[i] < y2){
            current_node_ids.push(i);
        }
    }

    for(i = 0; i<current_node_ids.length; i++){
        current_node = this.data.nodes[current_node_ids[i]];
        for(j = 0; j<current_node.objects.length; j++){
            this.current_row_ids.push(current_node.objects[j]);
        }
    }
    
    this.settings.current_row_ids_callback(this.current_row_ids);
}

InCHlib.prototype._hack_size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)){
            size++;
        }
    }
    return size;
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