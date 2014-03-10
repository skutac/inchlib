#coding: utf-8
import csv, json, copy, re, argparse

import numpy, scipy, hcluster, fastcluster, sklearn
from sklearn import preprocessing
from scipy import spatial

RAW_LNKAGES = ["ward", "centroid"]
NUMERIC_DISTANCES =  ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "euclidean", "mahalanobis", "minkowski", "seuclidean", "sqeuclidean"]
BINARY_DISTANCES = ["dice","hamming","jaccard","kulsinski","matching","rogerstanimoto","russellrao","sokalmichener","sokalsneath","yule"]
DATA_TYPES = {"numeric": NUMERIC_DISTANCES,
              "binary":BINARY_DISTANCES}

class Dendrogram():
    """Class which handles the generation of cluster heatmap format of clustered data. 
    As an input it takes a Cluster instance with clustered data."""

    def __init__(self, clustering):
        self.cluster_object = clustering
        self.data_type = clustering.data_type
        self.axis = clustering.clustering_axis
        self.clustering = clustering.clustering
        self.tree = hcluster.to_tree(self.clustering)
        self.data = clustering.data
        self.data_names = clustering.data_names
        self.header = clustering.header
        self.dendrogram = False

    def __get_dendrogram__(self, write_data=True):
        root, nodes = hcluster.to_tree(self.clustering, rd=True)
        node_id2node = {}
        dendrogram = {"nodes":{}}

        for node in nodes:
            node_id = "node@{0}".format(node.id)
            if node.count == 1:
                node_id2node[node_id] = {"count":1, "distance":0}

            else:
                node_left_child = "node@{0}".format(node.get_left().id)
                node_right_child = "node@{0}".format(node.get_right().id)
                node_id2node[node_id] = {"count":node.count, "distance":round(node.dist, 3), "left_child": node_left_child, "right_child": node_right_child}

        for n in node_id2node:
            node = node_id2node[n]
            if node["count"] != 1:
                node_id2node[node["left_child"]]["parent"] = n
                node_id2node[node["right_child"]]["parent"] = n

        for n in node_id2node:
            node = node_id2node[n]

            original_id = int(n.split("@")[-1])
            if node["count"] == 1:
                data = self.data[original_id]
                node["objects"] = [self.data_names[original_id]]
                node_id = self.data_names[original_id]

                while node_id in dendrogram["nodes"]:
                    node_id = self.__create_unique_id__(node_id)

                if node_id2node[node["parent"]]["left_child"] == n:
                    node_id2node[node["parent"]]["left_child"] = node_id
                else:
                    node_id2node[node["parent"]]["right_child"] = node_id

                if not write_data:
                    data = []

                node["features"] = data
                dendrogram["nodes"][node_id] = node

        for n in node_id2node:
             if node_id2node[n]["count"] != 1:
                dendrogram["nodes"][n] = node_id2node[n]

        return dendrogram

    def __get_column_dendrogram__(self):
        root, nodes = hcluster.to_tree(self.cluster_object.column_clustering, rd=True)
        node_id2node = {}
        dendrogram = {"nodes":{}}

        for node in nodes:
            node_id = "node@{0}".format(node.id)
            if node.count == 1:
                node_id2node[node_id] = {"count":1, "distance":0}

            else:
                node_left_child = "node@{0}".format(node.get_left().id)
                node_right_child = "node@{0}".format(node.get_right().id)
                node_id2node[node_id] = {"count":node.count, "distance":round(node.dist, 3), "left_child": node_left_child, "right_child": node_right_child}

        for n in node_id2node:
            node = node_id2node[n]
            if node["count"] != 1:
                node_id2node[node["left_child"]]["parent"] = n
                node_id2node[node["right_child"]]["parent"] = n

        for n in node_id2node:
             if not n in dendrogram["nodes"]:
                dendrogram["nodes"][n] = node_id2node[n]

        return dendrogram

    def create_dendrogram(self, contract_clusters=False, cluster_count=1000, write_data=True):
        """Creates cluster heatmap representation in inchlib format. By setting contract_clusters to True you can
        cut the dendrogram in a distance to decrease the row size of the heatmap to count specified by cluster count parameter.
        By setting write_data to False the data features won't be present in the resulting format."""
        self.dendrogram = {"data": self.__get_dendrogram__(write_data)}

        self.contract_clusters = contract_clusters
        self.contract_cluster_treshold = 0
        if self.contract_clusters:
            self.contract_cluster_treshold = self.__get_distance_treshold__(cluster_count)
            print "Distance treshold for contraction:", self.contract_cluster_treshold
            if self.contract_cluster_treshold >= 0:
                self.__contract_data__()

        if self.header and write_data:
            self.dendrogram["data"]["feature_names"] = [h for h in self.header]
        elif self.header and not write_data:
            self.dendrogram["data"]["feature_names"] = []
        
        if self.axis == "both" and len(self.cluster_object.column_clustering):
            column_dendrogram = hcluster.to_tree(self.cluster_object.column_clustering)            
            self.dendrogram["column_dendrogram"] = self.__get_column_dendrogram__()
        return

    def __contract_data__(self):
        nodes = {}
        to_remove = set()
        
        for n in self.dendrogram["data"]["nodes"]:
            node = self.dendrogram["data"]["nodes"][n]

            if node["count"] == 1:
                objects = node["objects"]
                data = node["features"]
                node_id = n

                while self.dendrogram["data"]["nodes"][node["parent"]]["distance"] <= self.contract_cluster_treshold:
                    to_remove.add(node_id)
                    node_id = node["parent"]
                    node = self.dendrogram["data"]["nodes"][node_id]

                if node["count"] != 1:

                    if not "objects" in self.dendrogram["data"]["nodes"][node_id]:
                        self.dendrogram["data"]["nodes"][node_id]["objects"] = []
                        self.dendrogram["data"]["nodes"][node_id]["features"] = []
                    
                    self.dendrogram["data"]["nodes"][node_id]["objects"].extend(objects)

                    if data:
                        self.dendrogram["data"]["nodes"][node_id]["features"].append(data)

        for node in to_remove:
            self.dendrogram["data"]["nodes"].pop(node)

        for k in self.dendrogram["data"]["nodes"]:
            node = self.dendrogram["data"]["nodes"][k]
            if "objects" in node and node["count"] != 1:
                self.dendrogram["data"]["nodes"][k]["distance"] = 0
                self.dendrogram["data"]["nodes"][k]["count"] = 1
                self.dendrogram["data"]["nodes"][k].pop("left_child")
                self.dendrogram["data"]["nodes"][k].pop("right_child")
                rows = zip(*self.dendrogram["data"]["nodes"][k]["features"])
                self.dendrogram["data"]["nodes"][k]["features"] = [round(numpy.median(row), 3) for row in rows]

        self.__adjust_node_counts__()

        return

    def __adjust_node_counts__(self):
        leaves = []

        for n in self.dendrogram["data"]["nodes"]:
            if self.dendrogram["data"]["nodes"][n]["count"] > 1:
                self.dendrogram["data"]["nodes"][n]["count"] = 0
            else:
                leaves.append(n)

        for n in leaves:
            node = self.dendrogram["data"]["nodes"][n]
            parent_id = node["parent"]

            while parent_id:
                node = self.dendrogram["data"]["nodes"][parent_id]
                self.dendrogram["data"]["nodes"][parent_id]["count"] += 1
                parent_id = False
                if "parent" in node:
                    parent_id = node["parent"]
        return

    def __create_unique_id__(self, node_id):
        if re.match(".*?#\d+", node_id):
            node_id, num = node_id.split("#")
            num = str(int(num)+1)
            node_id = "#".join([node_id, num])
        else:
            node_id = "#".join([node_id, "2"])
        return node_id

    def __get_distance_treshold__(self, cluster_count=100):
        print "Calculating distance treshold for cluster contraction..."
        if cluster_count >= self.tree.count:
            return -1
        
        i = 0
        count = cluster_count + 1
        test_step = self.tree.dist/2

        while test_step >= 0.1:
            count = len(set([c for c in hcluster.fcluster(self.clustering, i, "distance")]))
            if count < cluster_count:
                if i == 0:
                    return 0
                i = i - test_step
                test_step = test_step/2
            elif count == cluster_count:
                return i
            else:
                i += test_step

        return i+test_step*2

    def export_dendrogram_as_json(self, filename=None):
        """Returns cluster heatmap in a JSON format or exports it to the file specified by the filename parameter."""
        dendrogram_json = json.dumps(self.dendrogram, indent=4)
        if filename:
            output = open(filename, "w")
            output.write(dendrogram_json)
        return dendrogram_json

    def add_metadata_from_file(self, metadata_file, delimiter, header=True):
        """Adds metadata from csv file."""
        self.metadata, self.metadata_header = self.__read_metadata_file__(metadata_file, delimiter, header)
        self.__connect_metadata_to_data__()
        return

    def add_metadata(self, metadata, header=True):
        """Adds metadata in a form of list of lists (tuples)"""
        self.metadata, self.metadata_header = self.__read_metadata__(metadata, header)
        self.__connect_metadata_to_data__()
        return

    def __connect_metadata_to_data__(self):
        if len(set(self.metadata.keys()) & set(self.data_names)) == 0:
            raise Exception("Metadata objects must correspond with original data objects.")

        if not self.dendrogram:
            raise Exception("You must create dendrogram before adding metadata.")

        self.dendrogram["metadata"] = {"nodes":{}}

        if self.metadata_header:
            self.dendrogram["metadata"]["feature_names"] = self.metadata_header

        leaves = {n:self.dendrogram["data"]["nodes"][n] for n in self.dendrogram["data"]["nodes"] if self.dendrogram["data"]["nodes"][n]["count"] == 1}

        if not self.contract_clusters:
            
            for leaf in leaves:
                try:
                    self.dendrogram["metadata"]["nodes"][leaf] = self.metadata[leaf]
                except Exception, e:
                    continue
        else:

            for leaf in leaves:
                objects = []
                for item in leaves[leaf]["objects"]:
                    try:
                        objects.append(self.metadata[item])
                    except Exception, e:
                        continue

                cols = zip(*objects)
                row = []
                cols = [list(c) for c in cols]

                for col in cols:
                    col.sort()
                    
                    try:
                        value = round(numpy.median(col), 3)
                    except Exception, e:
                        value = col[int(len(col)/2)]
                    row.append(value)

                self.dendrogram["metadata"]["nodes"][leaf] = row

        return

    def __read_metadata__(self, metadata, header):
        metadata_header = []
        rows = metadata
        metadata = {}
        data_start = 0

        if header:
            metadata_header = rows[0][1:]
            data_start = 1
        
        for row in rows[data_start:]:
            metadata[str(row[0])] = [r for r in row[1:]]

        return metadata, metadata_header

        
    def __read_metadata_file__(self, metadata_file, delimiter, header):
        csv_reader = csv.reader(open(metadata_file, "r"), delimiter=delimiter)
        metadata_header = []
        rows = [row for row in csv_reader]
        metadata = {}
        data_start = 0

        if header:
            metadata_header = rows[0][1:]
            data_start = 1
        
        for row in rows[data_start:]:
            metadata[str(row[0])] = [r for r in row[1:]]

        return metadata, metadata_header


class Cluster():
    """Class for clustering"""

    def __init__(self):
        self.write_original = False

    def read_csv(self, filename, delimiter=",", header=False):
        """Reads data from the CSV file"""
        self.filename = filename
        self.delimiter = delimiter
        self.header = header
        csv_reader = csv.reader(open(self.filename, "r"), delimiter=delimiter)
        rows = [row for row in csv_reader]
        self.read_data(rows, header)

    def read_data(self, rows, header=False):
        """Reads data in a form of list of lists (tuples)"""
        self.header = header
        data_start = 0

        if self.header:
            self.header = rows[0][1:]
            data_start = 1
        
        self.data_names = [str(row[0]) for row in rows[data_start:]]
        self.data = [[round(float(value), 3) for value in row[1:]] for row in rows[data_start:]]
        return

    def __binarize_nominal_data__(self):
        nominal_data = zip(*self.data)
        pos2bits = {}
        self.binarized_data = []
        
        for x in xrange(len(nominal_data)):
            unique = list(set(nominal_data[x]))
            unique.sort()
            pos2bits[x] = unique

        for row in self.data:
            binarized_row = []
            for pos in xrange(len(row)):
                binarized_position = [0 for x in xrange(len(pos2bits[pos]))]
                binarized_position[pos2bits[pos].index(row[pos])] = 1
                binarized_row.extend(binarized_position)
            self.binarized_data.append(binarized_row)

        self.original_data = copy.deepcopy(self.data)
        self.data = self.binarized_data
        return

    def __integerize_nominal_data(self__):
        nominal_data = zip(*self.original_data)
        pos2categories = {}
        self.integerized_data = []
        
        for x in xrange(len(nominal_data)):
            unique = list(set(nominal_data[x]))
            unique.sort()
            pos2categories[x] = unique

        for row in self.original_data:
            integerized_row = []
            for pos in xrange(len(row)):
                integerized_row.append(pos2categories[pos].index(row[pos]))
            self.integerized_data.append(integerized_row)
        self.data = self.integerized_data
        return            

    def normalize_data(self, feature_range=(0,1), write_original=False):
        """Normalizes data to a scale from 0 to 1. When write_original is set to True, 
        the normalized data will be clustered, but original data will be written to the heatmap."""
        self.write_original = write_original
        min_max_scaler = preprocessing.MinMaxScaler(feature_range)
        self.original_data = copy.deepcopy(self.data)
        self.data = min_max_scaler.fit_transform(self.data)
        self.data = [[round(v, 3) for v in row] for row in self.data]
        return

    def cluster_data(self, data_type="numeric", row_distance="euclidean", row_linkage="single", axis="row", column_distance="euclidean", column_linkage="ward"):
        """Performs clustering according to the given parameters."""
        print "Cluster analysis setup:", data_type, row_distance, row_linkage, axis
        self.data_type = data_type
        self.clustering_axis = axis
        row_linkage = str(row_linkage)

        if data_type == "nominal":
            self.__binarize_nominal_data__()
        
        if row_linkage in RAW_LNKAGES:
            self.clustering = fastcluster.linkage(self.data, method=row_linkage, metric=row_distance)

        else:
            self.distance_vector = fastcluster.pdist(self.data, row_distance)

            if data_type == "numeric" and not row_distance in DATA_TYPES[data_type]:
                raise Exception("".join(["When clustering numeric data you must choose from these distance measures: ", ", ".join(DATA_TYPES[data_type])]))
            elif (data_type == "binary" or data_type == "nominal") and not row_distance in DATA_TYPES[data_type]:
                raise Exception("".join(["When clustering binary or nominal data you must choose from these distance measures: ", ", ".join(DATA_TYPES[data_type])]))
            elif not data_type in DATA_TYPES.keys():
                raise Exception("".join(["You can choose only from data types: ", ", ".join(DATA_TYPES.keys())]))

            self.clustering = fastcluster.linkage(self.distance_vector, method=str(row_linkage))

        self.column_clustering = []
        if axis == "both" and len(self.data[0]) > 2:
            self.__cluster_columns__(column_distance, column_linkage)

        if data_type == "nominal":
            self.__integerize_nominal_data__()

        if self.write_original:
            self.data = self.original_data

        return

    def __cluster_columns__(self, column_distance, column_linkage):
        columns = zip(*self.data)
        self.column_clustering = fastcluster.linkage(columns, method=column_linkage, metric=column_distance)
        self.data_order = hcluster.leaves_list(self.column_clustering)
        self.__reorder_data__()
        return

    def __reorder_data__(self):
        for i in xrange(len(self.data)):
            reordered_data = []
            for j in self.data_order:
                reordered_data.append(self.data[i][j])
            reordered_data.reverse()
            self.data[i] = reordered_data

        if self.header:
            data = self.header
            reordered_data = []
            for i in self.data_order:
                reordered_data.append(data[i])
            reordered_data.reverse()
            self.header = reordered_data
        return

def _process_(arguments):
    c = Cluster()
    c.read_csv(arguments.data_file, arguments.data_delimiter, arguments.data_header)
    
    if arguments.normalize:
        c.normalize_data(feature_range=(0,1), write_original=arguments.write_original)

    c.cluster_data(data_type=arguments.datatype, row_distance=arguments.row_distance, row_linkage=arguments.row_linkage, axis=arguments.axis, column_distance=arguments.column_distance, column_linkage=arguments.column_linkage)

    d = Dendrogram(c)
    d.create_dendrogram(contract_clusters=arguments.compress, cluster_count=arguments.compress, write_data= not arguments.dont_write_data)
    
    if arguments.metadata:
        d.add_metadata_from_file(metadata_file=arguments.metadata, delimiter=arguments.metadata_delimiter, header=arguments.metadata_header)
    
    if arguments.output_file:
        d.export_dendrogram_as_json(arguments.output_file)
    else:
        print json.dumps(d.dendrogram, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("data_file", type=str, help="csv(text) data file with delimited values")
    parser.add_argument("-o", "--output_file", type=str, help="the name of output file")
    parser.add_argument("-rd", "--row_distance", type=str, default="euclidean", help="set the distance to use for clustering rows")
    parser.add_argument("-rl", "--row_linkage", type=str, default="ward", help="set the linkage to use for clustering rows")
    parser.add_argument("-cd", "--column_distance", type=str, default="euclidean", help="set the distance to use for clustering columns (only when clustering by both axis -a parameter)")
    parser.add_argument("-cl", "--column_linkage", type=str, default="ward", help="set the linkage to use for clustering columns (only when clustering by both axis -a parameter)")
    parser.add_argument("-a", "--axis", type=str, default="row", help="define clustering axis (row/both)")
    parser.add_argument("-dt", "--datatype", type=str, default="numeric", help="specify the type of the data (numeric/binary)")
    parser.add_argument("-dd", "--data_delimiter", type=str, default=",", help="delimiter of values in data file")
    parser.add_argument("-m", "--metadata", type=str, default=None, help="csv(text) metadata file with delimited values")
    parser.add_argument("-md", "--metadata_delimiter", type=str, default=",", help="delimiter of values in metadata file")
    parser.add_argument("-dh", "--data_header", default=False, help="whether the first row of data file is a header", action="store_true")
    parser.add_argument("-mh", "--metadata_header", default=False, help="whether the first row of metadata file is a header", action="store_true")
    parser.add_argument("-c", "--compress", type=int, default=0, help="compress the data to contain maximum of specified count of rows")
    parser.add_argument("-dwd", "--dont_write_data", default=False, help="don't write clustered data to the inchlib data format", action="store_true")
    parser.add_argument("-n", "--normalize", default=False, help="normalize data to [0, 1] range", action="store_true")
    parser.add_argument("-wo", "--write_original", default=False, help="cluster normalized data but write the original ones to the heatmap", action="store_true")
    
    args = parser.parse_args()
    _process_(args)
    
