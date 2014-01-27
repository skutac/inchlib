#coding: utf-8

import csv, json, copy, re, argparse

import numpy, scipy, hcluster, fastcluster
from scipy import spatial

RAW_LNKAGES = ["ward", "centroid"]
NONBINARY_DISTANCES =  ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "euclidean", "mahalanobis", "minkowski", "seuclidean", "sqeuclidean"]
BINARY_DISTANCES = ["dice","hamming","jaccard","kulsinski","matching","rogerstanimoto","russellrao","sokalmichener","sokalsneath","yule"]
DATA_TYPES = {"nonbinary": NONBINARY_DISTANCES,
              "binary":BINARY_DISTANCES}

class Dendrogram():

    def __init__(self, clustering, heatmap=True):
        self.cluster_object = clustering
        self.data_type = clustering.data_type
        self.axis = clustering.clustering_axis
        self.clustering = clustering.clustering
        self.tree = hcluster.to_tree(self.clustering)
        self.data = clustering.data
        self.data_names = clustering.data_names
        self.header = clustering.header
        self.heatmap = heatmap
        self.dendrogram = False

    def get_dendrogram(self, write_data=True):
        root, nodes = hcluster.to_tree(self.clustering, rd=True)
        node_id2node = {}
        dendrogram = {"nodes":{}}

        for node in nodes:
            node_id = "node@{0}".format(node.id)
            if node.count == 1:
                node_id2node[node_id] = {"count":1, "distance":0}

            else:
                node_left_id = "node@{0}".format(node.get_left().id)
                node_right_id = "node@{0}".format(node.get_right().id)
                node_id2node[node_id] = {"count":node.count, "distance":round(node.dist, 3), "left_id": node_left_id, "right_id": node_right_id}

        for n in node_id2node:
            node = node_id2node[n]
            if node["count"] != 1:
                node_id2node[node["left_id"]]["parent"] = n
                node_id2node[node["right_id"]]["parent"] = n

        for n in node_id2node:
            node = node_id2node[n]

            original_id = int(n.split("@")[-1])
            if node["count"] == 1:
                data = self.data[original_id]
                node["items"] = [self.data_names[original_id]]
                node_id = self.data_names[original_id]

                while node_id in dendrogram["nodes"]:
                    node_id = self.create_unique_id(node_id)

                if node_id2node[node["parent"]]["left_id"] == n:
                    node_id2node[node["parent"]]["left_id"] = node_id
                else:
                    node_id2node[node["parent"]]["right_id"] = node_id

                if not write_data:
                    data = []

                node["data"] = data
                dendrogram["nodes"][node_id] = node

        for n in node_id2node:
             if node_id2node[n]["count"] != 1:
                dendrogram["nodes"][n] = node_id2node[n]

        return dendrogram

    def get_column_dendrogram(self):
        root, nodes = hcluster.to_tree(self.cluster_object.column_clustering, rd=True)
        node_id2node = {}
        dendrogram = {"nodes":{}}

        for node in nodes:
            node_id = "node@{0}".format(node.id)
            if node.count == 1:
                node_id2node[node_id] = {"count":1, "distance":0}

            else:
                node_left_id = "node@{0}".format(node.get_left().id)
                node_right_id = "node@{0}".format(node.get_right().id)
                node_id2node[node_id] = {"count":node.count, "distance":round(node.dist, 3), "left_id": node_left_id, "right_id": node_right_id}

        for n in node_id2node:
            node = node_id2node[n]
            if node["count"] != 1:
                node_id2node[node["left_id"]]["parent"] = n
                node_id2node[node["right_id"]]["parent"] = n

        for n in node_id2node:
             if not n in dendrogram["nodes"]:
                dendrogram["nodes"][n] = node_id2node[n]

        return dendrogram

    def create_dendrogram(self, contract_clusters=False, cluster_count=1000, write_data=True):
        print "Building dendrogram..."
        self.dendrogram = self.get_dendrogram(write_data)

        self.contract_clusters = contract_clusters
        self.contract_cluster_treshold = 0
        if self.contract_clusters:
            self.contract_cluster_treshold = self.get_distance_treshold(cluster_count)
            print "Distance treshold for contraction:", self.contract_cluster_treshold
            if self.contract_cluster_treshold >= 0:
                self.contract_data()

        if self.heatmap:
            self.dendrogram["heatmap"] = {}
            if self.header and write_data:
                self.dendrogram["heatmap"]["header"] = [h for h in self.header]
            elif self.header and not write_data:
                self.dendrogram["heatmap"]["header"] = []
        
        if self.axis == "both" and len(self.cluster_object.column_clustering):
            column_dendrogram = hcluster.to_tree(self.cluster_object.column_clustering)            
            self.dendrogram["column_dendrogram"] = self.get_column_dendrogram()

        return

    def contract_data(self):
        nodes = {}
        to_remove = set()
        
        for n in self.dendrogram["nodes"]:
            node = self.dendrogram["nodes"][n]

            if node["count"] == 1:
                items = node["items"]
                data = node["data"]
                node_id = n

                while self.dendrogram["nodes"][node["parent"]]["distance"] <= self.contract_cluster_treshold:
                    to_remove.add(node_id)
                    node_id = node["parent"]
                    node = self.dendrogram["nodes"][node_id]

                if node["count"] != 1:

                    if not "items" in self.dendrogram["nodes"][node_id]:
                        self.dendrogram["nodes"][node_id]["items"] = []
                        self.dendrogram["nodes"][node_id]["data"] = []
                    
                    self.dendrogram["nodes"][node_id]["items"].extend(items)

                    if data:
                        self.dendrogram["nodes"][node_id]["data"].append(data)

        for node in to_remove:
            self.dendrogram["nodes"].pop(node)

        for k in self.dendrogram["nodes"]:
            node = self.dendrogram["nodes"][k]
            if "items" in node and node["count"] != 1:
                self.dendrogram["nodes"][k]["distance"] = 0
                self.dendrogram["nodes"][k]["count"] = 1
                self.dendrogram["nodes"][k].pop("left_id")
                self.dendrogram["nodes"][k].pop("right_id")
                rows = zip(*self.dendrogram["nodes"][k]["data"])
                self.dendrogram["nodes"][k]["data"] = [round(numpy.median(row), 3) for row in rows]

        self.adjust_node_counts()

        return

    def adjust_node_counts(self):
        leaves = []

        for n in self.dendrogram["nodes"]:
            if self.dendrogram["nodes"][n]["count"] > 1:
                self.dendrogram["nodes"][n]["count"] = 0
            else:
                leaves.append(n)

        for n in leaves:
            node = self.dendrogram["nodes"][n]
            parent_id = node["parent"]

            while parent_id:
                node = self.dendrogram["nodes"][parent_id]
                self.dendrogram["nodes"][parent_id]["count"] += 1
                parent_id = False
                if "parent" in node:
                    parent_id = node["parent"]
        return

    def create_unique_id(self, node_id):
        if re.match(".*?#\d+", node_id):
            node_id, num = node_id.split("#")
            num = str(int(num)+1)
            node_id = "#".join([node_id, num])
        else:
            node_id = "#".join([node_id, "2"])
        return node_id

    def get_distance_treshold(self, cluster_count=100):
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

    def add_heatmap_settings(self):
        self.dendrogram["heatmap"] = {}
        if self.header:
            self.dendrogram["heatmap"]["header"] = [h for h in self.header]

    def export_dendrogram_as_json(self, filename=None):
        dendrogram_json = json.dumps(self.dendrogram, indent=4)
        if filename:
            output = open(filename, "w")
            output.write(dendrogram_json)
        return dendrogram_json

    def add_metadata_from_file(self, metadata_file, delimiter, header=True):
        self.metadata, self.metadata_header = self.read_metadata_file(metadata_file, delimiter, header)
        self.connect_metadata_to_data()
        return

    def add_metadata(self, metadata, header=True):
        self.metadata, self.metadata_header = self.read_metadata(metadata, header)
        self.connect_metadata_to_data()
        return

    def connect_metadata_to_data(self):
        if len(set(self.metadata.keys()) & set(self.data_names)) == 0:
            raise Exception("Metadata IDs must correspond with original data IDs.")

        if not self.dendrogram:
            raise Exception("You must create dendrogram before adding metadata.")

        self.dendrogram["metadata"] = {"nodes":{}}

        if self.metadata_header:
            self.dendrogram["metadata"]["header"] = self.metadata_header

        leaves = {n:self.dendrogram["nodes"][n] for n in self.dendrogram["nodes"] if self.dendrogram["nodes"][n]["count"] == 1}

        if not self.contract_clusters:
            
            for leaf in leaves:
                try:
                    self.dendrogram["metadata"]["nodes"][leaf] = self.metadata[leaf]
                except Exception, e:
                    continue
        else:

            for leaf in leaves:
                items = []
                for item in leaves[leaf]["items"]:
                    try:
                        items.append(self.metadata[item])
                    except Exception, e:
                        continue

                cols = zip(*items)
                self.dendrogram["metadata"]["nodes"][leaf] = {}
                self.dendrogram["metadata"]["nodes"][leaf] = [round(numpy.median(col), 3) for col in cols]

        return

    def read_metadata(self, metadata, header):
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

        
    def read_metadata_file(self, metadata_file, delimiter, header):
        csv_reader = csv.reader(open(metadata_file, "r"), delimiter=delimiter)
        metadata_header = []
        rows = [row for row in csv_reader]
        metadata = {}
        data_start = 0

        if header:
            metadata_header = rows[0][1:]
            data_start = 1
        
        for row in rows[data_start:]:
            metadata[str(row[0])] = [int(r) for r in row[1:]]

        return metadata, metadata_header


class Cluster():

    def __init__(self):
        pass

    def read_csv(self, filename, delimiter=",", header=False):
        self.filename = filename
        self.delimiter = delimiter
        self.header = header
        csv_reader = csv.reader(open(self.filename, "r"), delimiter=delimiter)
        rows = [row for row in csv_reader]
        self.read_data(rows, header)

    def read_data(self, rows, header=False):
        self.header = header
        data_start = 0

        if self.header:
            self.header = rows[0][1:]
            data_start = 1
        
        self.data_names = [str(row[0]) for row in rows[data_start:]]
        self.data = [[float(value) for value in row[1:]] for row in rows[data_start:]]
        return

    def binarize_nominal_data(self):
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

    def integerize_nominal_data(self):
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

            

    def cluster_data(self, data_type="nonbinary", distance_measure="euclidean", linkage="single", axis="row"):
        print "Cluster analysis setup:", data_type, distance_measure, linkage, axis
        self.data_type = data_type
        self.clustering_axis = axis
        linkage = str(linkage)

        if data_type == "nominal":
            self.binarize_nominal_data()
        
        if linkage in RAW_LNKAGES:
            self.clustering = fastcluster.linkage(self.data, method=linkage, metric=distance_measure)

        else:
            self.distance_vector = fastcluster.pdist(self.data, distance_measure)

            if data_type == "nonbinary" and not distance_measure in DATA_TYPES[data_type]:
                raise Exception("".join(["When clustering nonbinary data you must choose from these distance measures: ", ", ".join(DATA_TYPES[data_type])]))
            elif (data_type == "binary" or data_type == "nominal") and not distance_measure in DATA_TYPES[data_type]:
                raise Exception("".join(["When clustering binary or nominal data you must choose from these distance measures: ", ", ".join(DATA_TYPES[data_type])]))
            elif not data_type in DATA_TYPES.keys():
                raise Exception("".join(["You can choose only from data types: ", ", ".join(DATA_TYPES.keys())]))

            self.clustering = fastcluster.linkage(self.distance_vector, method=str(linkage))

        self.column_clustering = []
        if axis == "both" and len(self.data[0]) > 2:
            self.cluster_columns()

        if data_type == "nominal":
            self.integerize_nominal_data()

        return

    def remove_null_rows_from_data(self):
        print "Number of rows:", len(self.data)
        self.data = [row for row in self.data if sum(row) > 0]
        print "After removal of null rows:", len(self.data)
        return

    def remove_zero_variance_columns(self):
        columns = zip(*self.data)
        print "Number of columns:", len(columns)
        columns = [c for c in columns if len(set(c)) > 1]
        print "After removal of zero variance:", len(columns)
        self.data = zip(*columns)
        return

    def cluster_columns(self):
        columns = zip(*self.data)
        self.column_clustering = fastcluster.linkage(columns, method="ward", metric="euclidean")
        self.data_order = hcluster.leaves_list(self.column_clustering)
        self.reorder_data()
        return

    def reorder_data(self):
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

def process(arguments):
    c = Cluster()
    c.read_csv(arguments.file, arguments.delimiter, arguments.header)
    print "ok"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str, help="csv(text) file with delimited values")
    parser.add_argument("distance", type=str, help="set the distance to use for clustering")
    parser.add_argument("linkage", type=str, help="set the linkage to use for clustering")
    parser.add_argument("-d", "--delimiter", type=str, help="delimiter of values in datafile")
    parser.add_argument("--header", help="whether the first row of datafile is a header", action="store_true")
    

    args = parser.parse_args()
    process(args)
    
    # c = Cluster()
    # # c.read_csv(filename="../static/data/target_correlation_matrix.csv", delimiter=",", header=False)
    # # c.cluster_data(data_type="binary", distance_measure="jaccard", linkage="ward", axis="both")

    # # d = Dendrogram(c, heatmap=True)
    # # d.create_dendrogram(contract_clusters=False, cluster_count=1000, write_data=True)
    # # # d.add_metadata(metadata_file="../static/data/activity_imprints_agonist.csv", delimiter=",", header=True)
    # # d.export_dendrogram_as_json("../static/dendrograms/target_correlation_matrix.json")

    # # c = Cluster()
    # c.read_csv(filename="/home/ctibor/Desktop/steroid_matrices/chemogenomic_matrix_b_score.csv", delimiter=",", header=True)
    # # c.remove_null_rows_from_data()
    # # c.remove_zero_variance_columns()
    # c.cluster_data(data_type="nonbinary", distance_measure="euclidean", linkage="ward", axis="both")
    
    # d = Dendrogram(c, heatmap=True)
    # d.create_dendrogram(contract_clusters=False, cluster_count=20, write_data=True)
    # # metadata = [["id","metadata"],[1, "positive"],[2,"negative"]]
    # # d.add_metadata(metadata, header=True)
    # d.export_dendrogram_as_json("/home/ctibor/Desktop/steroid_matrices/bscore_steroid_matrix.csv")
