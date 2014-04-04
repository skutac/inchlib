import inchlib_clust

#instantiate the Cluster object
c = inchlib_clust.Cluster()

# read csv data file with specified delimiter, also specify whether there is a header row
c.read_csv(filename="../source_data/wine.csv", delimiter=",", header=True)
# c.read_data(data, header=bool) use read_data() for list of lists instead of a data file

# normalize data to (0,1) scale, but after clustering write the original data to the heatmap
c.normalize_data(feature_range=(0,1), write_original=True)

# cluster data according to the parameters
c.cluster_data(data_type="numeric", row_distance="euclidean", row_linkage="ward", axis="both", column_distance="euclidean", column_linkage="ward")

# instantiate the Dendrogram class with the Cluster instance as an input
d = inchlib_clust.Dendrogram(c)

# create the cluster heatmap representation and define whether you want to compress the data by defining the maximum number of heatmap rows, the resulted value of compressed (merged) rows and whether you want to write the features
d.create_cluster_heatmap(compress=100, compressed_value="median", write_data=True)

# read metadata file with specified delimiter, also specify whether there is a header row
d.add_metadata_from_file(metadata_file="../source_data/wine_metadata.csv", delimiter=",", header=True, metadata_compressed_value="frequency")

# export the cluster heatmap on the standard output or to the file if filename specified
d.export_cluster_heatmap_as_json("/home/ctibor/Desktop/to_delete.json")