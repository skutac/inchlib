import unittest
import random
import inchlib_clust

class Tests(unittest.TestCase):

	def set_cluster(self):
		self.cluster = inchlib_clust.Cluster()

	def test_read_csv(self):
		self.set_cluster()
		self.cluster.read_csv("tests/data/test_data.csv", header=True)
		self.assertEqual(len(self.cluster.data), 49)
		self.assertEqual(len(self.cluster.header), 3)

	def test_compression(self):
		self.set_cluster()
		self.cluster.read_csv("tests/data/test_data.csv", header=True)
		self.cluster.cluster_data(data_type="numeric", row_distance="euclidean", row_linkage="single", axis="both", column_distance="euclidean", column_linkage="ward")
		dendrogram = inchlib_clust.Dendrogram(self.cluster)
		compress = random.randrange(2, len(self.cluster.data) - 1)
		dendrogram.create_cluster_heatmap(compress=compress, compressed_value="median", write_data=True)
		root = [values for n, values in dendrogram.dendrogram["data"]["nodes"].items() if not "parent" in values][0]
		self.assertTrue(root["count"] <= compress)

	def test_normalization(self):
		self.set_cluster()
		self.cluster.read_csv("tests/data/test_data.csv", header=True)
		original_len = len(self.cluster.data)
		self.cluster.normalize_data()
		self.assertEqual(original_len, len(self.cluster.data))

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
	unittest.TextTestRunner(verbosity=2).run(suite)

