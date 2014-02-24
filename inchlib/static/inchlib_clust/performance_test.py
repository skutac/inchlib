import time

import numpy

import inchlib_clust

def generate_data(rows, features):
    data = []
    for i in xrange(rows):
        row = [i]
        row.extend(numpy.random.randint(0, 1000, features))
        data.append(row)

    return data

def test_rows():
    counts = [10, 100, 500, 1000, 2000, 5000, 10000, 20000, 30000]
    times = []

    for c in counts:
        data = generate_data(c, 20)
        cluster_time = get_cluster_time(data, "_".join(["test", str(c)]))
        times.append(round(cluster_time, 3))

    return times


def test_features():
    counts = [10, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 8000, 10000, 20000]
    times = []

    for c in counts:
        data = generate_data(1000, c)
        cluster_time = get_cluster_time(data)
        times.append(round(cluster_time, 3))

    return times

def get_cluster_time(data, filename):
    start = time.time()
    c = inchlib_clust.Cluster()
    c.read_data(data, header=False)
    c.cluster_data(data_type="numeric", distance_measure="euclidean", linkage="ward", axis="rows")

    d = inchlib_clust.Dendrogram(c)
    d.create_dendrogram(contract_clusters=False)
    end = time.time()
    if filename:
        d.export_dendrogram_as_json("".join(["../data/", filename]))
    return end-start

if __name__ == '__main__':
    times = test_rows()
    # times = test_features()
    print times