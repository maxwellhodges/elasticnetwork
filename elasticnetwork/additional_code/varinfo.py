import numpy as np
import scipy


class Varinfo(object):

    def __init__(self, louvain_ensemble):
        self.louvain_ensemble = np.array(louvain_ensemble)
        self.number_of_partitions = []
        self.n = []
        self.nodes = []
        self.ones = []
        self.VI_mat = []
        self.indices = []
        self.VI = True


    def remove_repeats(self):
        self.number_of_partitions = self.louvain_ensemble.shape[0]
        self.n = self.louvain_ensemble.shape[1]
        #VI_mat = np.zeros(number_of_partitions) - can get rid of this line
        VI = 0

        
        # If all the partitions are identifcal, VI = 0 and there is no need to do the
        # rest of the calculations which are computationally expensive

        if np.all( self.louvain_ensemble == np.tile(self.louvain_ensemble[0,:], (self.number_of_partitions, 1))):
            self.VI = 0 # np.zeros((number_of_partitions, number_of_partitions))) - this refers to VI matrix


        # select only partitions that are different, remove partitions that are repeated
        # need to copy this into new array to stop numpy converting tuples to arrays
        # otherwise 'np.unique' flattens result into single array rather than returning
        # unique rows
        louvain_ensemble = [tuple(row) for row in self.louvain_ensemble]
        temp_louvain_ensemble = np.empty((self.number_of_partitions, ), dtype=object)
        
        # enumerate doesn't work here, creates arrays
        for i in range(len(self.louvain_ensemble)):
            temp_louvain_ensemble[i] = louvain_ensemble[i]

        
        (self.louvain_ensemble, self.indices)  = np.unique(temp_louvain_ensemble, return_inverse = True) 
        self.number_of_partitions = len(self.louvain_ensemble)
        

        #VI_tot = 0
        self.nodes = np.linspace(0,self.n-1,self.n)
        self.ones = np.ones(self.n)

        self.VI_mat = np.zeros(((len(self.louvain_ensemble), self.number_of_partitions)))




        # parallelized
    def parallel_partitions(self, i):
        
        VI_mat_row = np.zeros((1, self.number_of_partitions))
        
        
        A_1 = scipy.sparse.csr_matrix((self.ones , (self.louvain_ensemble[i], self.nodes)))
        n_1_all = np.sum(A_1.toarray(), axis=1)                    

        
        for j in range(i):
            A_2 = scipy.sparse.csr_matrix((self.ones,(self.nodes, self.louvain_ensemble[j])))
            n_2_all = np.sum(A_2.toarray(), axis=0)#check
            #n_12_all = np.dot(A_1, A_2)
            n_12_all = A_1 * A_2

            (rows, cols, n_12) = find(n_12_all.toarray())

            n_1 = n_1_all[np.array(rows).astype(int)]
            n_2 = n_2_all[np.array(cols).astype(int)]

            VI = np.sum(n_12*np.log((n_12**2)/(n_1 * n_2))) # all element-wise operations
            VI = -1/(self.n*np.log(self.n))*VI

            VI_mat_row[0][j] = VI
            #VI_tot = VI_tot+VI

        self.VI_mat[i] = VI_mat_row[0]
        print(i)
        return (i, VI_mat_row)

    
    


def find(matrix):

    unweighted_partition_list = list(zip(np.where(matrix > 0))) #gives [row, col] without values

    edge_partition_list = np.where(matrix > 0)
    partition_list_of_indices = list(zip(edge_partition_list[0], edge_partition_list[1])) #i.e. [(0, 1), (1, 0), (1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]

    values = [matrix[element[0]][element[1]] for element in partition_list_of_indices] #extracts weighting values using indices above

    graph = np.array([unweighted_partition_list[0][0], unweighted_partition_list[1][0], values], dtype=np.float64) #bit ugly, should fix

    return graph


# test = varinfo(louvain_ensemble)

# pool = multiprocessing.Pool()
# pool_result  = pool.map(varinfo(louvain_ensemble).parallel_partitions, range(len(louvain_ensemble)))
# pool.close()
# pool.join()
# #result = [pool_result.get() for result in pool_result]
# print(pool_result)

            
        
        # VI_mat_full = np.zeros((number_of_partitions, len(indices)))

        # for i in range(number_of_partitions):
        #     VI_mat_full[i] = VI_mat[i, np.array(indices)]
        
        # VI_mat_full = VI_mat_full[np.array(indices)]

        # VI_mat = VI_mat_full + np.transpose(VI_mat_full)

        # VI = np.mean(ssd.squareform(VI_mat))

        # return VI # decide on VI_mat








# helper function - replicates Matlab's find function.  Takes ndarray as argument.

