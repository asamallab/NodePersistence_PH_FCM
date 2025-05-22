##################################################################################
# Extract the representative cycle of one-dimentional holes     
# This code runs in Jython   
##################################################################################

import os, sys, csv, time 
import edu.stanford.math.plex4.api
from edu.stanford.math.plex4.api import Plex4
from edu.stanford.math.plex4.streams.impl import ExplicitSimplexStream

def extract_1d_cycles(barcode1_intervals, list_gen, outpathlocal, order_udist, max_filtration_value):
    outfile = open(outpathlocal, "w")
    outfile.write('Birth\tDeath\t1DHole\n')
    for n,key in enumerate(list_gen):
        test=str(barcode1_intervals[n]).split(',');
        test[0]=test[0].strip(' [')
        test[1]=test[1].strip(' )')
        if test[1]=='infinity':
            test[1]=str(max_filtration_value);
        line=str(key).replace(' +', ',').replace('-', '');
        birth = format(order_udist[int(eval(test[0]))], ".16f") 
        death = format(order_udist[int(eval(test[1]))], ".16f")
        # print(type(test), test, birth, death)
        outfile.write(birth+'\t'+death+'\t['+line+']\n')	
        #break
    outfile.close()
    return 'Done Jython'

def cteate_complex(distance_matrix,unique_distances_order, max_filtration_value):
    stream = ExplicitSimplexStream()
    num_points = len(distance_matrix)
    for i in range(num_points):
        stream.addVertex(i)
            
    for i in range(num_points):
        for j in range(i + 1, num_points):
            if distance_matrix[i][j] <= max_filtration_value:
                filtration_value = distance_matrix[i][j]
                stream.addElement([i,j], unique_distances_order[filtration_value])
            for k in range(j + 1, num_points):
                max_dist = max(distance_matrix[i][j], distance_matrix[j][k], distance_matrix[i][k])
                if max_dist <= max_filtration_value:
                    filtration_value = max_dist  
                    stream.addElement([i, j, k], unique_distances_order[filtration_value])
    
    stream.finalizeStream()
    return stream


# Parameters
dataset = sys.argv[1]
path_file = sys.argv[2]
max_filtration_value = float(sys.argv[3])
max_dimension = int(sys.argv[4])

files_list = os.listdir(path_file)
#print(dataset, len(files_list))

for fileno in range(len(files_list)):
    t2 = time.time()
    SubID = files_list[fileno].split('.')[0].split('_')[-1]
    distance_matrix_path = path_file + files_list[fileno]
    with open(distance_matrix_path, "r") as file:
        distance_matrix = [[float(value) for value in row] for row in csv.reader(file)]
    unique_distances = sorted(set(distance_matrix[i][j] for i in range(len(distance_matrix)) for j in range(i + 1, len(distance_matrix))).union({0.0}))
    unique_distances_order = {unique_distances[i]: i for i in range(len(unique_distances))}
    order_udist = {i: unique_distances[i] for i in range(len(unique_distances))}
    # print('len of unique_distances:', len(unique_distances))
    stream = cteate_complex(distance_matrix, unique_distances_order, max_filtration_value)
    num_simplices = stream.getSize()
    #print('num_simplices ', num_simplices)
    persistence=Plex4.getModularSimplicialAlgorithm(max_dimension,2);       
    annotated_intervals = persistence.computeAnnotatedIntervals(stream)
    
    barcode0_intervals = list(annotated_intervals.getIntervalsAtDimension(0))  # For H0
    barcode1_intervals = list(annotated_intervals.getIntervalsAtDimension(1))  # For H1
    list_gen=list(annotated_intervals.getGeneratorsAtDimension(1)) 
    
    outpathlocal = '../OutputFiles/PosCorr/'+dataset+'/Holes_1D_Javaplex/Sub_'+SubID+'.txt' 
    cycle_details = extract_1d_cycles(barcode1_intervals, list_gen, outpathlocal, order_udist, max_filtration_value)
    
    N_simp = (len(barcode0_list[fileno][SubID]), len(barcode1_list[fileno][SubID]))
    print(str(fileno)+ '  Done for '+ str(SubID)+ '  len of unique_distances:', len(unique_distances), str(N_simp)+'  Time: ', time.time() -t2)
    #break


print("Done.")