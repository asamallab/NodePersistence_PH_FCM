import os, sys, csv, time 
import pickle as pkl
import jarray 
import edu.stanford.math.plex4.api
from edu.stanford.math.plex4.api import Plex4
# from edu.stanford.math.plex4.metric.impl import ExplicitMetricSpace
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
    
    for i in range(num_points):
        for j in range(i + 1, num_points):
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
RSN = sys.argv[3]
RSNs_details_path = sys.argv[4]
rsn = RSN.replace(' ','')
max_filtration_value = float(sys.argv[5])
max_dimension = int(sys.argv[6])
# num_divisions = int(sys.argv[7])

RSNs_details = []
with open(RSNs_details_path, "r") as file:
    reader = csv.DictReader(file)  # Use DictReader to handle column names
    for row in reader:
        RSNs_details.append(row)

# Create the result dictionary
result_dict = {}
for row in RSNs_details:
    key = row["RSN"]
    node_number = row["Node_number"]
    if key not in result_dict:
        result_dict[key] = []
    result_dict[key].append(node_number)

indices = [int(x) for x in result_dict[RSN]]
print(len(result_dict), rsn, len(indices), indices)
files_list = os.listdir(path_file)
#print(dataset, len(files_list))

dgms_list = list()         # This will store the persistence diagrams
barcode0_list = list()     # This will store H0 persistence intervals
barcode1_list = list()     # This will store H1 persistence intervals
# barcode2_list = list()     # This will store H2 persistence intervals


for fileno in range(len(files_list)):
    t2 = time.time()
    SubID = files_list[fileno].split('.')[0].split('_')[-1]
    distance_matrix_path = path_file + files_list[fileno]
    with open(distance_matrix_path, "r") as file:
        DisMat = [[float(value) for value in row] for row in csv.reader(file)]
    distance_matrix = [[DisMat[k][j] for j in indices] for k in indices]
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
    # barcode2_intervals = list(annotated_intervals.getIntervalsAtDimension(2))  # For H2
    All_barcodes = barcode0_intervals + barcode1_intervals #+ barcode2_intervals
    list_gen=list(annotated_intervals.getGeneratorsAtDimension(1)) 
    
    
    # Add the intervals to their respective lists
    dgms_list.append({SubID: [[order_udist[int(interval.getStart())], order_udist[int(interval.getEnd())] if interval.getEnd()!= None else max_filtration_value] for interval in All_barcodes]})
    barcode0_list.append({SubID: [[order_udist[int(interval.getStart())], order_udist[int(interval.getEnd())] if interval.getEnd()!= None else max_filtration_value] for interval in barcode0_intervals]})
    barcode1_list.append({SubID: [[order_udist[int(interval.getStart())], order_udist[int(interval.getEnd())] if interval.getEnd()!= None else max_filtration_value] for interval in barcode1_intervals]})
    
    outpathlocal = '../OutputFiles_Javaplex/PosCorr/'+dataset+'/Holes_1D_UniqDist/'+rsn+'/Sub_'+SubID+'.txt' 
    directory_name = '../OutputFiles_Javaplex/PosCorr/'+dataset+'/Holes_1D_UniqDist/'+rsn+'/'
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    cycle_details = extract_1d_cycles(barcode1_intervals, list_gen, outpathlocal, order_udist, max_filtration_value)
    
    N_simp = (len(dgms_list[fileno][SubID]), len(barcode0_list[fileno][SubID]), len(barcode1_list[fileno][SubID]))#, len(barcode2_list[fileno][SubID]))
    print(str(fileno)+ '  Done for '+ str(SubID)+ '  len of unique_distances:', len(unique_distances), str(N_simp)+'  Time: ', time.time() -t2)
    # if fileno == 5:
    break


outpath = '../OutputFiles_Javaplex/PosCorr/'+dataset+'/Holes_1D_UniqDist/'+rsn+'/'

file = open(outpath+ rsn+ "_dgms_list.pkl", "wb")
try:
    pkl.dump(dgms_list, file)
finally:
    file.close()

file0 = open(outpath+ rsn+ "_barcode0_list.pkl", "wb")
try:
    pkl.dump(barcode0_list, file0)
finally:
    file0.close()

file1 = open(outpath +rsn+ "_barcode1_list.pkl", "wb")
try:
    pkl.dump(barcode1_list, file1)
finally:
    file1.close()

print("Barcodes saved successfully as .pkl files.")