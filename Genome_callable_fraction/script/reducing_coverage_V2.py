import argparse
import gzip

parser = argparse.ArgumentParser(description='Process a gzipped file.')
parser.add_argument('filename', metavar='filename', type=str,
                    help='the name of the input gzipped file')

args = parser.parse_args()

filename = args.filename

output_filename = "/beegfs/data/mbastian/Enard_postbusco/coverage/"+filename[:-3]+"_reduced.gz"
#output_filename = "/home/mbastian/data/Enard_postbusco/coverage/"+filename[:-3]+"_reduced"

def read_blocks(filename):
    with gzip.open(filename, 'rt') as f:
        current_name = ""
        sum_cor = 0
        count_lines = 0
        for line in f:
            if len(line.split())==3:
                name, pos, cor = line.split()
                if name != current_name:
                    if count_lines > 0:
                        mean_cor = sum_cor / count_lines
                        yield current_name, start_pos, end_pos, mean_cor
                    current_name = name
                    sum_cor = 0
                    count_lines = 0
                    start_pos = int(pos)
                sum_cor += float(cor)
                count_lines += 1
                end_pos = pos  # initialisation de la variable end_pos avant le bloc if
                if count_lines == 100:
                    #end_pos = int(pos)
                    mean_cor = sum_cor / count_lines
                    yield current_name, start_pos, end_pos, mean_cor
                    current_name = ""
                    sum_cor = 0
                    count_lines = 0

        if count_lines > 0: #gestion dernier bloc fichier
            end_pos = int(pos)
            mean_cor = sum_cor / count_lines
            yield current_name, start_pos, end_pos, mean_cor

with open(output_filename, 'wt') as f:
    for name, start_pos, end_pos, mean_cor in read_blocks(filename):
        f.write(name + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(mean_cor) + "\n")
        #print(name + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(mean_cor) + "\n")


