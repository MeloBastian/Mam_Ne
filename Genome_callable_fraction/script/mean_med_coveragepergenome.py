import numpy as np
import time
import gzip
import zlib

#pathfile="/home/mbastian/data/Enard_postbusco/coverage/test_Aco2_reduced.gz"
listsp_file=open("/home/mbastian/data/Enard_postbusco/coverage/list_file","r")
listsp=listsp_file.readlines()


def open_depthfile(specie, error_file):
    try:
        file = gzip.open("/home/mbastian/data/Enard_postbusco/coverage/" + specie, "rt")
    except (gzip.BadGzipFile, EOFError):
        print(f"Le fichier {specie} est corrompu ou incomplet.")
        error_file.write(specie + "\tcorrompu\n")
        return None
    except FileNotFoundError:
        print(f"Le fichier {specie} n'existe pas.")
        return None
    else:
        return file

def read_file(specie, error_file):
    try:
        with gzip.open("/home/mbastian/data/Enard_postbusco/coverage/" + specie, 'rt') as f:
            sumcov = 0
            sites = 0
            for elmt in f:
                try:
                    cov = elmt.split()[3] #compte bien toutes les fenetres, pour la moyenne generale
                    sumcov += float(cov)
                    sites += 1
                except(IndexError):
                    print()
                    print(specie, "index error", elmt)
                    print()
                    error_file.write(specie + "\tIndexError\n")
                    return None
                except(ValueError):
                    print()
                    print(specie, "ValueError", elmt)
                    print()
                    error_file.write(specie + "\tValueError\n")
                    return None
            return sumcov, sites
    except (gzip.BadGzipFile, EOFError, zlib.error):
        error_file.write(specie+"\tcorrompu\n")
        print()
        print(f"Le fichier {specie} est corrompu ou incomplet.")
        print()
        return None

def stat_cov (sumcov, sites):
    mean_cov=float(sumcov/sites)
    #med_cov =np.median(list_coverage)
    return mean_cov

def count_out(borne_min, borne_max, sites, file, bed_mask):
    out_min=0
    out_max=0
    count=0
    elmt = file.readline()
    while elmt != "":
        count+=1
        cov=float(elmt.split()[3])
        if cov<borne_min:
            out_min+=1
            bed_mask.write(elmt)
        if cov>borne_max:
            out_max+=1
            bed_mask.write(elmt)
        elmt=file.readline()
    taux_out_min=out_min/sites
    taux_out_max=out_max/sites

    return (taux_out_min, taux_out_max)

def main ():
    t1_tot=time.perf_counter()
    error_file=open("/home/mbastian/data/Enard_postbusco/coverage/bed_mask/error_files","w")
    summary_results=open("/home/mbastian/data/Enard_postbusco/coverage/bed_mask/summary_results_V2","w")
    summary_results.write("##each line correspond to à mean on 100pb max. the coverage is a mean of a mean.\n")
    summary_results.write("specie\tnbsite(max100pb)\tmean_coverage\trate_out_inf\trate_out_sup\n")
    for elmt in listsp:
        elmt=elmt[:-1]#remove \n
        print(elmt)
        t1 = time.perf_counter()
        results=read_file(elmt, error_file)
        if results is not None:
            (sumcov, nbsites)=results #sumcov = sum covergage for all winbdows, nbsites=nb windows, to do the mean
            mean_cov=stat_cov(sumcov,nbsites)
            borne_mean_min = mean_cov*0.5
            borne_mean_max=mean_cov*2
            file2 = open_depthfile(elmt, error_file)
            bed_mask=open("/home/mbastian/data/Enard_postbusco/coverage/bed_mask/"+elmt[:-21]+"_mean_"+str(mean_cov),"w") #vrai moyenne
            (taux_out_mean_min, taux_out_mean_max)=count_out(borne_mean_min, borne_mean_max, nbsites, file2, bed_mask)

            print(str(nbsites)+" sites")
            print("mean = "+str(mean_cov)+", taux de sites inferieur et supérieur aux bornes :"+ str(taux_out_mean_min)+"  "+ str(taux_out_mean_max))
            summary_results.write(elmt[:-21]+"\t"+str(nbsites)+"\t"+str(mean_cov)+"\t"+ str(taux_out_mean_min)+"\t"+ str(taux_out_mean_max)+"\n")
            t2 = time.perf_counter()
            print('time in seconds ', t2 - t1)
    t2_tot=time.perf_counter()
    print("whole time in minuts ", (t2_tot-t1_tot)/60)

main()