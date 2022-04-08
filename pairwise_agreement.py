import os, glob, shutil
import argparse
import pandas as pd
import wide2long
import itertools
import subprocess as sp
import time

parser = argparse.ArgumentParser(description='check if the run secceeded')

parser.add_argument('-a', '--groupA', type=str, default = 'time01', 
                    help='first tract groups for pairwise agreement, \
                    possible values: time01 or com01/02/03')
parser.add_argument('-b', '--groupB', type=str, default = 'time02', 
                    help='second tract groups for pairwise agreement, \
                    possible values: time02, or com01/02/03')
parser.add_argument('-n', '--analysis', type=str, default = 'AL_07')
                    

pj='ThaTract'
ct='rtp-pipeline_4.4.1'
def checkResults(ct, a, b, n, pj):
 
    appended_data = []

    baseDir = '/scratch/lmx/ThaTract/Nifti/derivatives'
    codeDir = '/dipc/lmx/GIT/ThaTract'
 
    subseslist=os.path.join(codeDir,"subSesList.txt")
    tractlist = os.path.join("/dipc/lmx", "tmp_tractparams.csv")
    # READ SUBJECTID FILE
    dt = pd.read_csv(subseslist, sep=",", header=0)
    # READ THALAMIC TRACT NAME
    thalist = pd.read_csv(tractlist, sep=",", header=0)
    slabel = thalist.slabel
    if "time" in a and "time" in b:
        Ases = a[-2:]; Bses = b[-2:]
        analysis = n
        if not os.path.exists(f'{baseDir}/{ct}/analysis-{analysis}/pairwise_agreement'):
            os.makedirs(f'{baseDir}/{ct}/analysis-{analysis}/pairwise_agreement')
        for sub, tract in itertools.product(dt["sub"].unique(), slabel):
            T01 =  (f'{baseDir}/{ct}/analysis-{analysis}/sub-{sub}/' +
                    f'ses-T{Ases}/output/{tract}_clean.tck')
            T02 =  (f'{baseDir}/{ct}/analysis-{analysis}/sub-{sub}/' +
                    f'ses-T{Bses}/output/{tract}_clean.tck')
            ref =  (f'{baseDir}/{ct}/analysis-{analysis}/sub-{sub}/' +
                    f'ses-T{Ases}/output/wmMask.nii.gz')
            json = (f'{baseDir}/{ct}/analysis-{analysis}/' + 
                    f'pairwise_agreement/{tract}_time_{sub}.json')
            print(T01)
            if os.path.exists(T01) and os.path.exists(T02):
                print(f"{sub} {tract}")
                start_time = time.time()
                cmdstr = (f"python3 /scratch/lmx/GIT/scilpy/scripts/scil_evaluate_bundles_pairwise_agreement_measures.py " +   
                         f"--reference {ref} {T01} {T02} {json} -f ")
                sp.call(cmdstr, shell=True)
                print(time.time()-start_time)
            else:
                continue
        os.chdir(f'{baseDir}/{ct}/analysis-{analysis}/pairwise_agreement')
        print("merging all jsons....")
        cmdstr = (f" python3 /scratch/lmx/GIT/scilpy/scripts/scil_merge_json.py --keep_separate -f *_time_*.json {baseDir}/{ct}/analysis-{analysis}/all.json")
        sp.call(cmdstr, shell=True)
        alljs = pd.read_json(f"{baseDir}/{ct}/analysis-{analysis}/all.json")
        alljs = alljs.T
        for col in alljs.columns:
            alljs.loc[:, col] = alljs[col].map(lambda x: x if pd.isnull(x) else x[0])
        alljs['tract'], alljs['SUBID'] = alljs.index.str.split('_time_', 1).str
        alljs.to_csv(f"{baseDir}/{ct}/analysis-{analysis}/pairwise_agreement.csv", index=False)
    elif "com" in a and "com" in b:
        Aana = a[3:]; Bana = b[3:]
        json_dir = f'{baseDir}/{ct}/pairwise_agreement/{Aana}_{Bana}'
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        # sleep after submitting 100
        n = 0
        for sub, tract in itertools.product(dt["sub"].unique(), slabel):
            T01 =  (f'{baseDir}/{ct}/analysis-{Aana}/sub-{sub}/' +
                    f'ses-T01/output/{tract}_clean.tck')
            T02 =  (f'{baseDir}/{ct}/analysis-{Bana}/sub-{sub}/' +
                    f'ses-T01/output/{tract}_clean.tck')
            ref =  (f'{baseDir}/{ct}/analysis-{Aana}/sub-{sub}/' +
                    f'ses-T01/output/wmMask.nii.gz')
            json = (f'{baseDir}/{ct}/pairwise_agreement/' + 
                    f'{Aana}_{Bana}/{tract}_sep_{sub}_sep_{a}vs{b}.json')

            if os.path.exists(T01) and os.path.exists(T02):
                start_time = time.time()
                print(f"{sub} {tract}")
                cmdstr = (f" python3 /scratch/lmx/GIT/scilpy/scripts/scil_evaluate_bundles_pairwise_agreement_measures.py " + 
                             f"--reference {ref} {T01} {T02} {json} -f --keep_tmp ")
                #print(cmdstr)
                sp.call(cmdstr, shell=True)
                print(time.time()-start_time)
            else:
                continue
            n =+ 1
            # if  n % 100 == 0: sleep(20)
        #while not os.path.exists(json):
        #    pass
        os.chdir(f'{baseDir}/{ct}/pairwise_agreement/{Aana}_{Bana}')
        print("merging all jsons....")
        cmdstr = (f" python3 /scratch/lmx/GIT/scilpy/scripts/scil_merge_json.py --keep_separate -f *{a}vs{b}.json " +
              f"{baseDir}/{ct}/pairwise_agreement/{Aana}_{Bana}/all_{a}vs{b}.json")
        sp.call(cmdstr, shell=True)
        alljs = pd.read_json(f"{baseDir}/{ct}/pairwise_agreement/{Aana}_{Bana}/all_{a}vs{b}.json")
        alljs = alljs.T
        for col in alljs.columns:
            alljs.loc[:, col] = alljs[col].map(lambda x: x[0] if isinstance(x, list) else 'NaN')
        alljs['tract'], alljs['SUBID'], alljs['btw'] = alljs.index.str.split('_sep_', 2).str
        alljs.to_csv(f"{baseDir}/{ct}/pairwise_agreement/pairwise_agreement_{a}vs{b}.csv", index=False)


    
if __name__ == "__main__":
    args = parser.parse_args()
    checkResults( ct, args.groupA, args.groupB, args.analysis,  pj)


