import os, shutil, argparse
from multiprocessing import Pool

def get_args():
    parser = argparse.ArgumentParser(description='Rename locus names for fun_pipeline.py output folder')
    parser.add_argument('-i', '--indir', help='input fun_pipeline.py output folder', required=True)
    parser.add_argument('-t', '--threads', help='number of threads', default=1, type=int, required=False)
    parser.add_argument('-l', '--locus', help='locus tag', required=True)
    return parser.parse_args()

def read_file(path):
    with open(path, 'r') as f:
        return f.read()

def write_file(path, contents):
    with open(path, 'w') as f:
        f.write(contents)

def build_queue(path):
    try:
        for child_path in os.listdir(path):
            build_queue(f"{path}/{child_path}")
    except:
        queue.append(path)

def rename(path):
    try:
        print(f"Renaming {path}")
        contents=read_file(path)
        new_contents=contents.replace("FUN_", f"{locus}_")
        write_file(path, new_contents)
    except:
        print(f"could not remain{path}")

def get_predict_folder(path):
    for folder in os.listdir(path):
        if "predict_folder" in folder:
            return f"{path}/{folder}"

def main():
    args = get_args()
    path_predict_folder = get_predict_folder(args.indir)

    global queue
    global locus
    locus = args.locus
    queue = []
    build_queue(path_predict_folder)

    pool=Pool(args.threads)
    pool.map(rename, queue)



if __name__ == '__main__':
    main()