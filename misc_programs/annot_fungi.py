#classes for annotated fungi objects
import argparse, os

def arguments():
    parser = argparse.ArgumentParser(description="Testing classes")
    parser.add_argument("-i", "--input", help="input dir", required=True)
    return parser.parse_args()

def read(file):
    with open(file, "r") as f:
        return f.read()
class fungi:
    def __init__(self, fasta):
        self._fasta = fasta
        self._annot_dir = None  # Initialize annot_dir

    @property
    def fasta(self):
        """The radius property."""
        return self._fasta

    @fasta.setter
    def fasta(self, path):
        self._fasta = path
        self._annot_dir = None  

    @property
    def annot_dir(self):
        if self._annot_dir is None:  # Calculate annot_dir only if necessary
            self._annot_dir = (self._fasta.replace("tw_fun_assem", "tw_fun_assem_annotation")).replace(".fasta", "_annotation")
        return self._annot_dir
    
    @property
    def taxon(self):
        end=self._fasta.split("/")[-1]
        return ' '.join(end.split("_")[0:2])
    
    @property
    def strain(self):
        end=self._fasta.split("/")[-1]
        return (' '.join(end.split("_")[3:])).replace(".fasta", "")
    
    @property
    def basename(self):
        return self.fasta.split("/")[-1].replace(".fasta", "")
    
    @property
    def predict_dir(self):
        predict_dir=f"{self.annot_dir}/{self.basename}_predict_folder"
        return predict_dir
    
    @property
    def species(self):
        return self.basename.split("_")[1]
    
    @property
    def log_dir(self):
        return f"{self.predict_dir}/logfiles"
    
    
    def log_list(self):
        print(self.log_dir)
        for file in os.listdir(self.log_dir):
            print(file)

    def final_ext(self, ext):
        return f"{self.predict_dir}/{self.species}.{ext}"
    
    def attributes(self):
        for attr_name in dir(self.__class__):
            if isinstance(getattr(self.__class__, attr_name), property):
                print(f"{attr_name} = {getattr(self, attr_name)}")


def main():
    args = arguments()
    fungi_s=fungi(args.input)
    print(fungi_s.final_ext("gff3"))
    print(fungi_s.log_file('busco'))


if __name__ == "__main__":
    main()