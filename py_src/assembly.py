from os.path import basename, join, splitext, exists

from slugify import slugify


class Assembly:
    def __init__(self, fname=None, name=None, out_dir=None):
        self.fname = fname
        self.real_coords = None
        self.label = name or slugify(splitext(basename(fname))[0])
        self.name = slugify(splitext(basename(fname))[0])
        self.contig_name = "contig"
        self.chains_fname = join(out_dir, "%s.txt" % (self.name))
        self.bed_fname = join(out_dir, "%s.bed" % (self.name))
        self.sam_fname = join(out_dir, "%s.sam" % (self.name))
