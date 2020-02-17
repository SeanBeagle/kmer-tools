# IMPORT FROM PYTHON STANDARD LIBRARY
import sys
import os
import itertools
import threading  #Get Kmers Faster
from collections import Counter

import cv2
import random
import string

from types import SimpleNamespace

# EXTERNAL PACKAGES
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

# IMPORT FROM EXTERNAL MODULE:

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Seq
import numpy as np
from PIL import Image

# BLOOM FILTER:
import mmh3
from bitarray import bitarray, bitdiff
import math
from multiprocessing import Process, Pool, Queue




# CONFIGURE APPLICATION
class Config:
    SQLALCHEMY_DATABASE_URI = "sqlite3:///kmerdb.db"
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    DEBUG = True


app = Flask(__name__)
app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)


DATA = "/Users/beagles/Desktop/Data/ref_genomes"
FILES = {'mav': ['MAV.HOM.104.fasta', 'MAV.HOM.H87.fasta', 'MAV.PTB.K10.fasta']}

mavfiles = [os.path.join(DATA, file) for file in FILES.get('mav')]

class Sequence(db.Model):
    __tablename__ = 'records'

    id = db.Column(db.Integer, primary_key=True)

    # SEQUENCE
    file = db.Column(db.String)
    header = db.Column(db.String)
    index = db.Column(db.String)
    basepairs = db.Column(db.Integer)

    # IMAGES



def random_filename(ext, path=os.getcwd()):
    name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
    return f'{path}/{name}.{ext}'


class KmerIndex:
    alphabets = {'nucleotide': 'ATCG'}

    def __init__(self, k=8, alphabet='nucleotide', from_file=False, **kwargs):
        self.alphabet = self.alphabets.get(alphabet)
        self.k = k
        self._array = None
        self._common = None

        if from_file is False:
            self.kmers = k
            self.file = None
        else:
            self._kmers = kwargs.get('kmers')
            self.file = kwargs.get('file')

        self._kmerlookup = {kmer: i for i, kmer in enumerate(self.kmers)}

    def __eq__(self, kmer):
        return self._kmerlookup.get(kmer)

    def __getitem__(self, item):
        if type(item) == int:
            return self.kmers[item]
        elif type(item) == str:
            return self._kmerlookup.get(item)

    def __repr__(self):
        return f"<{self.__class__.__name__} object: k={self.k}>"

    def __len__(self):
        return len(self.kmers)

    @property
    def kmers(self):
        return self._kmers

    @kmers.setter
    def kmers(self, k):
        """Return a sorted list of all kmers of length(k) from alphabet"""
        kmers = [''.join(kmer) for kmer in set(itertools.combinations(self.alphabet*k, k))]
        kmers.sort()
        self._kmers = kmers

    def save(self, file='kmers.txt'):
        """Print every kmer to file on a newline"""
        with open(file, 'w') as fh:
            for i, kmer in enumerate(self.kmers):
                if i == 0:
                    fh.write(kmer)
                else:
                    fh.write(f"\n{kmer}")
        self.file = file

    @classmethod
    def open(cls, file):
        return cls(from_file=True, file=file, kmers=[kmer.strip() for kmer in open(file)])


class Records:
    def __init__(self):
        self.records = []
        self._array = None
        self._common = None
        self.common_barcode = None


    def add(self, file, index='all'):
        if index == 'all':
            n = len(list(SeqIO.parse(file, "fasta")))
            for i in range(n):
                self.records.append(Kmers(file=file, index=i))
        elif isinstance(index, int):
            self.records.append(Kmers(file=file, index=index))

    def __getitem__(self, item):
        return self.records[item]

    def __len__(self):
        return len(self.records)
    @property
    def name(self):
        return "Sean"
    @property
    def array(self):
        if not isinstance(self._array, np.ndarray):
            self._array = np.array([record.list for record in self.records])
        return self._array

    @property
    def common(self):
        if self._common is None:
            l = []

            for i in range(self.array.shape[1]):
                if 0 not in self.array[:, i]:
                    l.append(i)
                else:
                    l.append(0)

            file = random_filename('png')
            indices = np.indices((256, 256))
            x = np.reshape(indices[0], 256 * 256)
            y = np.reshape(indices[1], 256 * 256)
            coords = list(zip(x, y))
            im = Image.new("RGB", (256, 256), "white")

            for i, freq in enumerate(l):
                if freq > 0:
                    im.putpixel(coords[i], (0, 0, 0))
            im.save(file)
            self.common_barcode = file
            print(f"file:{file}")
            self._common = l

        return self._common

    def compare(self, record1, record2):
        image1 = cv2.imread(self[record1].barcode)
        image2 = cv2.imread(self[record2].barcode)

        difference = cv2.subtract(image1, image2)
        result = not np.any(difference)

        if result is True:
            print("Barcode1 = Barcode2")
        else:
            cv2.imwrite("result.png", difference)
            print("Barcode1 != Barcode2")


class Kmers:
    """ Create dictionary of kmers and frequencies from a sequence """
    def __init__(self, file, k=51, index='all', normalize=True, mode='expand'):
        self._kmers = None
        self._list = None
        self.total = None
        self.unique = None
        self.id = None
        self._barcode = None
        self.file = file
        self.index = index
        self.seq = self.from_fasta(file, index=index)
        self.k = k
        self.normalize = normalize
        self.mode = mode
        self.kmers = self.seq

    def __getitem__(self, item):
        if isinstance(item, str):
            if item in self.kmers:
                return item, self.kmers.get(item)
            elif Seq.reverse_complement(item) in self.kmers:
                return item, self.kmers.get(Seq.reverse_complement(item))
            else:
                return item, 0

        elif isinstance(item, int):
            try:
                return list(self._kmers.items())[item]
            except IndexError:
                return None
        elif isinstance(item, slice):
            l = list(self._kmers.items())[item]
            if len(l) > 0:
                print(f"kmers={len(l)}, proportion={sum([v for k, v in l])}")
                return l
        elif isinstance(item, float):
            n = 0
            output = []
            for kmer in list(self.kmers.items()):
                if kmer[1] + n <= item:
                    output.append(kmer)
                    n += kmer[1]
                else:
                    print(f"kmers={len(output)}, proportion={n}")
                    return output

    def __len__(self):
        return len(self.kmers)

    def __repr__(self):
        return f"<{self.__class__.__name__} object: k={self.k}, n={len(self.kmers)}>"

    def from_fasta(self, file, index):
        if not file:
            return None
        if index == 'all':
            return [record.seq for record in SeqIO.parse(file, "fasta")]
        else:
            return [record.seq for i, record in enumerate(SeqIO.parse(file, "fasta"))
                    if i == index]

            # for i, record in enumerate(SeqIO.parse(file, "fasta")):
            #     if i == index:
            #         self.id = record.id
            #         self.file = file
            #         return record.seq

    @property
    def n(self):
        return self.unique

    @property
    def kmers(self):
        return self._kmers

    @kmers.setter
    def kmers(self, seq):
        """ Returns dictionary of {kmer:frequency}"""
        kmers = {}
        # GATHER UNIQUE KMERS AND FREQUENCY
        for seq in seq:
            for i in range(len(seq) - self.k + 1):
                kmer = str(seq[i:i + self.k])
                kmers.setdefault(kmer, 0)
                kmers[kmer] += 1

        # COLLAPSE REVERSE COMPLEMENTS
        if self.mode in ['collapse', 'expand']:
            collapsed = {}
            for k, v in kmers.items():
                if Seq.reverse_complement(k) in collapsed:
                    collapsed[Seq.reverse_complement(k)] += v
                else:
                    collapsed[k] = v
            kmers = collapsed

        # EXPAND REVERSE COMPLEMENTS
        if self.mode == 'expand':
            expanded = {}
            for k, v in kmers.items():
                expanded[Seq.reverse_complement(k)] = v
                expanded[k] = v
            kmers = expanded

        self.total = sum(kmers.values())
        if self.normalize:
            kmers = {k: v / self.total for k, v in kmers.items()}

        # SORT KMERS
        kmers = {k: v for k, v in sorted(kmers.items(), key=lambda kv: kv[1], reverse=True)}
        #print(f"ADDING {self.id}: total={self.total}, unique={len(kmers)}")

        # SET VALUES
        self._kmers = kmers
        self.unique = (len(kmers))

    @property
    def list(self):
        if not self._list:
            self._list = [self[kmer][1] for kmer in KmerIndex(k=self.k)]
        return self._list

    @property
    def barcode(self):
        if self._barcode == None:
            self.barcode = 1
        return self._barcode

    @barcode.setter
    def barcode(self, value):
        file = random_filename('png')

        indices = np.indices((256, 256))
        x = np.reshape(indices[0], 256 * 256)
        y = np.reshape(indices[1], 256 * 256)
        coords = list(zip(x, y))
        im = Image.new("RGB", (256, 256), "white")

        for i, freq in enumerate(self.list):
            if freq > 0:
                im.putpixel(coords[i], (0, 0, 0))
        im.save(file)
        self._barcode = file


class MultiKmer:
    def __init__(self, file):
        self.k = 8
        self.seq = list(SeqIO.parse(file, "fasta"))[0].seq
        self.kmers = {}

    def get_kmers(self):
        threads = [threading.Thread(target=self.add_kmer, args=(self.seq, i, self.k)) for i in range(len(self.seq))]

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        self.kmers = {k: v for k, v in sorted(self.kmers.items(), key=lambda kv: kv[1], reverse=True)}

        print('done')


    def add_kmer(self, seq, i, k):
            kmer = str(seq[i:i + k])
            self.kmers.setdefault(kmer, 0)
            self.kmers[kmer] += 1


class BloomFilter(object):
    """ Class for Bloom filter, using murmur3 hash function """

    def __init__(self, n=3e7, p=0.05):
        self.p = p      # error
        self.n = n      # number of items
        self._m = None  # size of filter in bits
        self._k = None  # number of hash functions

        self.items = 0

        # Bit array of given size
        self.bitarray = bitarray(self.m)
        self.bitarray.setall(0)

        # initialize all bits as 0
        #self.bitarray.setall(0)

    @property
    def m(self):
        if self._m is None:
            self._m = int(-(self.n * math.log(self.p)) / (math.log(2) ** 2))
        return self._m

    @property
    def k(self):
        if self._k is None:
            self._k = int((self.m / self.n) * math.log(2))
        return self._k

    @property
    def usage(self):
        return self.bitarray.count(True) / self.m

    @property
    def size(self):
        return f"{round(self.bitarray.buffer_info()[1]/1000000, 3)} MB"

    def __getitem__(self, item):
        """Return True possible existence of an item in filter"""
        if isinstance(item, str):
            for digest in self.digest(item):
                return self.bitarray[digest]
        if isinstance(item, Kmers):
            return self[item.kmers]
        if isinstance(item, (list, dict)):
            return [self[kmer] for kmer in item].count(True) / len(item)
        if isinstance(item, self.__class__):
            return (self.bitarray & item.bitarray).count(True) / item.bitarray.count(True)

    # def check_in(i):
    #     k = Kmers('intact_phages.fasta', k=51, mode='expand', index=i)
    #     bloom = BloomFilter(3e7)
    #     bloom + k
    #     score1 = (fmint2.bitarray & bloom.bitarray).count(True) / bloom.bitarray.count(True)
    #     score2 = fmint2[k]
    #     if score1 > 0.3 or score2 > 0.3 or 1 == 1:
    #         scores[i] = {"BloomScore:": score1, "K-Score": score2, "Diff": abs(score1 - score2)}
    #
    # def multithread():
    #     scores = {}
    #     threads = [threading.Thread(target=kmer_list, args=(i,)) for i in range(chunks)]
    #     [thread.start() for thread in threads]
    #     [thread.join() for thread in threads]

    def __add__(self, other):
        """Add a string, list of strings, or dictionary to the filter"""
        if isinstance(other, (list, dict)):

            # def add_list(i, step):
            #     print(f"New Process. [{i}:{i+step}] of {len(other)}")
            #     for item in list(other)[i:i+step]:
            #         for i in self.digest(item):
            #             self.bitarray[i] += 1
            #     print("Process Done")
            #
            # step = math.ceil(len(other)/8)
            # print(f"step={step}")
            # workers = [Process(target=add_list, args=(i, step)) for i in range(0, len(other), step)]
            # for i, worker in enumerate(workers):
            #     print(f"Starting Worker: {i+1} of {len(workers)}")
            #     worker.start()
            # for i, worker in enumerate(workers):
            #     print(f"Joining Worker: {i + 1} of {len(workers)}")
            #     worker.join()
            #
            #     # for n, bf in q.get():
            #     #     self.items += n
            #     #     for i in set(bf):
            #     #         self.bitarray[i] = True

            for item in other:
                for i in self.digest(item):
                    self.bitarray[i] = True
                self.items += 1
        elif isinstance(other, str):
            for i in self.digest(other):
                self.bitarray[i] = True
            self.items += 1
        elif isinstance(other, Kmers):
            self.__add__(other.kmers)
        else:
            print(f"{type(other)} cannot be added to BloomFilter")

    def __repr__(self):
        return f"<{self.__class__.__name__}: m={self.m}, n={self.n}, p={self.p}, k={self.k}, items={self.items}>"

    def __contains__(self, item):
        if isinstance(item, str):
            for i in self.digest(item):
                return self.bitarray[i]
        if isinstance(item, self.__class__):
            return (self.bitarray & item.bitarray).count(True) / item.bitarray.count(True)

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return bitdiff(self.bitarray, other.bitarray) / self.m

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return 1 - bitdiff(self.bitarray, other.bitarray) / self.m

    def add(self, item):
        self.__add__(item)

    def digest(self, item):
        """ return a list of hash values """
        return [mmh3.hash(item, i) % self.m for i in range(self.k)]





# def check_in(i):
#     k = Kmers('intact_phages.fasta', k=51, mode='expand', index=i)
#     bloom = BloomFilter(3e7)
#     bloom + k
#     score1 = (fmint2.bitarray & bloom.bitarray).count(True) / bloom.bitarray.count(True)
#     score2 = fmint2[k]
#     if score1 > 0.3 or score2 > 0.3 or 1 == 1:
#         scores[i] = {"BloomScore:": score1, "K-Score": score2, "Diff": abs(score1 - score2)}
#
# def multithread():
#     scores = {}
#     threads = [threading.Thread(target=check_in, args=(i,)) for i in range(70)]
#     [thread.start() for thread in threads]
#     [thread.join() for thread in threads]
