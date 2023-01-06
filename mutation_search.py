import sys
import Bio.SeqIO
from collections import defaultdict
from itertools import product
from typing import Dict, List, OrderedDict, Tuple


class SNPSearcher:
    def __init__(self,
                 genom1_path: str,
                 genom2_path: str,
                 k: int = 30
                ) -> None:
        """
            genom1_path: path to file with reads 
                corresponded with first genom
            genom2_path: path to file with reads 
                corresponded with second genom
            k: kmers' length
        """
        self.genom1_path = genom1_path
        self.genom2_path = genom2_path
        self.k = k

        self.genom1_kmer_frequency_dct = self.generate_kmers_from_reads(
            self.genom1_path, self.k
        )
        self.genom2_kmer_frequency_dct = self.generate_kmers_from_reads(
            self.genom2_path, self.k
        )


    def find_snp(self):
        """
            Entire algorithm for snp searching 
        """
        
        genom1_kmer_frequency_dct_wo_errors = self.throw_away_kmers_with_errors(
            self.genom1_kmer_frequency_dct
        )
        genom2_kmer_frequency_dct_wo_errors = self.throw_away_kmers_with_errors(
            self.genom2_kmer_frequency_dct
        )

        snp_kmers_dct = self.find_kmers_with_snp(
            genom1_kmer_frequency_dct_wo_errors,
            genom2_kmer_frequency_dct_wo_errors
        )

        genom1_snps = self.merge_sequences(snp_kmers_dct['genom1_snp_kmers'])
        genom2_snps = self.merge_sequences(snp_kmers_dct['genom2_snp_kmers'])

        snp_pairs = self.get_snp_pairs(genom1_snps, genom2_snps)

        return snp_pairs


    @staticmethod
    def get_snp_pairs(
            genom1_snps: List[str],
            genom2_snps: List[str],
            threshold: int = 3
        ) -> List[Tuple[str]]:
        """
            Match parts of two genomes together.

            genom1_snps: list of genome1 parts with snps
            genom2_snps: list of genome2 parts with snps
            threshold: size of snp in nucleotides
        """
        def distance(str1: str, str2: str) -> int:
            distance = 0
            for i, j in zip(str1, str2):
                if i != j:
                    distance += 1
            return distance
        
        snp_pairs = []
        for g1, g2 in product(genom1_snps, genom2_snps):
            if len(g1) != len(g2):
                continue
            if distance(g1, g2) <= threshold:
                snp_pairs.append((g1, g2))
        
        return snp_pairs


    @staticmethod
    def throw_away_kmers_with_errors(
            kmer_frequency_dct: Dict[str, int],
            threshold: int = 5
        ) -> Dict[str, int]:
        """
            If frequency of k-mer is lower then threshold,
            we concider it as k-mer with error.
        """
        kmer_frequency_dct_wo_errors = {}
        for kmer, frequency in kmer_frequency_dct.items():
            if frequency > threshold:
                kmer_frequency_dct_wo_errors[kmer] = frequency
        
        return kmer_frequency_dct_wo_errors


    @staticmethod
    def find_kmers_with_snp(
            kmer1_frequency_dct: Dict[str, int],
            kmer2_frequency_dct: Dict[str, int],
        ) -> Dict[str, List[str]]:
        """
            If k-mer with high frequency is present only in one genome,
            then we assume that it has snp.
        """
        genom1_snp_kmers = []
        genom2_snp_kmers = []

        for kmer in kmer1_frequency_dct.keys():
            if kmer not in kmer2_frequency_dct.keys():
                genom1_snp_kmers.append(kmer)
        for kmer in kmer2_frequency_dct.keys():
            if kmer not in kmer1_frequency_dct.keys():
                genom2_snp_kmers.append(kmer)
        return dict(
            genom1_snp_kmers=genom1_snp_kmers, 
            genom2_snp_kmers=genom2_snp_kmers)
        


    @staticmethod
    def merge_sequences(sequences: List[str]):
        """
            Build longer part of genome from k-mers with snp.
        """
        stack = []
        max_kmer = []
        max_kmer_len = 0
        stack = [sequences[0]]
        history = defaultdict(int)
        while len(stack) > 0:
            cur = stack.pop()
            for kmer in sequences:
                overlap_pos = len(cur) - len(kmer) + 1
                if cur[overlap_pos:] == kmer[:-1]:
                    cur += kmer[-1]
                    if history[cur] == 0:
                        stack.append(cur)
                        history[cur] = 1
                    break
                if history[kmer] == 0:
                    stack.append(kmer)
                    history[kmer] = 1
            if len(cur) > max_kmer_len:
                max_kmer = [cur]
                max_kmer_len = len(cur)
            elif len(cur) == max_kmer_len:
                max_kmer.append(cur)
        return set(max_kmer)

    @staticmethod
    def _right_expand_sequence(
            original: str,
            sequences: List[str],
            expended_size: int = 50
        ) -> str:
        expanded = original
        while len(expanded) < expended_size:
            for kmer in sequences:
                overlap_pos = len(expanded) - len(kmer) + 1
                if expanded[overlap_pos:] == kmer[:-1]:
                    expanded += kmer[-1]
                    break
        return expanded


    @staticmethod
    def _left_expand_sequence(
            original: str,
            sequences: List[str],
            expended_size: int = 50
        ) -> str:
        expanded = original
        while len(expanded) < expended_size:
            for kmer in sequences:
                overlap_pos = len(expanded) - len(kmer) + 1
                if expanded[:overlap_pos] == kmer[1:]:
                    expanded = kmer[0] + expanded
                    break

        return expanded

    @staticmethod
    def generate_kmers_from_reads(
            path: str,
            k: int
        ) -> OrderedDict[str, int]:
        """
            Create dict with kmer frequency.

            path: path to file with reads
            k: length of each kmer
        """
        kmer_dct = defaultdict(int)
        for record in Bio.SeqIO.parse(path, "fasta"):
            record = str(record.seq)

            last_idx = len(record) - k + 1
            for idx in range(last_idx):
                kmer = record[idx:idx + k]
                kmer_dct[kmer] += 1

        return kmer_dct


def main(path1, path2):
    snp_searcher = SNPSearcher(path1, path2)
    print(snp_searcher.find_snp())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])